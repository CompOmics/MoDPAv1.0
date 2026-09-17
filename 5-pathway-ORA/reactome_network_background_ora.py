"""
Supplementary re-analysis: Reactome ORA re-tested against a network-restricted background.

Context
-------
Reactome's Analysis Service always tests against the full annotated human reference
proteome. This was verified directly against the live v3 REST API spec
(https://reactome.org/AnalysisService/v3/api-docs) and the reactome2py source: no
"background" parameter or endpoint exists anywhere in /identifiers or /token. A custom
background is therefore not an option the main pipeline (20260901-reactome-enrichment.ipynb)
failed to use -- it is not offered by the tool at all.

This script does NOT change the manuscript's reported results (the whole-proteome
background is Reactome's standard behaviour and is what is reported). It is a sensitivity
analysis: it answers, with numbers, how much the whole-proteome background matters.

Method
------
1. Submit the proteins present in the MoDPA association network (union of UniAcc across
   all clusters, any size; n=1800 for the 2026-09-13 keen_raman network) as a single
   Reactome ORA query. For every pathway this returns entities.found = how many of those
   network proteins are annotated to it. Call this K, the pathway's size restricted to
   the network background.
2. For EVERY pathway Reactome tested for a given cluster (not just the ones that were
   significant under the whole-proteome background), recompute a one-tailed
   hypergeometric test with:
     N = len(background_proteins)  network background size
     K = network-annotated count from step 1
     n = cluster size (n_proteins) number of cluster proteins submitted to Reactome
     k = entities_found            cluster proteins found in the pathway (from the
                                    per-cluster Reactome result)
   p = P(X >= k) = hypergeom.sf(k - 1, N, K, n)
3. Benjamini-Hochberg-correct these p-values within each cluster, over ALL pathways
   tested for that cluster (the true per-cluster candidate set), not just the ones that
   were significant under the whole-proteome background.

This requires 20260901-reactome-enrichment.ipynb's run_ora() to query Reactome with
p_value="1" (return every tested pathway) and to filter to FDR<FDR_CUTOFF only
afterwards, before the redundancy step -- otherwise the per-cluster cache
(per-cluster/clusterN-reac.csv) only contains the pathways that were already
significant, and step 3 above cannot be done correctly (see "Version history" below).

Reproducing this analysis
-------------------------
Only OUT_DIR below needs to change when the pipeline is rerun on a new clustering.

The single Reactome API call is cached in OUT_DIR/network-background-pathway-sizes.json
(37 KB). If that file is present, this script runs fully offline and deterministically,
with no network access and no reactome2py dependency. Ship it alongside the results so
the sensitivity analysis can be reproduced without depending on a live service whose
release may have moved on. Delete it to force a refresh, e.g. after a Reactome update.

On a cache miss the script needs to reach https://reactome.org. It prefers the
reactome2py client and falls back to a direct urllib call against the same REST
endpoints if reactome2py is not installed. The urllib path sets an explicit User-Agent
because the Analysis Service returns HTTP 403 to clients that do not send one. Both
paths issue exactly one identifiers + token request pair and produce the same K values.

Version history / two bugs this fixes
-------------------------------------
1. An earlier version BH-corrected only over the pathways that were already significant
   under the whole-proteome background (because that was all the old per-cluster cache
   contained). That is WRONG, not merely incomplete: restricting the correction to a
   pre-filtered, small-p-value subset uses a denominator (m) far smaller than the true
   number of pathways tested per cluster, and BH's adjusted p is proportional to m, so
   the old approach was ANTI-CONSERVATIVE. It systematically under-stated adjusted
   p-values, making the network background look more forgiving than it really is.
   Confirmed with a worked example (5-pathway restricted subset vs the true 20-pathway
   universe: adjusted p for the same raw p-value rose from 0.005 to 0.02, a 4x change).
   Fixed by reading the full per-cluster cache, so BH uses the true per-cluster m.
2. The script was pinned to the superseded output directory reactome-128_ref-0_6 and
   imported reactome2py unconditionally, so it did not run against the current results.
   Fixed by the OUT_DIR constant below and the urllib fallback described above.

Inputs (all read-only, produced by 20260901-reactome-enrichment.ipynb)
------------------------------------------------------------------
OUT_DIR/nicely-formatted-nodes.csv           cluster membership, PTM events
OUT_DIR/cluster-sizes.csv                    n_proteins per cluster
OUT_DIR/per-cluster/cluster<N>-reac.csv      ALL pathways Reactome tested for cluster N
                                              (p_value="1"), one file per tested cluster
                                              (>=20 proteins)

Outputs
-------
OUT_DIR/network-background-pathway-sizes.json   cache: {stId: K}, one Reactome API call
OUT_DIR/network-background-reanalysis.csv       per pathway-cluster pair, for EVERY
    pathway tested (not just originally-significant ones): whole-proteome FDR,
    network-restricted raw p and properly-corrected FDR, side by side
"""
import os
import time
import json
import urllib.parse
import urllib.request

import pandas as pd
import numpy as np
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "2026-09-13-keen_raman-Reactome")
PER_CLUSTER_DIR = os.path.join(OUT_DIR, "per-cluster")

CACHE_PATH = os.path.join(OUT_DIR, "network-background-pathway-sizes.json")
OUT_PATH = os.path.join(OUT_DIR, "network-background-reanalysis.csv")

WHOLE_PROTEOME_FDR_CUTOFF = 0.01  # a pathway is "originally significant" if its whole-proteome FDR is below this
FDR_CUTOFF = 0.01  # significance threshold applied to the recomputed (properly BH-corrected) network-background FDR

REACTOME_BASE = "https://reactome.org/AnalysisService"
# The Analysis Service returns HTTP 403 to clients that send no User-Agent.
HTTP_HEADERS = {"User-Agent": "MoDPA-network-background-ora/1.0", "Accept": "application/json"}
TOKEN_PARAMS = {"species": "Homo sapiens", "pageSize": "-1", "page": "-1",
                "resource": "TOTAL", "pValue": "1"}


def _query_via_reactome2py(background_proteins):
    """Preferred path: the reactome2py client, as used by the main notebook."""
    from reactome2py import analysis

    print(f"  using reactome2py (Reactome db version {analysis.db_version()})")
    result = None
    for attempt in range(5):
        try:
            time.sleep(3 * (attempt + 1))
            result = analysis.identifiers(ids=",".join(background_proteins))
            break
        except Exception as e:
            print(f"  attempt {attempt + 1}/5 failed: {e}")
    if result is None:
        raise RuntimeError("Could not reach Reactome Analysis Service")

    token = result["summary"]["token"]
    token_result = analysis.token(token, species="Homo sapiens", page_size="-1", page="-1",
                                  resource="TOTAL", p_value="1",
                                  min_entities=None, max_entities=None)
    return token_result


def _query_via_urllib(background_proteins):
    """Fallback path: the same two REST endpoints, called directly.

    Used when reactome2py is not installed. Issues the identical identifiers + token
    request pair and returns the same structure, so K values match the client path.
    """
    print("  reactome2py not available, calling the REST API directly")
    post_headers = dict(HTTP_HEADERS, **{"Content-Type": "text/plain"})
    url = f"{REACTOME_BASE}/identifiers/?interactors=false&pageSize=1&page=1"

    result = None
    for attempt in range(5):
        try:
            time.sleep(3 * (attempt + 1))
            req = urllib.request.Request(url, data=",".join(background_proteins).encode(),
                                         headers=post_headers, method="POST")
            result = json.load(urllib.request.urlopen(req, timeout=300))
            break
        except Exception as e:
            print(f"  attempt {attempt + 1}/5 failed: {e}")
    if result is None:
        raise RuntimeError("Could not reach Reactome Analysis Service")

    token = result["summary"]["token"]
    token_url = f"{REACTOME_BASE}/token/{token}?{urllib.parse.urlencode(TOKEN_PARAMS)}"
    req = urllib.request.Request(token_url, headers=HTTP_HEADERS)
    return json.load(urllib.request.urlopen(req, timeout=300))


def get_network_background_sizes(background_proteins, cache_path):
    """Return {Reactome stId: number of background_proteins annotated to it}.

    Cached to disk after the first (and only) Reactome API call this script makes.
    With the cache present the function is offline and deterministic. Delete the cache
    file to force a refresh, e.g. after a Reactome database update.
    """
    if os.path.isfile(cache_path):
        with open(cache_path) as fh:
            bg_found = json.load(fh)
        print(f"Using cached network background sizes: {len(bg_found)} pathways "
              f"({os.path.relpath(cache_path, BASE)}); no Reactome query needed")
        return bg_found

    print(f"Submitting {len(background_proteins)} network proteins to the Reactome "
          f"Analysis Service...")
    try:
        token_result = _query_via_reactome2py(background_proteins)
    except ImportError:
        token_result = _query_via_urllib(background_proteins)

    bg_found = {p["stId"]: p["entities"]["found"] for p in token_result["pathways"]}
    with open(cache_path, "w") as fh:
        json.dump(bg_found, fh)
    print(f"  {len(bg_found)} pathways annotated by >=1 network protein; cached to {cache_path}")
    return bg_found


def bh_adjust(pvals):
    """Benjamini-Hochberg adjustment. NaN entries (untestable pathways) pass through as NaN."""
    pvals = pvals.copy()
    mask = pvals.notna()
    if mask.sum() == 0:
        return pvals
    adj = np.full(len(pvals), np.nan)
    p = pvals[mask].values
    order = np.argsort(p)
    ranked = p[order]
    m = len(p)
    bh = ranked * m / np.arange(1, m + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0, 1)
    out = np.empty(m)
    out[order] = bh
    adj[np.where(mask)[0]] = out
    return pd.Series(adj, index=pvals.index)


def load_full_per_cluster_results(cluster_sizes):
    """Read every pathway Reactome tested (p_value="1") for each tested cluster."""
    tested_clusters = cluster_sizes.index[cluster_sizes.tested]
    frames = []
    for cluster in tested_clusters:
        cache = os.path.join(PER_CLUSTER_DIR, f"cluster{cluster}-reac.csv")
        if not os.path.isfile(cache):
            print(f"  WARNING: no cache for cluster {cluster} ({cache}), skipping")
            continue
        df = pd.read_csv(cache)
        df["cluster"] = cluster
        frames.append(df)
    full = pd.concat(frames, ignore_index=True)
    return full


def main():
    nodes = pd.read_csv(os.path.join(OUT_DIR, "nicely-formatted-nodes.csv"))
    cluster_sizes = pd.read_csv(os.path.join(OUT_DIR, "cluster-sizes.csv"), index_col=0)

    full = load_full_per_cluster_results(cluster_sizes)
    m_per_cluster = full.groupby("cluster").size()
    print(f"Loaded the full per-cluster candidate set: {len(full)} pathway-cluster pairs "
          f"across {full.cluster.nunique()} clusters (every pathway Reactome tested, not just "
          f"the ones significant under the whole-proteome background)")
    print(f"  per-cluster candidate-set size (m) used for BH correction: "
          f"{m_per_cluster.min()} to {m_per_cluster.max()}")

    background_proteins = sorted(set(nodes.UniAcc))
    N = len(background_proteins)
    print(f"Network background: {N} distinct proteins (all clusters, any size)")

    bg_found = get_network_background_sizes(background_proteins, CACHE_PATH)

    def recompute(row):
        K = bg_found.get(row.stId, 0)
        if K == 0:
            return np.nan  # not annotated to any network protein: untestable vs this background
        n = int(cluster_sizes.loc[row.cluster, "n_proteins"])
        k = min(int(row.entities_found), K)  # cluster proteins are a subset of network proteins
        return stats.hypergeom.sf(k - 1, N, K, n)

    full["p_network_raw"] = full.apply(recompute, axis=1)
    # BH correction here uses the TRUE per-cluster candidate set (every pathway Reactome
    # tested for that cluster) as the denominator, not just the originally-significant subset.
    full["FDR_network"] = full.groupby("cluster")["p_network_raw"].transform(bh_adjust)
    full["K_network_annotated"] = full.stId.map(bg_found).fillna(0).astype(int)
    full["N_network"] = N

    full.rename(columns={"FDR": "FDR_wholeproteome"}, inplace=True)
    cols = ["cluster", "stId", "name", "entities_found", "entities_tot",
            "FDR_wholeproteome", "K_network_annotated", "N_network",
            "p_network_raw", "FDR_network"]
    full[cols].to_csv(OUT_PATH, index=False)
    print(f"Wrote {OUT_PATH} ({len(full)} rows, the full per-cluster candidate set)")

    # Headline comparison: of the pathways that were significant under the whole-proteome
    # background, how many survive under the network background with a PROPERLY corrected FDR.
    orig_sig = full[full.FDR_wholeproteome < WHOLE_PROTEOME_FDR_CUTOFF].copy()
    n_pairs = len(orig_sig)
    n_testable = orig_sig.p_network_raw.notna().sum()
    n_still_sig = (orig_sig.FDR_network < FDR_CUTOFF).sum()
    print(f"\nOf the {n_pairs} pathway-cluster pairs significant at whole-proteome FDR<"
          f"{WHOLE_PROTEOME_FDR_CUTOFF} ({n_pairs - n_testable} untestable vs the network "
          f"background: pathway not annotated by any network protein):")
    print(f"  still significant at FDR<{FDR_CUTOFF} under the network background, PROPERLY "
          f"BH-corrected over the full per-cluster candidate set: "
          f"{n_still_sig} of {n_pairs} ({n_still_sig / n_pairs * 100:.1f}%)")
    mean_wp = (-np.log10(orig_sig.FDR_wholeproteome)).mean()
    mean_net = (-np.log10(orig_sig.FDR_network.clip(lower=1e-300))).mean()
    print(f"  mean -log10(adjusted p): {mean_wp:.2f} (whole-proteome) -> {mean_net:.2f} (network)")

    per_cluster = orig_sig.groupby("cluster").agg(
        n_orig=("stId", "size"),
        n_still_sig=("FDR_network", lambda x: (x < FDR_CUTOFF).sum()),
    )
    per_cluster["pct_surviving"] = (100 * per_cluster.n_still_sig / per_cluster.n_orig).round(1)
    per_cluster = per_cluster.join(m_per_cluster.rename("m_true_candidate_set"))
    lost_all = per_cluster.index[(per_cluster.n_orig > 0) & (per_cluster.n_still_sig == 0)].tolist()
    kept_any = per_cluster.index[per_cluster.n_still_sig > 0].tolist()
    print(f"  clusters that lose ALL originally significant pathways under the (properly "
          f"corrected) network background: {len(lost_all)} of {len(per_cluster)}: {lost_all}")
    print(f"  clusters retaining at least one: {kept_any}")
    print(f"\nPer-cluster summary:")
    print(per_cluster.to_string())


if __name__ == "__main__":
    main()
