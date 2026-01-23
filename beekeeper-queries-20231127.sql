
select count(distinct(pxd_identifier)) from project;

--rawfile2pride
--\copy (SELECT project_id, file_name, number_msms, number_peptides, number_peptidoforms FROM assay JOIN run ON assay.id = run.assay_id) TO 'D:/20241107_rawfile2pride_v0113_metamORF.csv' DELIMITER ',' CSV HEADER;
--\copy (SELECT project_id, file_name, searchdatabase.name as search_db, number_msms, number_peptides, number_peptidoforms FROM run JOIN searchdatabase on searchdatabase.id=run.searchdatabase_id) TO 'D:/20241107_rawfile2pride_v0113_metamORF.csv' DELIMITER ',' CSV HEADER;
\copy (select project_id as PXD,file_name,searchdatabase.name as search_db,number_msms,number_peptides,number_peptidoforms from run join searchdatabase on searchdatabase.id = run.searchdatabase_id join assay on assay.id = run.assay_id) TO 'D:/20250217_rawfile2pride_v0113.csv' DELIMITER ',' CSV HEADER;

--IDs
\copy (SELECT peptidoform.id AS peptidoform_id, peptide_id, peptidoform.is_modified, peptide.sequence,modification.ionbot_name AS PTM_name, modificationsite.residue_location AS PTM_loc, modification.residue AS PTM_res, modification.classification FROM modificationsite join modification on modification.id = modificationsite.modification_id full outer join peptidoform on peptidoform.id = modificationsite.peptidoform_id join peptide on peptide.id = peptidoform.peptide_id) TO 'D:/20240120_Peptidoforms_IDs_v0110.csv' DELIMITER ',' CSV HEADER;

--COUNTS
--\copy (SELECT run.file_name,peptidoform.id as peptidoform_id,peptidoform.is_modified,peptidoform.has_unexpected_modification,COUNT(DISTINCT (psm.id)) as psm_counts,searchdatabase.name as searchdb FROM psm JOIN run ON run.id = psm.run_id FULL OUTER JOIN peptidoform ON peptidoform.id = psm.peptidoform_id join assay on assay.id = run.assay_id join peptide on peptide.id=peptidoform.peptide_id join searchdatabase on run.searchdatabase_id = searchdatabase.id GROUP BY peptidoform.id,run.file_name,peptidoform.is_modified,peptidoform.has_unexpected_modification,searchdatabase.name) TO 'D:/20240120_Peptidoforms_counts_v0110.csv' DELIMITER ',' CSV HEADER;
\copy (SELECT run.file_name, peptidoform.id as peptidoform_id, COUNT(DISTINCT(psm.id)) as psm_counts FROM psm JOIN run ON run.id = psm.run_id FULL OUTER JOIN peptidoform ON peptidoform.id = psm.peptidoform_id join assay on assay.id = run.assay_id join peptide on peptide.id=peptidoform.peptide_id join searchdatabase on run.searchdatabase_id = searchdatabase.id GROUP BY peptidoform.id,run.file_name,peptidoform.is_modified,peptidoform.has_unexpected_modification,searchdatabase.name) TO 'D:/20240120_Peptidoforms_counts_v0110.csv' DELIMITER ',' CSV HEADER;
\copy (SELECT run.file_name, peptidoform.id as peptidoform_id, COUNT(DISTINCT(psm.id)) as psm_counts FROM psm JOIN run ON run.id = psm.run_id FULL OUTER JOIN peptidoform ON peptidoform.id = psm.peptidoform_id join peptide on peptide.id=peptidoform.peptide_id join searchdatabase on run.searchdatabase_id = searchdatabase.id GROUP BY peptidoform.id,run.file_name,peptidoform.is_modified,peptidoform.has_unexpected_modification,searchdatabase.name) TO 'D:/20241122_Peptidoforms_counts_v0113_metamORF.csv' DELIMITER ',' CSV HEADER;

--MAPS (non mi servono piú)
--\copy (select peptide_id,peptide.sequence,is_unique,uniprot_accession,uniprot_name,searchdatabase.name as searchdb,peptide_start from peptide_protein join protein on protein.id = peptide_protein.protein_id join peptide on peptide.id=peptide_protein.peptide_id join protein_searchdatabase on peptide_protein.protein_id=protein_searchdatabase.protein_id join searchdatabase on protein_searchdatabase.searchdatabase_id = searchdatabase.id) TO 'D:/20240120_Peptides_mappings_v0110.csv' DELIMITER ',' CSV HEADER;
--\copy (select peptidoform_id,peptide_id,peptide.sequence,searchdatabase.name as searchdb from psm join peptidoform on peptidoform.id=psm.peptidoform_id join peptide on peptide.id=peptidoform.peptide_id join run on run.id=psm.run_id join searchdatabase on searchdatabase.id=run.searchdatabase_id where searchdatabase.name='human') TO 'D:/20240223_Peptides_mappings_v0110.csv' DELIMITER ',' CSV HEADER;

select
  psm.id as psm_id,
  peptidoform_id,
  precursor_mass,
  precursor_mass_error,
  precursor_mass - precursor_mass_error as corrected_mass
from
  psm
limit 1000

\copy (select psm.id as PSM_id,run_id,peptidoform_id,scan_number,usi,precursor_mass,precursor_charge,retention_time from psm where substr(usi, 8, 9) = 'PXD000498') TO 'D:/link-mgfquant-to-ionbot/PSMs_PXD000498_v0110.csv' DELIMITER ',' CSV HEADER;
\copy (select psm.id as psm_id,run_id,peptidoform_id,scan_number,usi,precursor_mass,precursor_charge,retention_time from psm order by psm.id) TO 'D:/20241218_PSMs_v0113_metamORF.csv' DELIMITER ',' CSV HEADER;


-- partial
\copy (SELECT run.file_name, peptidoform.id as peptidoform_id, COUNT(DISTINCT(psm.id)) as psm_counts FROM psm JOIN run ON run.id = psm.run_id FULL OUTER JOIN peptidoform ON peptidoform.id = psm.peptidoform_id join assay on assay.id = run.assay_id join peptide on peptide.id=peptidoform.peptide_id join searchdatabase on run.searchdatabase_id = searchdatabase.id WHERE assay.project_id='PXD005833' GROUP BY peptidoform.id,run.file_name,peptidoform.is_modified,peptidoform.has_unexpected_modification,searchdatabase.name) TO 'D:/20251003_Peptidoforms_counts_v0113.csv' DELIMITER ',' CSV HEADER;

