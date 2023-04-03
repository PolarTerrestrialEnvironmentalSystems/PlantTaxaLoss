# PlantTaxaLoss

Code repository for project entitled 'Potential plant extinctions with the loss of the Pleistocene mammoth-steppe' 

## General content
The repository contains 8 files:
- 1_merge_and_arrange_all_cores.R
- 2_community_investigation.R
- 3_Resampling.R
- 4_Figures-preparation.R
- 5_Extirpation_and_extinction_curves.R
- 6_Correlation_explorations.R
- 7_Assessment_candidates_absent-modern-databases.R
- 8_Simulate_GBIF_coverage.R

## Details
To run the code you need to have the following data (future repository link)
Each R script should be run in the numerical order as they use outputs from previous scripts.

Please, make folders in your environment: "output", "plant_types", "Megafauna_data", "Temperature_estimates", "check_against_databases", "Simulate_sampling_GBIF".

For script: 1_merge_and_arrange_all_cores.R you need for each of the 8 investigated cores:
- core_identitylevel_90_merged_replicates.csv (wide dataframe with identified DNA sequenced and their reads count per sample in columns after OBITools3 as well as columns "best_identity", with the best identity score of the DNA sequence agaisnt the database; "scientific_name", "best_family", "best_genus", "best_species", holding information on the best taxonomic match against the database)
- core_agefile.xlsx (metadata file linking the sample names of the .csv file to metadata information of the samples e.g PCR_number, Extraction_number, depth, age, etc...)

For script: 2_community_investigation.R you need:
- only outputs from 1_merge_and_arrange_all_cores.R

For script: 3_Resampling.R you need:
- the outputs from previous scripts
- to have in "plant_types", a dataframe called plant_types_list.csv with 3 columns: "scientific_name", "type2" and "type3". type2 defines if the scientific name if a tree, schrub, herb, grass, moss or fern and type3 define if it is terrestrial or aquatic.

For script: 4_Figures-preparation.R you need:
- the outputs from previous scripts
- a type_family.csv dataframe with 2 columns: family (plant family) and a plant functional type associated ("shrub"; "shrub-tree"; "herb"; ...)

For script: 5_Extirpation_and_extinction_curves.R you need:
- the outputs from previous scripts
- megafauna.csv: a dataframe that list megafauna extinction with 8 columns: Common name; Scientific name; Area (Alaska/Yukon or Eurasia); Period of extinction; Approximate extinction time; timeslice_disappear (estimated timeslice of the approximate extinction time); Comment; link to source file. Here the list was created from Stuart et al. 2015.

For script: 6_Correlation_explorations.R you need:
- the outputs from previous scripts
- 2023_simulated_aggregated_timeslice.csv: a dataframe with 9 columns. Downloaded data from Dallmeyer et al., 2022. Here, for each timeslice we recovered from the  mean and median July temperature and precipitation anomalies between our studied time slices (from 23,000 to 1,000 cal. yrs BP). 
- Temperature_differences_Alaska_Siberia_pollen_based.csv: a dataframe with 9 columns. Downloaded from Herzschuh et al., 2022a, b. Here, we used pollen-based reconstructed temperature records from nine sites in the study area: Alut Lake; Elikchan 4 Lake; Hanging Lake; Joe Lake; Kaiyak Lake; Smorodinovoye Lake; Tukuto Lake; Zagoskin Lake; Lake Billyakh. 
- precipitation_differences_alaska_sibeia_pollen_based.csv: a dataframe with 9 columns.Downloaded from Herzschuh et al., 2022a, b. Here, we used pollen-based reconstructed temperature records from nine sites in the study area: Alut Lake; Elikchan 4 Lake; Hanging Lake; Joe Lake; Kaiyak Lake; Smorodinovoye Lake; Tukuto Lake; Zagoskin Lake; Lake Billyakh.
For the pollen-based reconstructed climate changes, both MAT and weighted average-partial least square regression.

For script: 7_Assessment_candidates_absent-modern-databases.R you need:
- the outputs from previous scripts
- 3 fasta files with trnlg/h databases built with OBITools 3: EMBL143 (Kanz et al., 2005); arctborbryo (Soininen et al., 2015; Sønstebø et al., 2010; Willerslev et al., 2014); PhyloNorway (Alsos et al., 2022). See details in Stoof-Leichsenring et al. in prep.

For script: 8_Simulate_GBIF_coverage.R you need:
- the outputs from previous scripts
- a dataframe GBIF_DB_coverage.csv available from Stoof-Leichsenring et al. in prep.

## References:
Alsos et al. 2022: "Postglacial species arrival and diversity build-up of northern ecosystems took millennia", https://doi.org/10.1126/sciadv.abo7434
Dallmeyer et al. 2022: "The deglacial forest conundrum", https://doi.org/10.1038/s41467-022-33646-6
Herzschuh et al. 2022a: "LegacyClimate 1.0: A dataset of pollen-based climate reconstructions from 2594 Northern Hemisphere sites covering the late Quaternary", https://doi.org/10.5194/essd-2022-38
Herzschuh et al. 2022b: "LegacyPollen 1.0: a taxonomically harmonized global late Quaternary pollen dataset of 2831 records with standardized chronologies", https://doi.org/10.5194/essd-14-3213-2022
Kanz et al. 2005: "The EMBL nucleotide sequence database", https://doi.org/10.1093/nar/gki098
Soininen et al. 2015: "Highly overlapping winter diet in two sympatric lemming species revealed by DNA metabarcoding", https://doi.org/10.1371/journal.pone.0115335 
Sønstebø et al. 2010: "Using next-generation sequencing for molecular reconstruction of past Arctic vegetation and climate", https://doi.org/10.1111/j.1755-0998.2010.02855.x
Stuart et al. 2015: "Late Quaternary megafaunal extinctions on the continents: a short review", https://doi.org/10.1002/gj.2633 
Willerslev et al. 2014: "Fifty thousand years of Arctic vegetation and megafaunal diet",  https://doi.org/10.1038/nature12921




  

