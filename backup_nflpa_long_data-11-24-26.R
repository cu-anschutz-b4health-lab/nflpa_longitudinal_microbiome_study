#File to create RData for publication 


#Loads in meta data for all timepoints data
meta_mindsource = read.csv('/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/Mindsource/Data/Mindsource_metadata.csv')

#Adding LIBRA scores to meta data
#LIBRA Data
libra_meta = read.csv('/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/Mindsource/Data/LIBRA_MINDSOURCE.csv')
libra_meta_bl = read.csv('/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/Mindsource/Data/LIBRA_MINDSOURCE_bl48hrs.csv')

participant_list_mindsource = read.csv('/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/Mindsource/Data/Participant_list.csv')

taxa_table_1 = read.table('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0001_DeepSeq_Jul_2023_Diversigen/filtered taxonomy tables/taxatables-aggregated-absolute/taxatable-species-absolute.tsv', header =T)
taxa_table_2 = read.table('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0003_DeepSeq_May_2024_Diversigen/filtered taxonomy tables/taxatables-by-level/taxatable-species-absolute.tsv', header=T)


path_table_1 = read.csv('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0001_DeepSeq_Jul_2023_Diversigen/filtered functional tables/uniref-pathway-filtered-absolute.csv', header =T, row.names=1)
path_table_2 = read.csv('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0003_DeepSeq_May_2024_Diversigen/filtered functional tables/uniref-pathway-filtered-absolute.csv', header =T,row.names = 1)

ec_table_1 = read.table('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0001_DeepSeq_Jul_2023_Diversigen/filtered functional tables/sample_ec4.tsv', header =T,sep="\t", check.names = F, row.names = 1)

ec_table_2 = read.table('/Volumes/B4_Backup/Metagemomic-sequences/MQ3201_0003_DeepSeq_May_2024_Diversigen/filtered functional tables/sample_ec4.tsv', header =T,sep="\t", check.names = F, row.names = 1)

filter_ec_table = read.csv('/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/Mindsource/Data/ec_number_list_final-09-29-26.csv')

ec_csv = read.table('/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/Mindsource/Analysis/MaAsLin3/MaAsLin3_EC_SCFA_2026-09-29/significant_results.tsv', header=T)

d <- as.matrix(read.table("rpca_dist/distance-matrix.tsv", header = TRUE,
                          row.names = 1, sep = "\t", check.names = FALSE))

txt  <- readLines(bzfile("metaphlan2_tree.nwk.bz2"), warn = FALSE)

save.image("/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/Mindsource/Data/NFLPA_Longitudinal_Microbiome-10-26-26.RData")
