###LIBRA Mndsource Microbiome 
library(vegan)
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ANCOMBC)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(vegan)
library(knitr)
library(forcats)
library(tidyverse)
library(DT)
library(outliers)
library(beeswarm)
library(WGCNA)
library(flextable)
library(magrittr)
library(dplyr)
library(tidyheatmaps)
library(devtools)
library(ggvegan)
library(janitor)
library(stringr)
library(gtsummary)
#microbiomestat object
library(microbiome)
library(MicrobiomeStat)
library(GUniFrac)
BiocManager::install("biomformat")

library(biomformat)
#library(devtools)
#devtools::install_github(
#  "biobakery/perma-brioche",
#  build_vignettes = TRUE
#)
library(permabrioche)

#save and load data
load("/Path/To/NFLPA_Longidutinal_Microbiome.RData")

my_colors = c("#A77903", 
              "#1B3558",
              "#3DC8E0",
              "#5A595A")
options(
  ggplot2.discrete.fill = my_colors,
  ggplot2.discrete.colour = my_colors
)
theme_set(theme_bw())

#ste up path for later use

#Loads in meta data for all timepoints data
meta_mindsource

#Adding LIBRA scores to meta data
#LIBRA Data
libra_meta 
libra_meta_bl
#merge Libra information together
libra_meta_bl2 = libra_meta_bl %>%
  mutate(Ever_Homeless = case_when(homelessn >= 1 ~ "Yes",
                                   homelessn == 0 ~ "No",
                                   T~"Unknown") ) %>%
  mutate(LIBRA_Half = case_when(Total_LIBRA_Score >= median(Total_LIBRA_Score) ~"Upper",
                                Total_LIBRA_Score < median(Total_LIBRA_Score) ~"Lower")) %>%
  dplyr::rename(Participant = participant) 
libra_meta_bl2$Education[is.na(libra_meta_bl2$Education)] <- "Unknown"

libra_meta_bl2$Education = factor(libra_meta_bl2$Education, levels=c("High school diploma or equivalent", "Some college, no degree", "Associate's degree", 
                                                                     "Bachelor's degree", "Master's degree", "Unknown"))

libra_meta_bl2$Ever_Homeless = factor(libra_meta_bl2$Ever_Homeless, levels=c("No", "Yes", "Unknown"))
libra_meta_bl2$sex = factor(libra_meta_bl2$sex)

#Filters mindsource meta data to only those in the study
participant_list_mindsource

meta_mindsource_filter = meta_mindsource %>% filter(Participant %in% participant_list_mindsource$Participants)
meta_mindsource_filter_update_ID = meta_mindsource_filter %>% 
  mutate(Sample_ID = str_replace_all(SampleID, fixed("."), "_")) %>%
  mutate(Sample_ID = str_replace_all(Sample_ID, fixed("-"), "_"))


meta_and_libra = meta_mindsource_filter_update_ID %>% mutate(Time_num = case_when(Time_Point == "BL" ~ 0,
                                                                                  Time_Point == "48h" ~ 0.07,
                                                                                  Time_Point == "2m" ~ 2,
                                                                                  Time_Point == "3m" ~ 3,
                                                                                  Time_Point == "4m" ~ 4,
                                                                                  Time_Point == "5m" ~ 5,
                                                                                  Time_Point == "6m" ~ 6,
                                                                                  Time_Point == "7m" ~ 7,
                                                                                  Time_Point == "8m" ~ 8,
                                                                                  Time_Point == "9m" ~ 9,
                                                                                  Time_Point == "10m" ~ 10,
                                                                                  Time_Point == "11m" ~ 11,
                                                                                  Time_Point == "12m" ~ 12))


meta_and_libra$Time_Point = factor(meta_and_libra$Time_Point, levels = c("BL", "48h","2m", "3m", "4m", "5m", "6m", "7m","8m","9m","10m","11m","12m" ))
meta_and_libra_final = full_join(meta_and_libra, libra_meta_bl2, by="Participant") %>% filter(!Sample_ID %in% c('19_1307_P24_8m_g', '19_1307_P34_2m_g'))

meta_and_libra_final %>% group_by(Participant) %>% summarize(count=n())
######Table 1 Demographics#####
#Creates table1 and saves to the working directory
table1 <-
  tbl_summary(
    libra_meta_bl2,
    include = c(Total_LIBRA_Score, age, sex, race, Ever_Homeless, BMI_Current, Education),
    by = LIBRA_Half, # split table by group
    missing = "ifany", # don't list missing data separately
    statistic = all_continuous() ~ "{mean} (± {sd})" # Custom format for mean and SD
  ) |> 
  add_n() |> # add column with total number of non-missing observations
  add_p() |> # test for a difference between groups
  add_overall() %>%
  modify_header(label = "**Variable**", stat_1 = "**Lower Brain Health Risk Group**  \nN = {n}", stat_2 = "**Higher Brain Health Risk Group**  \nN = {n}") |> # update the column header
  bold_labels() |>
  modify_spanning_header( ~ "**Table 1 Summary of Demographics for Mindsource Microbiome Study**") |>
  as_flex_table() |> 
  flextable::save_as_docx(path = paste0(path, "Table_1_Summary_Mindsource.docx"))

summary(lm(Total_LIBRA_Score ~ age, data = libra_meta_bl2))
summary(lm(Total_LIBRA_Score ~ sex, data = libra_meta_bl2))
libra_meta_bl2$race[is.na(libra_meta_bl2$race)] <- "Unknown"
libra_meta_bl2$race = factor(libra_meta_bl2$race, levels=c("White", "Black", "Other", "Unknown"))
summary(lm(Total_LIBRA_Score ~ race, data = libra_meta_bl2))
summary(lm(Total_LIBRA_Score ~ Ever_Homeless, data = libra_meta_bl2))
summary(lm(Total_LIBRA_Score ~ BMI_Current, data = libra_meta_bl2))
summary(lm(Total_LIBRA_Score ~ Education, data = libra_meta_bl2))






libra_meta_bl2$Healthy_Diet =factor(libra_meta_bl2$Healthy_Diet)
libra_meta_bl2$Low_Alcohol =factor(libra_meta_bl2$Low_Alcohol)
libra_meta_bl2$Physical_Inactivity =factor(libra_meta_bl2$Physical_Inactivity)
libra_meta_bl2$Obesity =factor(libra_meta_bl2$Obesity)
libra_meta_bl2$Depression =factor(libra_meta_bl2$Depression)
libra_meta_bl2$Heart_Disease =factor(libra_meta_bl2$Heart_Disease)
libra_meta_bl2$CKD =factor(libra_meta_bl2$CKD)
libra_meta_bl2$Diabetes =factor(libra_meta_bl2$Diabetes)
libra_meta_bl2$Hypercholesterolemia =factor(libra_meta_bl2$Hypercholesterolemia)
libra_meta_bl2$Smoking =factor(libra_meta_bl2$Smoking)
libra_meta_bl2$Hypertension =factor(libra_meta_bl2$Hypertension)


#Creates table 2 for LIBRA Score characteristics
table2 <-
  tbl_summary(
    libra_meta_bl2,
    include = c(Total_LIBRA_Score, Healthy_Diet, Low_Alcohol, Physical_Inactivity, Obesity, Depression,Heart_Disease,CKD,Diabetes,Hypercholesterolemia,Smoking,Hypertension),
    by = LIBRA_Half, # split table by group
    missing = "ifany", # don't list missing data separately
    statistic = all_continuous() ~ "{mean} (± {sd})" # Custom format for mean and SD
  ) |> 
  add_n() |> # add column with total number of non-missing observations
  add_p() |> # test for a difference between groups
  add_overall() %>%
  modify_header(label = "**Variable**", stat_1 = "**Lower Brain Health Risk Group**  \nN = {n}", stat_2 = "**Higher Brain Health Risk Group**  \nN = {n}") |> # update the column header
  bold_labels() |>
  modify_spanning_header( ~ "**Table 2 Summary of Characteristics for Mindsource Microbiome Study**") |>
  as_flex_table() |> 
  flextable::save_as_docx(path = paste0(path, "Table_2_Summary_Mindsource.docx"))
######


####build out the phyloseq object
taxa_table_1 
trans_taxa_table_1 = t(data.frame(taxa_table_1))
trans_taxa_table_1 = trans_taxa_table_1 %>% row_to_names(row_number = 1)
trans_taxa_table_1 = data.frame(trans_taxa_table_1)
trans_taxa_table_1$Sample_ID = rownames(trans_taxa_table_1)
rownames(trans_taxa_table_1) = NULL
trans_taxa_table_1 = trans_taxa_table_1 %>% dplyr::select(Sample_ID, everything())


taxa_table_2
trans_taxa_table_2 = t(data.frame(taxa_table_2))
trans_taxa_table_2 = trans_taxa_table_2 %>% row_to_names(row_number = 1)
trans_taxa_table_2 = data.frame(trans_taxa_table_2)
trans_taxa_table_2$Sample_ID = rownames(trans_taxa_table_2)
rownames(trans_taxa_table_2) = NULL
trans_taxa_table_2 = trans_taxa_table_2 %>% dplyr::select(Sample_ID, everything())


taxa_table_all = bind_rows(trans_taxa_table_1,trans_taxa_table_2)
taxa_table_all[is.na(taxa_table_all)] <- 0
taxa_table_all_ID_Update = taxa_table_all %>% mutate(Sample_ID = str_replace(Sample_ID, "X", ""))
taxa_table_all_ID_Update$Sample_ID
taxa_table_all_filtered = taxa_table_all_ID_Update %>% dplyr::filter(Sample_ID %in% meta_mindsource_filter_update_ID$Sample_ID & (!Sample_ID %in% c('19_1307_P24_8m_g', '19_1307_P34_2m_g')))
rownames(taxa_table_all_filtered) = taxa_table_all_filtered$Sample_ID
taxa_table_all_filtered_mat = taxa_table_all_filtered %>% select(-Sample_ID) %>% select(where(~ any(. != 0)))
taxa_table_all_filtered_mat[] = lapply(taxa_table_all_filtered_mat, as.numeric)
taxa_table_all_filtered_mat = as.matrix(taxa_table_all_filtered_mat)

#otu table
otu_table = t(taxa_table_all_filtered_mat)
#taxa names
tax_names = rownames(otu_table)
#Update
tax_names_2 = data.frame(tax_names) %>%
  separate(tax_names, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\.")

rownames(tax_names_2) = tax_names

tax_names_2_mat = as.matrix(tax_names_2)

OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(tax_names_2_mat)
sample = data.frame(meta_and_libra_final)
rownames(sample) = sample$Sample_ID
sampledata = sample_data(sample)

physeq = phyloseq(OTU, TAX, sampledata)
library("ape")
physeq_object_mindsource_bac <- subset_taxa(physeq, Kingdom == "k__Bacteria")

#random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
#physeq_object_mindsource = phyloseq(OTU, TAX, sampledata, random_tree)
#physeq_object_mindsource_bac <- subset_taxa(physeq_object_mindsource, Kingdom == "k__Bacteria")
txt  <- readLines(bzfile("metaphlan2_tree.nwk.bz2"), warn = FALSE)
tree <- read.tree(text = txt)
tree$tip.label <- sub("\\|GCF_[0-9]+$", "", tree$tip.label)
sp_of <- function(x) sub("^.*[.|](s__.*)$", "\\1", x)        # "...s__Species" -> "s__Species"

tip_sp <- sp_of(tree$tip.label)
tree   <- drop.tip(tree, which(duplicated(tip_sp)))          # one tip per species

# species in your data -> your taxa name(s) (a species can have >1 row in your table)
lookup <- split(taxa_names(physeq_object_mindsource_bac), sp_of(taxa_names(physeq_object_mindsource_bac)))

# keep only tips whose species is in your data, then rename them
tree   <- keep.tip(tree, which(tip_sp %in% names(lookup)))
tip_sp <- sp_of(tree$tip.label)
tree$tip.label <- vapply(lookup[tip_sp], `[`, character(1), 1)

# Read the Newick tree
# Keep only taxa that are in the phyloseq object
tree_pruned <- drop.tip(
  tree,
  setdiff(tree$tip.label, taxa_names(physeq_object_mindsource_bac))
)

# Add the tree to the phyloseq object
ps_tree <- merge_phyloseq(physeq_object_mindsource_bac, phy_tree(tree_pruned))

# Check that the tree was added
ps
#random_tree2 = rtree(ntaxa(physeq_object_mindsource_bac), rooted=TRUE, tip.label=taxa_names(physeq_object_mindsource_bac))
physeq_object_mindsource_final = phyloseq(OTU, TAX, sampledata)

####alpha diversity
?estimate_richness
alpha_meta = estimate_richness(physeq_object_mindsource_final,measures=c("Shannon","Observed", "Chao1", "InvSimpson")) %>%
  rownames_to_column(var = "Sample_ID") 
alpha_meta2 = alpha_meta %>% mutate(Sample_ID = str_replace(Sample_ID, "X", ""))
alpha_meta2$Sample_ID

alpha_meta_final = left_join(alpha_meta2, as.data.frame(sample_data(physeq_object_mindsource_final)), by = "Sample_ID") 

library(DescTools)
mean(table(alpha_meta_final$Participant))
Mode(table(alpha_meta_final$Participant))
sd(table(alpha_meta_final$Participant))

#colvec <- c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A")

#Spaghetti plots
shannon_line = ggplot(alpha_meta_final, aes(x=Time_num, y=Shannon,group=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = alpha_meta_final)+
  labs(tag = "A", x="Time (Months)", y="Shannon Index", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #scale_y_continuous(limits = c(0,max(alpha_meta_final$Observed)), breaks =c(0,1, 2, 3, max(alpha_meta_final$Observed)))+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  ) 


shannon_linetest = ggplot(alpha_meta_final, aes(x=Time_num, y=Shannon,group=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = alpha_meta_final)+
  labs(tag = "A", x="Time (Months)", y="Shannon Index", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 14, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  facet_wrap(~Participant)+
  #scale_y_continuous(limits = c(0,max(alpha_meta_final$Observed)), breaks =c(0,1, 2, 3, max(alpha_meta_final$Observed)))+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  ) 


alpha_meta_final %>% 
  filter((Time_Point == "BL" | Time_Point == "48H")) %>% # & LIBRA_Half == "Lower") %>% 
  summarise(value = mean(Shannon))
alpha_meta_final %>% 
  filter((Time_Point == "BL" | Time_Point == "48H")) %>% # & LIBRA_Half == "Lower") %>% 
  summarise(value = sd(Shannon))

t.test(Shannon ~ LIBRA_Half, data = alpha_meta_final %>%  filter((Time_Point == "BL" | Time_Point == "48H") ))

t.test(Shannon ~ LIBRA_Half, data = alpha_meta_final %>%  filter((Time_num == 12 |Time_num == 11)  ))


alpha_meta_final %>% 
  filter((Time_Point == "BL" | Time_Point == "48H")  & LIBRA_Half == "Upper") %>% 
  summarise(value = mean(Shannon))
alpha_meta_final %>% 
  filter((Time_Point == "BL" | Time_Point == "48H")  & LIBRA_Half == "Upper") %>% 
  summarise(value = sd(Shannon))


alpha_meta_final %>% 
  filter((Time_num == 12 |Time_num == 11) )%>%# & LIBRA_Half == "Lower") %>% 
  summarise(value = mean(Shannon))
alpha_meta_final %>% 
  filter((Time_num == 12 |Time_num == 11) )%>%# & LIBRA_Half == "Lower") %>% 
  summarise(value = sd(Shannon))


alpha_meta_final %>% 
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Upper") %>% 
  summarise(value = mean(Shannon))
alpha_meta_final %>% 
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Upper") %>% 
  summarise(value = sd(Shannon))


observed_line = ggplot(alpha_meta_final, aes(x=Time_num, y=Observed,group=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = alpha_meta_final)+
  labs(tag = "A", x="Time (Months)", y="Observed OTUs", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 12, face="bold"),
        plot.tag = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #scale_y_continuous(limits = c(0,max(alpha_meta_final$Observed)), breaks =c(0, 50, 100, max(alpha_meta_final$Observed)))+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Better Score", "Worse Score")
  )


spaghetti_plot = ggarrange(observed_line, shannon_line)

ggsave(paste0(path, "alpha_spaghetti_plot_", Sys.Date(), ".jpeg"), spaghetti_plot, height = 8, width = 15, units="in",device="jpeg")


###LIBRA Score
#meta_and_alpha$Time_bin = factor(meta_and_alpha$Time_bin, levels = c("Base", "2-4m", "5-8m","9-12m"))
shannon_libra_bin = ggplot(alpha_meta_final, aes(x=reorder(Participant, Total_LIBRA_Score), y=Shannon, color=LIBRA_Half))+
  geom_boxplot()+
  geom_point()+
  labs(tag = "B", x="Participants \n(Ordered by LIBRA Score)", y="Shannon Index", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_x_reverse()+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )


observed_libra_bin = ggplot(alpha_meta_final, aes(x=reorder(Participant, Total_LIBRA_Score), y=Observed, color=LIBRA_Half))+
  geom_boxplot()+
  geom_point()+
  labs(tag = "B", x="Participants \n(Ordered by LIBRA Score)", y="Observed OTUs", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 12, face="bold"),
        plot.tag = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  #scale_x_reverse()+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Better Score", "Worse Score")
  )

alpha_libra_bin_plot = ggarrange(shannon_libra_bin, observed_libra_bin)

ggsave(paste0(path, "alpha_libra_bin_plot_", Sys.Date(), ".jpeg"), alpha_libra_bin_plot, height = 6, width = 10, units="in",device="jpeg")

shannon_merged_spaghetti_box = ggarrange(shannon_line, shannon_libra_bin)
otu_merged_spaghetti_box = ggarrange(observed_line, observed_libra_bin)

ggsave(paste0(path, "shannon_merged_spaghetti_box_", Sys.Date(), ".jpeg"), shannon_merged_spaghetti_box, height = 6, width = 14, units="in",device="jpeg")
ggsave(paste0(path, "otu_merged_spaghetti_box_", Sys.Date(), ".jpeg"), otu_merged_spaghetti_box, height = 8, width = 14, units="in",device="jpeg")


#alpha lmer
library(broom.mixed)
library(modelsummary)
library(lme4)
library(lmerTest)

alpha_meta_final$Ever_Homeless = factor(alpha_meta_final$Ever_Homeless, levels=c("No", "Unknown", "Yes"))


model_score_shan <- lmer(Shannon ~ Time_num*Total_LIBRA_Score + age+ as.factor(sex) + BMI_Current +as.factor(Ever_Homeless)+ (1 | Participant), data = alpha_meta_final)
#model_half_shan <- lmer(Shannon ~ Time_num*LIBRA_Half + age+sex + Ever_Homeless+ (1 | Participant), data = alpha_meta_final)
summary(model_score_shan)
anova(model_score_shan)
#summary(model_half_shan)
modelsummary(list(model_score_shan, model_half_shan), output = paste0(path, "shannon_model_summary.docx"))

sup_table_shan_model = flextable(tidy(model_score_shan))
sup_table_shan_model = set_caption(sup_table_shan_model, caption = "Supplementary Table x. Mixed Effects Model Shannon Index") %>% autofit()
save_as_docx(sup_table_shan_model, path = paste0(path, "Supplementary_Table_shan_mod._", Sys.Date(), ".docx"))


library(gamm4)

model_score_shan_gamm <- gamm4(
  Shannon ~ t2(Time_num, Total_LIBRA_Score) + age + as.factor(sex) + 
    BMI_Current + as.factor(Ever_Homeless),
  random = ~ (1 | Participant),
  data = alpha_meta_final
)

summary(model_score_shan_gamm$gam)   # smooth-term significance, deviance explained, etc.
summary(model_score_shan_gamm$mer)   # the underlying lmer fit, if you want the random-effect variance
plot(model_score_shan_gamm$gam, pages = 1)  # visualize the fitted smooth surface


library(lme4)
library(lmerTest)
model_score_ob <- lmer(Observed ~ Time_num*Total_LIBRA_Score + age+sex +BMI_Current+ Ever_Homeless+ (1 | Participant), data = alpha_meta_final)
model_half_ob <- lmer(Observed ~ Time_num*LIBRA_Half + age+sex + Ever_Homeless+ (1 | Participant), data = alpha_meta_final)
summary(model_score_ob)
summary(model_half_ob)
modelsummary(list(model_score_shan, model_half_shan), output = paste0(path, "obs_model_summary.docx"))

sup_table_obs_model = flextable(tidy(model_score_ob))
sup_table_obs_model = set_caption(sup_table_obs_model, caption = "Supplementary Table x. Mixed Effects Model Observed OTUs") %>% autofit()
save_as_docx(sup_table_obs_model, path = paste0(path, "Supplementary_Table_obs_mod._", Sys.Date(), ".docx"))




model_score_obs_gamm <- gamm4(
  Observed ~ t2(Time_num, Total_LIBRA_Score) + age + as.factor(sex) + 
    BMI_Current + as.factor(Ever_Homeless),
  random = ~ (1 | Participant),
  data = alpha_meta_final
)

summary(model_score_obs_gamm$gam)   # smooth-term significance, deviance explained, etc.
summary(model_score_obs_gamm$mer)   # the underlying lmer fit, if you want the random-effect variance
plot(model_score_obs_gamm$gam, pages = 1)  # visualize the fitted smooth surface


##Coefficient of variation
calculate_cv <- function(x, na.rm = TRUE) {
  (sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)) * 100
}


cv_by_libra_time <- alpha_meta_final %>%
  group_by(Time_num) %>%
  summarize(
    Mean_Value = mean(Shannon, na.rm = TRUE),
    SD_Value = sd(Shannon, na.rm = TRUE),
    CV_Percentage = calculate_cv(Shannon)
  )
cv_by_libra_time   %>% summarize(mean(CV_Percentage))

cv_by_libra_half <- alpha_meta_final %>%
  group_by(LIBRA_Half) %>%
  summarize(
    Mean_Value = mean(Shannon, na.rm = TRUE),
    SD_Value = sd(Shannon, na.rm = TRUE),
    CV_Percentage = calculate_cv(Shannon)
  )

# View the results
print(cv_by_group)


#####Form microbiomestat object for volatility tests
#removed from final analysis
data.obj <- mStat_convert_phyloseq_to_data_obj(physeq_object_mindsource_final)

alpha_volatility_test_results <- generate_alpha_volatility_test_long(
  data.obj = data.obj,
  alpha.obj = NULL,
  alpha.name = c("shannon","observed_species"),
  time.var = "Time_num",
  subject.var = "Participant",
  group.var = "LIBRA_Half",
  adj.vars ="Total_LIBRA_Score"
)


alpha_shannon_vol = flextable(alpha_volatility_test_results$shannon)
save_as_docx(alpha_shannon_vol, path = paste0(path,"Supplementary_Table_alpha_shannon_vol_results.docx"))
alpha_obs_vol = flextable(alpha_volatility_test_results$observed_species)
save_as_docx(alpha_obs_vol, path = paste0(path,"Supplementary_Table_alpha_obs_vol_results.docx"))


####Add additional tests for additive models

#beta
#####Beta Diversity from Phyloseq




meta_and_alpha_libra_beta_long = meta_and_libra_final %>% 
  group_by(Participant) %>%
  mutate(Time_num = case_when(Time_Point == "BL" ~ 0,
                              Time_Point == "48h" & any(Time_Point =="BL") ~ 0.07,
                              Time_Point == "48h" ~0,
                              Time_Point == "2m" ~ 2,
                              Time_Point == "3m" ~ 3,
                              Time_Point == "4m" ~ 4,
                              Time_Point == "5m" ~ 5,
                              Time_Point == "6m" ~ 6,
                              Time_Point == "7m" ~ 7,
                              Time_Point == "8m" ~ 8,
                              Time_Point == "9m" ~ 9,
                              Time_Point == "10m" ~ 10,
                              Time_Point == "11m" ~ 11,
                              Time_Point == "12m" ~ 12))


print(meta_and_alpha_libra_beta_long %>% select(Participant, Time_num, LIBRA_Half), n=Inf)

sample = data.frame(meta_and_alpha_libra_beta_long)
rownames(sample) = sample$Sample_ID
sampledata = sample_data(sample)

physeq3 = phyloseq(OTU, TAX, sampledata, tree_pruned)


d 
# Convert matrix to dist object
rpca_dist_obj <- as.dist(d)

# PCoA
pcoa <- ape::pcoa(rpca_dist_obj)
pcoa_df <- as.data.frame(pcoa$vectors[, 1:2])

pcoa_df$Sample_ID <- rownames(pcoa_df)
colnames(pcoa_df)[1:2] <- c("PCoA1", "PCoA2")
pc1 <- round(pcoa$values$Relative_eig[1] * 100, 2)
pc2 <- round(pcoa$values$Relative_eig[2] * 100, 2)
pcoa_df <- pcoa_df %>%
  left_join(alpha_meta_final, by = "Sample_ID")

phylo_rpca = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = LIBRA_Half)) +
  geom_point(size = 3) +
  labs(
    x = paste0("PCoA1 (", pc1, "%)"),
    y = paste0("PCoA2 (", pc2, "%)"),
    color = "LIBRA_Half"
  ) +
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw() +
  theme(text = element_text(size = 16, face="bold"))+
  labs(color="Brain Health Group \n(LIBRA Score)") +
  #theme(plot.tag.position = c(.885, .98), plot.tag = element_text(size=7.25))+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )




phylo_rpca_part = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Participant)) +
  geom_point(size = 3) +
  labs(
    x = paste0("PCoA1 (", pc1, "%)"),
    y = paste0("PCoA2 (", pc2, "%)"),
    color = "Participant"
  ) +
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw() +
  theme(text = element_text(size = 16, face="bold"))



ord <- read_qza("biplot.qza")$data$Vectors   # SampleID, PC1, PC2, ...
plt <- left_join(ord, meta, by = c("SampleID" = "sample")) |> arrange(subject, days)

ggplot(plt, aes(PC1, PC2, color = group, group = subject)) +
  geom_path(alpha = .4, arrow = arrow(length = unit(.1, "cm"))) +
  geom_point(aes(size = days_since_base == 0)) +
  scale_size_manual(values = c(1.5, 3), guide = "none")



data.obj <- mStat_convert_phyloseq_to_data_obj(physeq3)
otu_mat <- data.obj$feature.tab          # taxa in rows, samples in columns
aitch <- vegdist(t(otu_mat), method = "robust.aitchison")   # pseudocount for zeros

# or: vegdist(t(otu_mat), method = "robust.aitchison")   # zero-robust, no pseudocount
dist.obj = mStat_calculate_beta_diversity(data.obj,
                                          dist.name = c("BC", "Jaccard", "UniFrac", "GUniFrac", "WUniFrac", "JS"))
dist.obj$Aitchison <- as.matrix(aitch)
dist.obj$phylo_rpca <- as.matrix(d)


ps <- prune_samples(sample_sums(ps) > 0, ps)

# samples x taxa count matrix
otu <- (as(otu_table(physeq3), "matrix"))
#if (taxa_are_rows(ps)) otu <- t(otu)
write.csv(otu, "counts.csv")

write.table(
  otu,
  "counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = T
)
# metadata: needs a subject ID and a NUMERIC time (e.g. days)
meta <- data.frame(sample_data(physeq3))
write.csv(meta, "meta.csv")

# taxonomy and tree (for interpretation later)
write.csv(as.data.frame(as(tax_table(physeq3), "matrix")), "taxa.csv")
write.tree(phy_tree(physeq3), "tree.nwk")


feature_ids <- rownames(otu)

# Convert the "." separators to ";" separators
taxonomy <- gsub("\\.", ";", feature_ids)

# Create QIIME 2 taxonomy table
taxonomy_table <- data.frame(
  `Feature ID` = feature_ids,
  Taxon = taxonomy,
  check.names = FALSE
)

# Write QIIME 2-compatible TSV
write.table(
  taxonomy_table,
  "taxonomy.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

bad_id <- "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Dermatophilaceae.g__Dermatophilaceae_unclassified.s__Other"

bad_id %in% taxonomy_table$`Feature ID`
bad_id %in% tree_pruned$tip.label



otu_all <- as(otu_table(physeq_object_mindsource_final), "matrix")
if (taxa_are_rows(physeq_object_mindsource_final)) otu_all <- t(otu_all)
robust.aitchison.dist <- vegdist(otu_all, method = "robust.aitchison")
#changes data to abundance counts and moves 0 to a low value
phylo_mindsource_beta_filter_tranform = transform_sample_counts(ps_tree, function(x) x/sum(x))

ps_transform = transform_sample_counts(ps_tree, function(x) 1E6 * x/sum(x))
ps.meta_mindsource = meta(ps_tree)


# ---------------------------------------------------------------------------
# 0. Inputs: swap in whichever distance / metadata you're using
# ---------------------------------------------------------------------------
dist_obj2 <- robust.aitchison.dist        # any `dist` object (UniFrac, Aitchison, ...)
meta     <- ps.meta_mindsource           # sample metadata

candidates <- c("Time_num", "Total_LIBRA_Score", "BMI_Current",
                "age", "sex", "Ever_Homeless", "Participant")


meta$sex           <- as.factor(meta$sex)
meta$Ever_Homeless <- as.factor(meta$Ever_Homeless)

# Line up metadata rows with the distance matrix, then keep complete cases only
# (adonis2 / stepwise comparisons need the same samples in every model)
meta <- meta[labels(dist_obj), , drop = FALSE]
keep <- complete.cases(meta[, c(candidates, "Participant")])
meta_cc <- meta[keep, , drop = FALSE]
dist_cc <- as.dist(as.matrix(dist_obj)[keep, keep])


norm_term <- function(x) {
  vapply(strsplit(x, ":"), function(p) paste(sort(p), collapse = ":"), character(1))
}


# ---------------------------------------------------------------------------
# Backward elimination: all covariates as main effects, plus ONLY the
# interaction(s) you name in `interactions`.
# ---------------------------------------------------------------------------
backward_permanova2 <- function(D, data, candidates,
                                interactions = "Time_num:Total_LIBRA_Score",
                                alpha = 0.1, permutations = 999,
                                strata = NULL, seed = 64, max_steps = 50) {
  current <- c(candidates, interactions)
  history <- data.frame()
  
  fit_margin <- function(terms) {
    set.seed(seed)
    adonis2(reformulate(terms, response = "D"), data = data,
            permutations = permutations, by = "margin", strata = strata)
  }
  
  for (i in seq_len(max_steps)) {
    a <- fit_margin(current)
    tab <- data.frame(term = rownames(a), R2 = a$R2, F = a$F, p = a$`Pr(>F)`,
                      stringsAsFactors = FALSE)
    tab <- tab[!tab$term %in% c("Residual", "Total"), ]
    if (nrow(tab) == 0) break
    
    tab$p_rank <- ifelse(is.na(tab$p), Inf, tab$p)
    worst <- tab[order(-tab$p_rank, tab$R2), ][1, ]
    
    if (!is.na(worst$p) && worst$p < alpha) break     # all remaining terms significant
    
    # Match the reported label back to the term in `current`
    idx <- which(norm_term(current) == norm_term(worst$term))
    if (length(idx) != 1) {
      stop("Could not match '", worst$term, "' to a term in the model: ",
           paste(current, collapse = ", "))
    }
    
    history <- rbind(history, data.frame(step = nrow(history) + 1,
                                         removed = current[idx],
                                         R2 = worst$R2, p = worst$p))
    message("Step ", nrow(history), ": removed ", current[idx],
            " (R2 = ", round(worst$R2, 4), ", p = ", worst$p, ")")
    
    current <- current[-idx]
    if (length(current) == 0) break
  }
  
  final <- if (length(current) > 0) fit_margin(current) else NULL
  list(retained = current, removed = history, final_model = final)
}

res_back2 <- backward_permanova2(
  D            = dist_cc,
  data         = meta_cc,
  candidates   = candidates,                       # all covariates as main effects
  interactions = c("Time_num:Total_LIBRA_Score"),     # the only interaction tested
  alpha        = 0.1,
  permutations = 999,
  strata       = NULL
)

res_back2$retained
res_back2$removed
res_back2$final_model

final_terms <- res_back2$retained    # from the backward run

set.seed(123)
overall <- adonis2(
  reformulate(final_terms, response = "dist_cc"),
  data         = meta_cc,
  permutations = 999,
  by           = NULL      # overall model test (this is the default)
)

overall
overall$`Pr(>F)`[1]   # overall model p-value
overall$R2[1]         # total variance explained by the model
overall$F[1]          # overall pseudo-F





#beta long
beta_long = generate_beta_change_spaghettiplot_long(
  data.obj = data.obj,
  dist.obj = dist.obj,
  subject.var = "Participant",
  time.var = "Time_num", 
  t0.level = sort(unique(data.obj$meta.dat$Time_num))[1],
  ts.levels = sort(unique(data.obj$meta.dat$Time_num))[-1],
  group.var = "LIBRA_Half",
  strata.var = NULL,
  #dist.name = c("UniFrac","WUniFrac","Aitchison","phylo_rpca"),
  dist.name = c("phylo_rpca"),
  base.size = 20,
  theme.choice = "bw",
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11, 
  pdf.hei = 8.5
)


ggsave(paste0(path, "beta_long_plot_", Sys.Date(), ".jpeg"), beta_long, height = 10, width = 12, units="in",device="jpeg")


####plyo-rpca
rpca_data_long = data.frame(layer_data(beta_long$phylo_rpca, 1)) %>% select(group, colour, x,y) %>% mutate(LIBRA_Half = case_when(colour=="#E31A1C" ~ "Lower",
                                                                                                                                      T~"Upper"))

rpca_distance = ggplot(rpca_data_long, aes(x=x, y=y,group=LIBRA_Half))+
  geom_line(aes(group=group),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = rpca_data_long)+
  labs(tag = "A", x="Time (Months)", y="Phylo-RPCA Distance\n from Baseline", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c(  "#A77903", 
                 "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )


#'#0047AB' '#006666' '#002D72' '#FFB81C' "#1F3A5F"


lookup <- tibble(
  group = 1:17,
  Participant = c("P11","P14","P17","P19","P24","P25","P27","P28","P30",
                  "P31","P34","P42","P45","P48","P51","P52","P61")
)

rpca_data_long2 <- rpca_data_long %>% left_join(lookup, by = "group")
meta_and_alpha_libra_beta_long_filter = meta_and_alpha_libra_beta_long %>%
  select(Participant, Total_LIBRA_Score, age, sex, BMI_Current, Ever_Homeless) %>% distinct()
rpca_merge = left_join(rpca_data_long2, meta_and_alpha_libra_beta_long_filter, by="Participant")


#beta lmer
library(lme4)
library(lmerTest)
library(gamm4)
rpca_merge$Ever_Homeless = factor(rpca_merge$Ever_Homeless, levels=c("No", "Unknown", "Yes"))

rpca_model <- lmer(y ~ Total_LIBRA_Score*x+ age + as.factor(sex) + BMI_Current + factor(Ever_Homeless) + (1 | Participant), data = rpca_merge)
summary(rpca_model)
(anova(rpca_model))

rpca_model_gamm <- gamm4(
  y ~ t2(x, Total_LIBRA_Score) + age + as.factor(sex) + 
    BMI_Current + as.factor(Ever_Homeless),
  random = ~ (1 | Participant),
  data = rpca_merge
)

summary(rpca_model_gamm$gam)   # smooth-term significance, deviance explained, etc.
summary(rpca_model_gamm$mer)   # the underlying lmer fit, if you want the random-effect variance
plot(rpca_model_gamm$gam, pages = 1)  # visualize the fitted smooth surface



sup_table_uni_dist_model = flextable(tidy(uni_model))
sup_table_uni_dist_model = set_caption(sup_table_uni_dist_model, caption = "Supplementary Table x. Mixed Effects Model Unweighted Unifrac Distance") %>% autofit()
save_as_docx(sup_table_uni_dist_model, path = paste0(path, "Supplementary_Table_uni_mod._", Sys.Date(), ".docx"))


######Create figure 1 for all diversity 

#alpha spaghetti, alpha box, w coor, w distance, uw coord, uw distane
#coord_beta_w_fig1 = p_weight_group + labs(tag="C") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())
#coord_beta_uw_fig1 = p_unweight_group + labs(tag="E") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())
#dist_beta_w_fig1 = wuni_distance + labs(tag="D") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())
#dist_beta_uw_fig1 = uni_distance+ labs(tag="F") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())

coord_beta_aitchison_fig1 = p_aitch_group_plot + labs(tag="C") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())
dist_aitchison_fig1 = aitchison_distance + labs(tag="D") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())


coord_beta_rpca_fig1 = phylo_rpca + labs(tag="C") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())
dist_rpca__fig1 = rpca_distance + labs(tag="D") +  theme(plot.tag = element_text(size = 16, face = "bold"),plot.title = element_blank())

#beta_dista = ggarrange(wuni_distance,uni_distance, ncol=2,nrow=1)

figure_1 = ggarrange(shannon_line, coord_beta_rpca_fig1, shannon_libra_bin,dist_rpca__fig1, ncol=2,nrow=2)
ggsave(paste0(path, "figure_1_merged_", Sys.Date(), ".tiff"), figure_1, height = 10, width = 18, units="in",device="tiff")
ggsave(paste0(path, "beta_dist_figure_pres_merged_", Sys.Date(), ".jpeg"), beta_dista, height = 6, width = 12, units="in",device="jpeg")




########Differential abundance work
#maaslin
meta_and_libra_final = as.data.frame(meta_and_libra_final)
meta_and_libra_final$Ever_Homeless = factor(meta_and_libra_final$Ever_Homeless, levels=c("No", "Unknown", "Yes"))
rownames(meta_and_libra_final) = meta_and_libra_final$Sample_ID

ps_genus <- tax_glom(physeq_object_mindsource_final, "Genus")
#pull out taxa table into separate data frame
otu_ps <- otu_table(physeq_object_mindsource_final)
otu_ps <- t(otu_ps)
df_maslin <- as.data.frame(otu_ps)
colnames(df_maslin) <- as.data.frame(tax_table(physeq_object_mindsource_final))$Species



de_graph_fun = function(data_set, title, meta){
  de_data = data_set %>% 
    filter(qval_individual<0.2 & name == meta)%>%
    mutate(Up_Down = case_when(coef > 0 ~ "Up",
                               coef<  0  ~ "Down")) %>%
    arrange((qval_individual))%>% 
    mutate(feature = str_remove(feature, ".*?(?=s__)")) %>% 
    mutate(feature = str_replace(feature, "s__", "")) %>% 
    mutate(feature = str_replace(feature, "_", " "))
  
  de_plot = ggplot(de_data, aes(x=reorder(feature,coef), y=coef, fill=Up_Down)) +
    geom_bar(stat='identity') +
    theme_bw()+
    theme(legend.position = "none",
          text = element_text(size = 16, face="bold"))+
    #scale_x_discrete(position = "top")+
    theme(axis.text.y = element_text(face = "bold.italic")) + # Assign specific colors to groups
    coord_flip()+
    labs(#title="Species LFC by LIBRA Score", 
      x="Species", 
      y="Log Fold Change",
      title= title) +
    scale_fill_manual(values = c("Up" = "#1B3558", "Down" = "#A77903"))
  
  #if(length(de_data$name > 1)){
  # de_plot = de_plot + facet_wrap(~name, scales="free", nrow=3)
  #}else
  #{
  #  de_plot = de_plot
  #}
  
  return(de_plot)
}


physeq_object_mindsource_final_relative <- transform_sample_counts(physeq_object_mindsource_final, function(x) x / sum(x))

otu_ps <- otu_table(physeq_object_mindsource_final_relative)
otu_ps <- t(otu_ps)
df_maslin_rel <- as.data.frame(otu_ps)
colnames(df_maslin_rel) <- as.data.frame(tax_table(physeq_object_mindsource_final_relative))$Species

df_maslin_filter <- df_maslin_rel %>%
  select(where(~ !is.numeric(.x) || mean(.x >= 0.001, na.rm = TRUE) >= 0.1))

meta_and_libra_final$Education


library(maaslin3)
fit_data = maaslin3(
  input_data = df_maslin,
  input_metadata = meta_and_libra_final,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Projects/2026/Mindsource/Analysis/MaAsLin3/MaAsLin3_final_",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~Time_num*Total_LIBRA_Score + Age+Sex + group(Ever_Homeless) + BMI_Current+(1|Participant)",
  min_prevalence = 0.25,
  min_abundance = 10E-6,
  plot_associations = F,
  small_random_effects=F)



fit_data_filter = maaslin3(
  input_data = df_maslin_filter,
  input_metadata = meta_and_libra_final,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/Mindsource/Analysis/MaAsLin3/MaAsLin3_filtered_model_",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~Time_num*Total_LIBRA_Score + Age+Sex + group(Ever_Homeless) + BMI_Current+(1|Participant)",
  min_prevalence = 0,
  min_abundance = 0,
  plot_associations = F,
  small_random_effects=F)


fit_data_abund = maaslin3(
  input_data = taxa_table_all_filtered_mat,
  input_metadata = meta_and_libra_final,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("/path/to/output", Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~Time_num*LIBRA_Score_Total + Age+Sex+ any_homeless +(1|Participant)",
  min_prevalence = 0.25,
  min_abundance = 0,
  plot_associations = F,
  small_random_effects=F,                      
  plot_summary_plot = F,
  correction="BH",
  transform = "PLOG",
  evaluate_only = "abundance",
  warn_prevalence = F,
  zero_threshold = -1,
  median_comparison_abundance = F)


de_abund_full_time = de_graph_fun(fit_data$fit_data_abundance$results, "Differential Expression by Long Time for Mindsource, Abundance", "Time_num")
de_prevfull_time = de_graph_fun(fit_data$fit_data_prevalence$results, "Differential Expression by Long Time for Mindsource, Prevelence", "Time_num")
de_abund_full_libra = de_graph_fun(fit_data$fit_data_abundance$results, "Differential Expression by Long Time for Mindsource, Abundance", "Total_LIBRA_Score")
de_prevfull_libra = de_graph_fun(fit_data$fit_data_prevalence$results, "Differential Expression by Long Time for Mindsource, Prevelence", "Total_LIBRA_Score")
de_abund_full_time_libra = de_graph_fun(fit_data$fit_data_abundance$results, "Differential Expression of Time x LIBRA Score, Abundance", "Time_num:Total_LIBRA_Score")
de_prevfull_time_libra = de_graph_fun(fit_data$fit_data_prevalence$results, "Differential Expression of Time x LIBRA Score, Prevelence", "Time_num:Total_LIBRA_Score")


ggsave(paste0(path, "DE_long_mindsource_time_abundance_full_plot_", Sys.Date(), ".jpeg"), de_abund_full_time, height = 10, width = 10, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_time_prev_full_plot_", Sys.Date(), ".jpeg"), de_prevfull_time, height = 10, width = 10, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_libra_abundance_full_plot_", Sys.Date(), ".jpeg"), de_abund_full_libra, height = 10, width = 10, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_libra_prev_full_plot_", Sys.Date(), ".jpeg"), de_prevfull_libra, height = 10, width = 10, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_time_libra_abundance_full_plot_", Sys.Date(), ".jpeg"), de_abund_full_time_libra, height = 6, width = 6, units="in",device="jpeg")
ggsave(paste0(path, "DE_long_mindsource_time_libra_prev_full_plot_", Sys.Date(), ".jpeg"), de_prevfull_time_libra, height = 6, width = 6, units="in",device="jpeg")


abund_list = c("Clostridium clostridioforme", 'Akkermansia muciniphila', "Blautia hydrogenotrophica", "Sutterella wadsworthensis")

prev_list = c('Bilophila wadsworthia', 'Roseburia intestinalis', 'Roseburia inulinivorans', 'Eubacterium eligens')

scfa_list = c("Faecalibacterium prausnitzii", "Roseburia intestinalis", "Eubacterium rectale", "Eubacterium hallii",
              "Subdoligranulum variabile", "Roseburia hominis", "Anaerobutyricum soehngenii", "Coprococcus comes",
              "Coprococcus catus", "Butyricicoccus pullicaecorum", "Agathobacter rectalis", "Bacteroides thetaiotaomicron",
              "Bacteroides vulgatus", "Bacteroides uniformis", "Bacteroides fragilis", "Parabacteroides distasonis",
              "Megasphaera elsdenii", "Veillonella parvula", "Veillonella dispar", "Ruminococcus obeum", "Blautia obeum",
              "Bifidobacterium longum", "Bifidobacterium adolescentis", "Bifidobacterium breve", "Lactobacillus rhamnosus",
              "Lactobacillus plantarum", "Akkermansia muciniphila", "Ruminococcus gnavus")

# Create a regex pattern using the pipe '|' (which means "OR" in regex)
scfa_pattern <- paste(scfa_list, collapse = "|")

#Abundance plots
libra_filter_abd = fit_data$fit_data_abundance$results %>%
  filter(qval_individual<0.2 & name == "Time_num:Total_LIBRA_Score")%>%
  mutate(Up_Down = case_when(coef > 0 ~ "Up",
                             coef<  0  ~ "Down")) %>%
  arrange((qval_individual))%>%
  mutate(feature = str_remove(feature, ".*?(?=s__)")) %>%
  mutate(feature = str_replace(feature, "s__", "")) %>%
  mutate(feature = str_replace(feature, "_", " "))

libra_filter_prev = fit_data$fit_data_prevalence$results %>%
  filter(qval_individual<0.2 & name == "Time_num:Total_LIBRA_Score")%>%
  mutate(Up_Down = case_when(coef > 0 ~ "Up",
                             coef<  0  ~ "Down")) %>%
  arrange((qval_individual))%>%
  mutate(feature = str_remove(feature, ".*?(?=s__)")) %>%
  mutate(feature = str_replace(feature, "s__", "")) %>%
  mutate(feature = str_replace(feature, "_", " "))


sup_table_mal_abd = flextable(libra_filter_abd %>% dplyr::select(feature, name, coef, stderr,pval_individual, qval_individual))
sup_table_mal_abd = set_caption(sup_table_mal_abd, caption = "Supplementary Table x. Significant Species Abundance Associated with LIBRA x Time") %>% autofit()
save_as_docx(sup_table_mal_abd, path = paste0(path, "Supplementary_Table_libra_time_abnd._", Sys.Date(), ".docx"))


sup_table_mal_prev = flextable(libra_filter_prev %>% dplyr::select(feature, name, coef, stderr,pval_individual, qval_individual))
sup_table_mal_prev = set_caption(sup_table_mal_prev, caption = "Supplementary Table x. Significant Species Prevalence Associated with LIBRA x Time") %>% autofit()
save_as_docx(sup_table_mal_prev, path = paste0(path, "Supplementary_Table_libra_time_prev._", Sys.Date(), ".docx"))



libra_filter_prev$feature
#####creaate relative abundance 
physeq_object_mindsource_final_relative <- transform_sample_counts(physeq_object_mindsource_final, function(x) x / sum(x))
OTU_df <- as.data.frame(as(otu_table(physeq_object_mindsource_final_relative), "matrix"))


OTU_df$feature = rownames(OTU_df)
taxa_table_all_df2 = OTU_df %>%
  mutate(feature = str_remove(feature, ".*?(?=s__)")) %>%
  mutate(feature = str_replace(feature, "s__", "")) %>%
  mutate(feature = str_replace(feature, "_", " ")) %>%
  mutate(feature = str_remove(feature, "\\..*"))
#taxa_table_all_df_filter = taxa_table_all_df2 %>% filter(feature %in% libra_filter_abd$feature)
#taxa_table_all_df_filter = taxa_table_all_df2 %>% filter(str_detect(feature, scfa_pattern))


orgo_list = data.frame(taxa_table_all_df2$feature)
#rownames(taxa_table_all_df_filter) = taxa_table_all_df_filter$feature
relative_taxa_adbund_filtered_update2 <- data.frame(t(taxa_table_all_df2[-140]))
#colnames(relative_taxa_adbund_filtered_update2) <- taxa_table_all_df_filter[, 1]
relative_taxa_adbund_filtered_update3 <- relative_taxa_adbund_filtered_update2 %>%
  rownames_to_column(var = "Sample_ID")
relative_taxa_adbund_filtered_update3_melt = reshape2::melt(relative_taxa_adbund_filtered_update3, id.var="Sample_ID", variable.name="Species", value.name = "Abundance")
relative_taxa_adbund_filtered_update3_melt2 = relative_taxa_adbund_filtered_update3_melt %>%
  mutate(Species = str_remove(Species, ".*?(?=s__)")) %>%
  mutate(Species = str_replace(Species, "s__", "")) %>%
  mutate(Species = str_replace(Species, "_", " ")) %>%
  mutate(Species = str_remove(Species, "\\..*")) %>%
  mutate(Sample_ID = str_replace(Sample_ID, "X", "")) %>%
  dplyr::filter(Sample_ID %in% meta_and_libra_final$Sample_ID)

relative_taxa_adbund_filtered_update3_melt3 = relative_taxa_adbund_filtered_update3_melt2 %>%
  group_by(Sample_ID, Species) %>%
  mutate(total_abd = sum(as.numeric(Abundance))) %>%
  ungroup() %>%
  select(-Abundance) %>%
  distinct()
all_taxa_meta_data = merge(relative_taxa_adbund_filtered_update3_melt3, meta_and_libra_final, by="Sample_ID")
all_taxa_meta_data_update = all_taxa_meta_data %>% mutate(Prev = case_when(total_abd >0 ~ 1,
                                                                           T~0))

length(unique(meta_and_libra_final$Sample_ID))


sig_abnd_df = all_taxa_meta_data_update %>% dplyr::filter(Species %in% abund_list)

abund_fig = ggplot(sig_abnd_df, aes(y=total_abd,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, method="lm",data = sig_abnd_df)+
  labs(tag = "A", x="Time (Months)", y="Relative Abundance", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  facet_wrap(~Species, scales="free") +
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)

sig_abnd_df_supp = all_taxa_meta_data_update %>% dplyr::filter(Species %in% libra_filter_abd$feature)

abund_fig_supp = ggplot(sig_abnd_df_supp, aes(y=total_abd,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, method="lm",data = sig_abnd_df_supp)+
  labs( x="Time (Months)", y="Relative Abundance", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  facet_wrap(~Species, scales="free") +
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)


sig_prev_df = all_taxa_meta_data_update %>% dplyr::filter(Species %in% prev_list)

prev_plot = ggplot(sig_prev_df, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  #scale_size_continuous(limits=c(0, 1), breaks=seq(0, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  labs(tag = "B", x="Time (Months)", y="Participant", color="Brain Health Group \n(LIBRA Score)", size="Prevalence") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  scale_size(range = c(1, 3))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(group=LIBRA_Half)) +
  theme(strip.text.x = element_text(face = "italic"))+ # For x-axis facet titles)
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )


sig_prev_df_supp = all_taxa_meta_data_update %>% dplyr::filter(Species %in% libra_filter_prev$feature)

prev_plot_supp = ggplot(sig_prev_df_supp, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  #scale_size_continuous(limits=c(0, 1), breaks=seq(0, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  labs( x="Time (Months)", y="Participant", color="Brain Health Group \n(LIBRA Score)", size="Prevalence") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  scale_size(range = c(1, 3))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(group=LIBRA_Half)) +
  theme(strip.text.x = element_text(face = "italic"))+ # For x-axis facet titles)
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )



####Figure 3
de_abund_fig3 = de_abund_full_time_libra + labs(tag="A") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank(),    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
de_prev_fig3 = de_prevfull_time_libra + labs(tag="C") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank(),    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
species_abund_fig3 = abund_fig + labs(tag="B") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
prev_abund_fig3 = prev_plot+ labs(tag="D") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
ggsave(paste0(path, "figure_3b_abund_", Sys.Date(), ".tiff"), abund_fig, height = 4, width = 6, units="in",device="jpeg")
ggsave(paste0(path, "figure_3d_prev_", Sys.Date(), ".tiff"), prev_plot, height = 6, width = 10, units="in",device="jpeg")


species_abund_supp = abund_fig_supp +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
prev_abund_supp = prev_plot_supp +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
ggsave(paste0(path, "figure_supp_abund_", Sys.Date(), ".tiff"), species_abund_supp, height = 10, width = 20, units="in",device="tiff")
ggsave(paste0(path, "figure_supp_prev_", Sys.Date(), ".tiff"), prev_abund_supp, height = 12, width = 18, units="in",device="tiff")


figure_3 = ggarrange(de_abund_fig3,species_abund_fig3, de_prev_fig3, prev_abund_fig3, ncol=2,nrow=2)
ggsave(paste0(path, "figure_3_merged_", Sys.Date(), ".tiff"), figure_3, height = 14, width = 18, units="in",device="jpeg")


####supplemental figure for time and libra seperate
####Figure 3
de_abund_time_supp = de_abund_full_time + labs(tag="A") +  theme(aspect.ratio= 1, plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank()) +
  scale_x_discrete(labels = label_wrap(15)) # Wraps text after ~15 characters1)
de_prev_time_supp = de_prevfull_time + labs(tag="B") +  theme(aspect.ratio= 1,plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())+
  scale_x_discrete(labels = label_wrap(15)) # Wraps text after ~15 characters1)
de_abund_libra_supp = de_abund_full_libra + labs(tag="C") +  theme(aspect.ratio= 1,plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())+
  scale_x_discrete(labels = label_wrap(15)) # Wraps text after ~15 characters1)
de_prev_libra_supp = de_prevfull_libra+ labs(tag="D") +  theme(aspect.ratio= 1,plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())+
  scale_x_discrete(labels = label_wrap(15)) # Wraps text after ~15 characters1)



figure_supp_time_and_libra = ggarrange(de_abund_time_supp, de_prev_time_supp,de_abund_libra_supp, de_prev_libra_supp, ncol=2,nrow=2)
ggsave(paste0(path, "figure_supplemental_time_and_libra_merged_", Sys.Date(), ".tiff"), figure_supp_time_and_libra, height = 20, width = 20, units="in",device="tiff")




#####Start of SCFA Bacteria abundance figures
all_taxa_meta_data_update_scfa = all_taxa_meta_data_update %>% filter(Species %in% scfa_list) %>% group_by(Participant, Time_num) %>% mutate(scfa_abund = sum(total_abd))

all_taxa_meta_data_update_scfa_2 = all_taxa_meta_data_update_scfa %>% select(Participant, 
                                                                             age, 
                                                                             sex, 
                                                                             Ever_Homeless, 
                                                                             Time_num, 
                                                                             scfa_abund, 
                                                                             LIBRA_Half,
                                                                             Total_LIBRA_Score, 
                                                                             BMI_Current, 
                                                                             Education) %>% distinct()

# 1. Turn NA into a valid factor level (it defaults to literal "<NA>")
all_taxa_meta_data_update_scfa_2$Education <- addNA(all_taxa_meta_data_update_scfa_2$Education)

# 2. Rename that specific "<NA>" level to "unknown"
levels(all_taxa_meta_data_update_scfa_2$Education)[is.na(levels(all_taxa_meta_data_update_scfa_2$Education))] <- "Unknown"



all_taxa_meta_data_update_scfa_2 %>%
  filter(Time_num == 12 | Time_num == 11) %>%
  group_by(LIBRA_Half) %>%
  summarise(sd_scfa_abund = sd(scfa_abund, na.rm = TRUE))

scfa_abund_fig = ggplot(all_taxa_meta_data_update_scfa_2, aes(y=scfa_abund,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  geom_line(aes(group=Participant),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  #geom_smooth(aes(color=LIBRA_Half),se=F, method="lm",data = all_taxa_meta_data_update)+
  geom_smooth(aes(color=LIBRA_Half),se=F,data = all_taxa_meta_data_update_scfa_2)+
  labs(tag = "A", x="Time (Months)", y="Summed Relative Abundance of \nSCFA Associated Bacteria", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #facet_wrap(~Species, scales="free") +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
ggsave(paste0(path, "scfa_spaghetti_plot_", Sys.Date(), ".jpeg"), scfa_abund_fig, height = 8, width = 15, units="in",device="jpeg")

scfa_abund_fig2 = ggplot(all_taxa_meta_data_update_scfa_2, aes(y=scfa_abund,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  geom_line(aes(group=Participant,color=LIBRA_Half))+
  #geom_point()+
  #facet_wrap(~Participant)+
  #geom_smooth(aes(color=LIBRA_Half),se=F, method="lm",data = all_taxa_meta_data_update)+
  geom_smooth(aes(color=LIBRA_Half),se=F)+
  labs(tag = "A", x="Time (Months)", y="Summed Relative Abundance of \nSCFA Associated Bacteria", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 14, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  #facet_wrap(~Species, scales="free") +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)





scfa_abund_bin = ggplot(all_taxa_meta_data_update_scfa_2, aes(x=reorder(Participant, Total_LIBRA_Score), y=scfa_abund, color=LIBRA_Half))+
  geom_boxplot()+
  geom_point()+
  labs(tag = "B", x="Participants \n(Ordered by LIBRA Score)", y="Summed Relative Abundance of \nSCFA Associated Bacteria", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  #scale_x_reverse()+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )

ggsave(paste0(path, "scfa_box_plot_", Sys.Date(), ".jpeg"), scfa_abund_bin, height = 8, width = 15, units="in",device="jpeg")

figure_2= ggarrange(scfa_abund_fig, scfa_abund_bin)
ggsave(paste0(path, "figure_2_merged_", Sys.Date(), ".tiff"), figure_2, height = 6, width = 14, units="in",device="tiff")

library(lme4)
library(lmerTest)
library(car)
model_score_scfa <- lmer(scfa_abund ~ Time_num*Total_LIBRA_Score + 
                           age+
                           as.factor(sex) +
                           BMI_Current+ 
                           as.factor(Ever_Homeless) +
                           (1 | Participant), data = all_taxa_meta_data_update_scfa_2)
summary(model_score_scfa)
anova(model_score_scfa)

chi_test <- chisq.test(all_taxa_meta_data_update_scfa_2$Education, all_taxa_meta_data_update_scfa_2$Ever_Homeless)
print(chi_test)

# Install package if needed: install.packages("vcd")
library(vcd)

# Generate a contingency table and calculate association stats
assocstats(table(all_taxa_meta_data_update_scfa_2$Education, all_taxa_meta_data_update_scfa_2$Ever_Homeless))


all_taxa_meta_data_update_scfa_2$Education <- droplevels(all_taxa_meta_data_update_scfa_2$Education)

model_score_scfa_gamm <- gamm4(
  scfa_abund ~ t2(Time_num,Total_LIBRA_Score) + age + as.factor(sex) + 
     as.factor(Ever_Homeless),
  random = ~ (1 | Participant),
  data = all_taxa_meta_data_update_scfa_2
)

summary(model_score_scfa_gamm$gam)   # smooth-term significance, deviance explained, etc.
summary(model_score_scfa_gamm$mer)   # the underlying lmer fit, if you want the random-effect variance
plot(model_score_scfa_gamm$gam, pages = 1)  # visualize the fitted smooth surface




all_taxa_meta_data_update_scfa_2_Test = all_taxa_meta_data_update_scfa_2 %>% dplyr::filter(!(Participant=="P24" & Time_num == 12))
model_score_scfa <- lmer(scfa_abund ~ Time_num*Total_LIBRA_Score + age+sex +BMI_Current+ Ever_Homeless+ (1 | Participant), data = all_taxa_meta_data_update_scfa_2)
summary(model_score_scfa)
anova(model_score_scfa)


sup_table_scfab_model = flextable(tidy(model_score_scfa))
sup_table_scfab_model = set_caption(sup_table_scfab_model, caption = "Supplementary Table x. Mixed Effects Model SCFA-Bacteria Abundance Association LIBRA x Time") %>% autofit()
save_as_docx(sup_table_scfab_model, path = paste0(path, "Supplementary_Table_wuni_mod._", Sys.Date(), ".docx"))




model_half_scfa <- lmer(scfa_abund ~ Time_num*LIBRA_Half + age+sex + Ever_Homeless+ (1 | Participant), data = all_taxa_meta_data_update_scfa_2)
summary(model_half_scfa)
modelsummary(model_score_scfa, output = paste0(path, "scfa_model_summary.docx"))


cv_by_libra_time <- all_taxa_meta_data_update_scfa_2 %>%
  ungroup() %>%
  group_by(Time_num) %>%
  summarize(
    Mean_Value = mean(scfa_abund, na.rm = TRUE),
    SD_Value = sd(scfa_abund, na.rm = TRUE),
    CV_Percentage = calculate_cv(scfa_abund)
  )

cv_by_libra_time  %>% summarize(mean(CV_Percentage))

cv_by_libra_half <- all_taxa_meta_data_update_scfa_2 %>%
  group_by(LIBRA_Half) %>%
  summarize(
    Mean_Value = mean(scfa_abund, na.rm = TRUE),
    SD_Value = sd(scfa_abund, na.rm = TRUE),
    CV_Percentage = calculate_cv(scfa_abund)
  )


all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 0| Time_num == 0.07 )) %>% # & LIBRA_Half == "Lower") %>% 
  summarise(value = mean(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 0| Time_num == 0.07 ) ) %>% #& LIBRA_Half == "Lower") %>% 
  summarise(value = sd(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 0 | Time_num == 0.07)  & LIBRA_Half == "Upper") %>% 
  summarise(value = mean(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 0 | Time_num == 0.07)  & LIBRA_Half == "Upper") %>% 
  summarise(value = sd(scfa_abund))


all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 12 |Time_num == 11))%>% #  & LIBRA_Half == "Lower") %>% 
  summarise(value = mean(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 12 |Time_num == 11)) %>% # & LIBRA_Half == "Lower") %>% 
  summarise(value = sd(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Upper") %>% 
  summarise(value = mean(scfa_abund))

all_taxa_meta_data_update_scfa_2 %>% ungroup() %>%
  filter((Time_num == 12 |Time_num == 11)  & LIBRA_Half == "Upper") %>% 
  summarise(value = sd(scfa_abund))




abund_fig_test = ggplot(all_taxa_meta_data_update_scfa_2 , aes(y=scfa_abund,x=BMI_Current))+
  #geom_line(aes(group=Participant),color="grey")+
  geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(se=F, method="lm",data = all_taxa_meta_data_update_scfa_2) +#%>% filter(Time_num < 1  ))+
  labs(tag = "A", x="LIBRA Score", y="Relative Abundance")+#, color="Brain Health Group \n(LIBRA Score)") +
  theme(plot.tag = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  stat_cor(label.y.npc="top", label.x = 30, size = 3)

#facet_wrap(~Species, scales="free") +
#scale_colour_manual(
# values = c( "#4F94CD",  "#CD3333"),
#labels= c("Better Score", "Worse Score")
#)+
#theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)


all_taxa_meta_data_update_scfa_2_test = all_taxa_meta_data_update_scfa_2 %>% filter(Time_num > 6 & Total_LIBRA_Score < 10)
model_score_scfa_test <- lm(scfa_abund ~ Total_LIBRA_Score + age+sex +BMI_Current+ Ever_Homeless, data = all_taxa_meta_data_update_scfa_2)
summary(model_score_scfa_test)
(anova(model_score_scfa_test))



#####Pathways
path_table_1
path_table_1 <- path_table_1[-c(1, 2), ]
path_table_2
path_table_2 <- path_table_2[-c(1, 2), ]
path_table_all = merge(path_table_1,path_table_2, by=0, all=T)
path_table_all[is.na(path_table_all)] <- 0
rownames(path_table_all) = path_table_all$Row.names
path_table_all = path_table_all %>% dplyr::select(-Row.names)
# Remove the first two rows
#taxa_table_all <- taxa_table_all[-c(1, 2), ]
path_table_all_df = data.frame(path_table_all)
names(path_table_all_df) <- gsub("X", "", names(path_table_all_df))
path_table_all_df_filter <- path_table_all_df[, names(path_table_all_df) %in% meta_and_libra_final$Sample_ID]
#meta_and_alpha %>% filter(SampleID == "19-1307.P24.8m-g")
#dat[-1] <-sweep(dat[-1], 2, colSums(dat[,-1]), `/`) * 10000
#Round data so that all the data are integers

library(edgeR)
library(DESeq2)
library(splines)




pathway_df = data.frame(path_table_all_df_filter)
dat_pathway <- pathway_df %>%
  dplyr::mutate_if(is.numeric, round)
names(dat_pathway) <- gsub("X", "", names(dat_pathway))
keep_genes <- colSums( dat_pathway > 100 ) >= 70
dat_pathway <- dat_pathway[keep_genes,]
dat_pathway_2 <- as.matrix(dat_pathway)
library(maaslin3)
#' @return `x` with only the rows that pass the filter.
filter_rows <- function(x, min_value = 0.25, min_prop_samples = 0.5) {
  m <- as.matrix(x)
  if (!is.numeric(m)) stop("`x` must contain only numeric values.")
  
  n_samples <- ncol(m)
  # Small tolerance guards against floating-point error (e.g., 0.3 * 10 = 3.0000000000000004)
  min_samples <- ceiling(min_prop_samples * n_samples - 1e-9)
  
  n_passing <- rowSums(m >= min_value, na.rm = TRUE)
  keep <- n_passing >= min_samples
  
  message("Requiring >= ", min_value, " in at least ", min_samples, " of ",
          n_samples, " samples: ", sum(keep), " of ", nrow(m), " rows kept.")
  x[keep, , drop = FALSE]
}

run_maaslin <- function(df, meta, output, formula, min_prevalence, min_abundance,
                        plot_associations = FALSE, ...) {
  maaslin3(
    input_data = df, input_metadata = meta, output = output, formula = formula,
    min_prevalence = min_prevalence, min_abundance = min_abundance,
    plot_associations = plot_associations, ...
  )
}

r <- filter_rows(dat_pathway_2, min_value = 100, min_prop_samples = 0.5)

dat_pathway_filter <- filter_rows(dat_pathway_2, min_value = 100, min_prop_samples = 0.1)
dat_pathway_filter2 <- filter_rows(dat_pathway_2, min_value = 1000, min_prop_samples = 0.25)


fit_con_path <- run_maaslin(dat_pathway_filter, meta_and_libra_final, paste0(path,"MaAsLin3_path_confoudners_", Sys.Date()),
                            "~Time_num*Total_LIBRA_Score + Age+Sex + group(Ever_Homeless) + BMI_Current+(1|Participant)", min_prevalence = 0, min_abundance = 0,max_significance=0.2)

fit_con_path2 <- run_maaslin(dat_pathway_filter2, meta_and_libra_final, paste0(path,"MaAsLin3_path_confoudners2_", Sys.Date()),
                            "~Time_num*Total_LIBRA_Score + Age+Sex + group(Ever_Homeless) + BMI_Current+(1|Participant)", min_prevalence = 0, min_abundance = 0,max_significance=0.2)


set.seed(123)

d0 <- DGEList(dat_pathway_2)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
meta_and_libra_final

X <- ns(meta_and_libra_final$Time_num, df=3)
X2 = ns(meta_and_libra_final$Total_LIBRA_Score, df=3)
design <- model.matrix(~Time_num*Total_LIBRA_Score+age+as.factor(sex)+BMI_Current + as.factor(Ever_Homeless), (meta_and_libra_final))

y <- voom(d, design, plot = T)

y <- voom(d, design)
corfit <- duplicateCorrelation(y, design, block = meta_and_libra_final$Participant)
y <- voom(d, design, block = meta_and_libra_final$Participant, correlation =
            corfit$consensus, plot=T)

fit <- lmFit(y, design,block = meta_and_libra_final$Participant, correlation =
               corfit$consensus )
fit <- eBayes(fit)
top.table = topTable(fit, coef= "Time_num:Total_LIBRA_Score", n=Inf)
top.table2 = top.table %>% filter(adj.P.Val<0.2)
which(top.table$adj.P.Val < 0.2)


#This will be empty if no genes were significant
signif_genes <- as.character(rownames(top.table2))
#Below is old code to export sig-taxa
sigtab = cbind(as(top.table, "data.frame"), all_genes)


##Function to create volcano plot
volcplot <- function(data, padj_threshold = 0.2, fc = 1, plot_title = 'Volcano Plot', plot_subtitle = NULL, genelist_vector = genes$all_genes, genelist_filter = FALSE) {
  # Set the fold-change thresholds
  neg_log2fc <- -log2(fc)
  pos_log2fc <- log2(fc)
  # Make a dataset for plotting, add the status as a new column
  plot_ready_data <- data %>%
    mutate_at('adj.P.Val', ~replace(.x, is.na(.x), 1)) %>%
    mutate_at('logFC', ~replace(.x, is.na(.x), 0)) %>%
    mutate(
      log2fc_threshold = ifelse(logFC >= pos_log2fc & adj.P.Val <= padj_threshold, 'up',
                                ifelse(logFC <= neg_log2fc & adj.P.Val <= padj_threshold, 'down', 'ns')
      )
    ) %>%
    mutate(all_genes = replace_na(all_genes, 'none'))
  if (genelist_filter) {
    plot_ready_data <- plot_ready_data %>% dplyr::filter(all_genes %in% genelist_vector)
  }
  if(!is.null(genelist_vector)) {
    plot_ready_data <- plot_ready_data %>% mutate(all_genes = ifelse(all_genes %in% genelist_vector & adj.P.Val <= padj_threshold & log2fc_threshold != 'ns', all_genes, ''))
  }
  # Get the number of up, down, and unchanged genes
  up_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'up') %>% nrow()
  down_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'down') %>% nrow()
  unchanged_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'ns') %>% nrow()
  # Make the labels for the legend
  legend_labels <- c(
    str_c('Up: ', up_genes),
    str_c('NS: ', unchanged_genes),
    str_c('Down: ', down_genes)
  )
  # Set the x axis limits, rounded to the next even number
  x_axis_limits <- DescTools::RoundTo(
    max(abs(plot_ready_data$logFC)),
    0.1,
    ceiling
  )
  
  # Set the plot colors
  plot_colors <- c(
    'up' = "#1B3558",
    'ns' = 'gray',
    'down' = "#A77903"
  )
  # Make the plot, these options are a reasonable strting point
  plot <- ggplot(plot_ready_data) +
    geom_point(
      alpha = 0.25,
      size = 1.5
    ) +
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      color = log2fc_threshold,
      label = all_genes
    ) +
    geom_vline(
      xintercept = c(neg_log2fc, pos_log2fc),
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = -log10(padj_threshold),
      linetype = 'dashed'
    ) +
    scale_x_continuous(
      'log2(FC)',
      limits = c(-x_axis_limits, x_axis_limits)
    ) +
    scale_y_continuous('-log10(padj)',
                       limits = c(-0.25,1.5)
    )+
    scale_color_manual(
      values = plot_colors#,
      #labels = legend_labels
    ) +
    #labs(
    #  color = str_c(fc, '-fold, padj ≤', padj_threshold),
    #  title = plot_title,
    #  subtitle = plot_subtitle
    #) +
    theme_bw(base_size = 24) +
    theme(text = element_text(size = 14, face="bold"),
          legend.position = "none",
          #aspect.ratio = 1,
          axis.text = element_text(color = 'black'),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),  # Reduces dead area around legend
          legend.spacing.x = unit(0.2, 'cm'),
          legend.title = element_blank()
    )
  # Add gene labels if needed
  if (!is.null(genelist_vector)) {
    plot <- plot +
      geom_label_repel(
        size = 3.5,
        force = 1,
        max.overlaps = 100000,
        nudge_x = 0,
        nudge_y = 0,
        segment.color = 'black',
        min.segment.length = 0,
        show.legend = FALSE
      )
  }
  plot
}

top.table

sigtab2 = na.omit(sigtab)
rownames(sigtab2) <- NULL
#pulls out genes for labeling
genes = sigtab2 %>% dplyr::filter(logFC > log2(1) | logFC< -log2(1)) %>%
  dplyr::filter(adj.P.Val < 0.2 & P.Value < 0.05) %>%
  arrange(desc(abs(logFC))) %>% slice_head(n=6)
out_table = sigtab2 %>% dplyr::filter(adj.P.Val < 0.2 & P.Value < 0.05) %>%
  arrange(desc(logFC))
#creates volcano plot and prints
library(ggrepel)
sigtab2_udpate = sigtab2 %>% rename(padj = adj.P.Val, log2FoldChange = logFC)
volc = volcplot(sigtab2)
volc


ggsave(paste0(path, "figure_4_volc_", Sys.Date(), ".jpeg"), volc, height = 7.5, width = 7.5, units="in",device="jpeg")

#volc = volc + scale_y_continuous(limits=c(0))
out_table_pathways = genes %>% dplyr::select(all_genes, logFC, P.Value, adj.P.Val) %>% arrange(desc(logFC))
out_table_pathways

out_table_pathways %>% filter(log2FoldChange < 0)
genes
sup_table_3 = flextable(out_table_pathways)
sup_table_3 = set_caption(sup_table_3, caption = "Supplementary Table 3. Metabolic Pathways associated with LIBRA Scores") %>% autofit()
save_as_docx(sup_table_3, path = paste0(path, "Supplementary_Table_3._", Sys.Date(), ".docx"))

##############
#####Start EC analysis
#############

ec_table_1

ec_table_2 

ec_table_1 <- ec_table_1[-c(1, 2), ]
ec_table_2 <- ec_table_2[-c(1, 2), ]
ec_table_all = merge(ec_table_1,ec_table_2, by=0, all=T)
ec_table_all[is.na(ec_table_all)] <- 0
rownames(ec_table_all) = ec_table_all$Row.names
ec_table_all = ec_table_all %>% dplyr::select(-Row.names)
# Remove the first two rows
#taxa_table_all <- taxa_table_all[-c(1, 2), ]
ec_table_all_df = data.frame(ec_table_all)
names(ec_table_all_df) <- gsub("X", "", names(ec_table_all_df))
ec_table_all_df_filter <- ec_table_all_df[, names(ec_table_all_df) %in% meta_and_libra_final$Sample_ID]
#meta_and_alpha %>% filter(SampleID == "19-1307.P24.8m-g")
#dat[-1] <-sweep(dat[-1], 2, colSums(dat[,-1]), `/`) * 10000
#Round data so that all the data are integers

ec_df = data.frame(ec_table_all_df_filter)
dat_ec <- ec_df %>%
  dplyr::mutate_if(is.numeric, round)
names(dat_ec) <- gsub("X", "", names(dat_ec))
keep_genes <- colSums( dat_ec > 100 ) >= 70
dat_ec <- dat_ec[keep_genes,]
dat_ec_2 <- as.matrix(dat_ec)


dat_ec_2



trans_dat_ec = t(data.frame(dat_ec))


#trans_dat_pathway = trans_dat_pathway %>% row_to_names(row_number = 1)
trans_dat_ec = data.frame(trans_dat_ec)
trans_dat_ec$Sample_ID = rownames(trans_dat_ec)


trans_dat_ec_ID = trans_dat_ec %>% mutate(Sample_ID = str_replace(Sample_ID, "X", ""))
rownames(trans_dat_ec_ID) = trans_dat_ec_ID$Sample_ID
trans_ec_update = trans_dat_ec_ID %>% select(-Sample_ID) %>% select(where(~ any(. != 0)))
names(trans_ec_update) <- gsub("X", "", names(trans_ec_update))



filter_ec_table 
trans_ec_update_filter = trans_ec_update %>% filter()
trans_ec_update_filter <- trans_ec_update %>% 
  select(any_of(filter_ec_table$EC_Number))


fit_data_ec = maaslin3(
  input_data = trans_ec_update_filter,
  input_metadata = meta_and_libra_final,
  #output = paste0(new_dirc, "/", name, "/", name, "_MaAsLin3_output"),
  output= paste0("Path/SCFA/",Sys.Date()),
  #formula = "~group_cat+age+group(race_cat)+group(Picogreen_cat)+BMI",
  formula = "~Time_num*Total_LIBRA_Score + Age+Sex + group(Ever_Homeless) + BMI_Current+(1|Participant)",
  min_prevalence = 0,
  min_abundance = 0,
  plot_associations = F,
  small_random_effects=F,
  max_significance = 0.2)



trans_ec_update_filter2 <- trans_ec_update_filter %>%
  rownames_to_column(var = "Sample_ID")
trans_ec_update_filter2_melt = reshape2::melt(trans_ec_update_filter2, id.var="Sample_ID", variable.name="EC", value.name = "Abundance")
trans_ec_update_filter2_melt2 = trans_ec_update_filter2_melt %>%
  mutate(Sample_ID = str_replace(Sample_ID, "X", "")) %>%
  dplyr::filter(Sample_ID %in% meta_and_libra_final$Sample_ID)

trans_ec_update_filter2_melt3 = trans_ec_update_filter2_melt2 %>%
  group_by(Sample_ID, EC) %>%
  mutate(total_abd = sum(as.numeric(Abundance))) %>%
  ungroup() %>%
  select(-Abundance) %>%
  distinct()
all_taxa_meta_data_ec = merge(trans_ec_update_filter2_melt3, meta_and_libra_final, by="Sample_ID")
all_taxa_meta_data_ec_update = all_taxa_meta_data_ec %>% mutate(Prev = case_when(total_abd >0 ~ 1,
                                                                           T~0))




ec_abund_bin = ggplot(all_taxa_meta_data_ec_update, aes(x=reorder(Participant, Total_LIBRA_Score), y=total_abd, color=LIBRA_Half))+
  geom_boxplot()+
  geom_point()+
  labs(tag = "B", x="Participants \n(Ordered by LIBRA Score)", y="Summed Relative Abundance of \nSCFA Associated Bacteria", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~EC)+
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )




ec_fig_supp = ggplot(all_taxa_meta_data_ec_update %>% filter(EC == "2.7.2.7"), aes(y=total_abd,x=Time_num,group=LIBRA_Half,color=LIBRA_Half))+
  #geom_line(aes(group=Participant,color=LIBRA_Half))+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F)+
  labs( x="Time (Months)", y="Relative Abundance", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 16, face = "bold"))+
  facet_wrap(~Participant, scales="free") +
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c( "#A77903", 
                "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )+
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)

ec_csv

fit_data_ec2 = fit_data_ec
filter_ec_table$Pathway...Reaction
ec_csv_abund$feature <- filter_ec_table$Enzyme[match(ec_csv_abund$feature, filter_ec_table$EC_Number)]
ec_csv_abund$feature <- paste0(ec_csv_abund$feature, ": ", filter_ec_table$EC_Number[match(ec_csv_abund$feature, filter_ec_table$Enzyme)])


fit_data_ec2$fit_data_abundance$results$feature <- filter_ec_table$Enzyme[match(fit_data_ec2$fit_data_abundance$results$feature, filter_ec_table$EC_Number)]
fit_data_ec2$fit_data_abundance$results$feature <- paste0(fit_data_ec2$fit_data_abundance$results$feature, ": ", filter_ec_table$EC_Number[match(fit_data_ec2$fit_data_abundance$results$feature, filter_ec_table$Enzyme)])


de_abund_full_libra = de_graph_fun(fit_data_ec2$fit_data_abundance$results, "Differential Expression of EC Numbers", "Time_num:Total_LIBRA_Score")
de_abund_full_libra2 = de_graph_fun(ec_csv_abund, "Differential Expression of EC Numbers", "Time_num:Total_LIBRA_Score")

de_abund_full_libra2 = de_abund_full_libra +
  labs(x="Enzyme Commission Number")

de_abund_full_libra3 = de_abund_full_libra2 +
  labs(x="Enzyme Commission Number")


library(scales)
####Figure 4
pathway_volcano = volc + labs(tag="B") +  theme(aspect.ratio=12/16,
                                                plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
de_ec_fig = de_abund_full_libra2 + labs(tag="A") +  
  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())+
  scale_x_discrete(labels = label_wrap(15)) # Wraps text after ~15 characters1)

de_ec_fig = de_abund_full_libra3 + labs(tag="A") +  
  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())+
  scale_x_discrete(labels = label_wrap(15)) # Wraps text after ~15 characters1)


pathway_fig_merged = ggarrange(de_ec_fig,pathway_volcano, ncol=1,nrow=2) + 
  theme(plot.margin = margin(t = 0, r = 10, b = 0, l = 0, "mm"))
ggsave(paste0(path, "figure_4_merged_", Sys.Date(), ".tiff"), pathway_fig_merged, height = 12, width = 10, units="in",device="tiff")


# Optional: Arrange with the highest counts at the top
meta_and_libra_final %>% dplyr::count(Participant, sort = TRUE)



?maaslin3
#MaAslin Extra ###################################################
## Supplemental#################
######Extra for MaAsLin
#MaAsLin analysis for each libra factor
mal_fun= function(effect,ref){
  fit = Maaslin2(input_data = input_data, input_metadata = meta, min_prevalence = .5, 
                 min_abundance = .01, normalization = "TSS", 
                 output = paste0("MaAsLin_Out_species-level_",effect,"_",Sys.Date()), 
                 fixed_effects = effect,
                 reference = ref,
                 plot_scatter = F,
                 plot_heatmap = F,                          
                 correction="BH",
                 analysis_method = "LM",
                 transform = "LOG",
                 max_significance = 0.3)#,
  #reference = "Current_alcohol,No")}
  return(fit)
}


#


meta$Physical.inactivity_lib.weight = as.factor(meta$Physical.inactivity_lib.weight)
#head(meta$Healthy_diet)
depressed_mal = mal_fun("Depressed", "Depressed,No")

obese_mal = mal_fun("Obese", "Obese,No")

alc_mal = mal_fun("Current_alcohol","Current_alcohol,No")

chd_mal = mal_fun("Coronary.heart.disease", "Coronary.heart.disease,No")

cdd_mal = mal_fun("Chronic_kidney_disease", "Chronic_kidney_disease,No")

diabetes_mal = mal_fun("Diabetes","Diabetes,No")

Hypercholesterolemia_mal = mal_fun("Hypercholesterolemia", "Hypercholesterolemia,No")

Lifetime_smoking_mal = mal_fun("Lifetime_smoking", "Lifetime_smoking,No")

Hypertension_mal = mal_fun("Hypertension","Hypertension,No")

Physical.inactivity_lib.weight_mal = mal_fun("Physical.inactivity_lib.weight", "Physical.inactivity_lib.weight,0")

Healthy_diet_mal = mal_fun("Healthy_diet", "Healthy_diet,Yes")




abundance_plot = ggplot(all_taxa_meta_data_update, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  scale_size_continuous(limits=c(1, 1), breaks=seq(1, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(color=Total_LIBRA_Score)) +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
abundance_plot
all_taxa_meta_data_update

abundance_plot = ggplot(all_taxa_meta_data_update, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  scale_size_continuous(limits=c(1, 1), breaks=seq(1, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(color=Total_LIBRA_Score)) +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
abundance_plot

abundance_plot = ggplot(all_taxa_meta_data_update, aes(y=reorder(Participant, Total_LIBRA_Score), x=Time_num,color=LIBRA_Half)) +
  geom_point(aes(size=Prev)) +
  scale_size_continuous(limits=c(0, 1), breaks=seq(0, 1, by=1))+
  facet_wrap(~Species, scales="free") +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(color=Total_LIBRA_Score)) +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
abundance_plot

prev_plot = ggplot(all_taxa_meta_data_update, aes(y=total_abd, x=Time_num,group=LIBRA_Half, color=LIBRA_Half)) +
  #geom_point(aes(size=Prev)) +
  stat_summary(fun = "mean", geom = "point")+
  #scale_size_continuous(limits=c(10E-10, 0.2), breaks=seq(10E-10, 0.2, by=0.01))+
  facet_wrap(~Species, scales="free") +
  #geom_smooth(method="lm",  aes(color=LIBRA_Half)) +
  #geom_smooth(method="glm", method.args = list(family = "binomial"), aes(group=LIBRA_Half)) +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
prev_plot

#Function to create beta diversity measure and figure from phyloseq object
#inputs include phyloseq object, distance measure, weighted true or false, the coloring group, title for plot, outpath, and the adonis p-value
weighted_fun = function (data_set,da_method, weight_T_F, color_group, title_name,path_out){
  weight_ordination = ordinate(data_set, da_method, "unifrac", weighted=weight_T_F)
  p_weight <- plot_ordination(data_set,weight_ordination, type="samples", color=color_group, title=title_name)+ theme_bw()
  p_weight2 = p_weight + 
    stat_ellipse(type = "norm", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 16, face="bold"))+
    labs(color="Brain Health Group \n(LIBRA Score)") +
    #theme(plot.tag.position = c(.885, .98), plot.tag = element_text(size=7.25))+
    scale_colour_manual(
      values = c( "#A77903", 
                  "#1B3558"),
      labels= c("Lower Risk Group", "Higher Risk Group")
    )
  ggsave(path_out, p_weight2,device="jpeg", height = 6, width = 8, units = "in")
  return(p_weight2)
}

?ordinate
weighted_fun_part = function (data_set,da_method, weight_T_F, color_group, title_name,path_out){
  weight_ordination = ordinate(data_set, da_method, "unifrac", weighted=weight_T_F)
  p_weight <- plot_ordination(data_set,weight_ordination, type="samples", color=color_group, title=title_name)+ theme_bw()
  p_weight2 = p_weight + 
    stat_ellipse(type = "norm", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 12, face="bold"))+
    labs(color="Participant")
  #theme(plot.tag.position = c(.885, .98), plot.tag = element_text(size=7.25))+
  ggsave(path_out, p_weight2,device="jpeg", height = 6, width = 8, units = "in")
  return(p_weight2)
}


weighted.unifrac.dist = UniFrac(ps_transform ,  weighted = T, normalized = TRUE, parallel = FALSE,fast = TRUE)
Adonis_group_control_one_category <- adonis2(weighted.unifrac.dist~Time_num*Total_LIBRA_Score*BMI_Current*age*as.factor(sex)*as.factor(Ever_Homeless)*Participant, data = ps.meta_mindsource, permutations = 9999)
adonis = Adonis_group_control_one_category$`Pr(>F)`[1]
adonis


unweighted.unifrac.dist = UniFrac(ps_transform ,  weighted = F, normalized = TRUE, parallel = FALSE,fast = TRUE)
unAdonis_group_control_one_category <- adonis2(unweighted.unifrac.dist~Time_num*Total_LIBRA_Score*BMI_Current*age*as.factor(sex)*as.factor(Ever_Homeless)*Participant, data = ps.meta_mindsource, permutations = 9999)
unadonis = unAdonis_group_control_one_category$`Pr(>F)`[1]
unadonis

p_weight_group = weighted_fun(phylo_mindsource_beta_filter_tranform,"PCoA", T, "LIBRA_Half",paste0("Weighted UniFrac LIBRA Groups") ,paste0(path, "/weighted_mindsource_libra_half_",Sys.Date(),".jpeg"))
p_unweight_group = weighted_fun(phylo_mindsource_beta_filter_tranform,"PCoA", F, "LIBRA_Half",paste0("Unweighted UniFrac LIBRA Groups") ,paste0(path, "/unweighted_mindsource_libra_half_",Sys.Date(),".jpeg"))

p_weight_part = weighted_fun_part(phylo_mindsource_beta_filter_tranform,"PCoA", T, "Participant",paste0("Weighted UniFrac Participant") ,paste0(path, "/weighted_mindsource_Participant_",Sys.Date(),".jpeg"))
p_unweight_part = weighted_fun_part(phylo_mindsource_beta_filter_tranform,"PCoA", F, "Participant",paste0("Unweighted UniFrac Participant") ,paste0(path, "/unweighted_mindsource_Participant",Sys.Date(),".jpeg"))



#Functions for beta, heatmap, and stack bar plot
p_weight_group = weighted_fun(phylo_mindsource_beta_filter_tranform,"PCoA", T, "LIBRA_Half",paste0("Weighted UniFrac LIBRA Groups") ,paste0(path, "/weighted_mindsource_libra_half_",Sys.Date(),".jpeg"))
p_unweight_group = weighted_fun(phylo_mindsource_beta_filter_tranform,"PCoA", F, "LIBRA_Half",paste0("Unweighted UniFrac LIBRA Groups") ,paste0(path, "/unweighted_mindsource_libra_half_",Sys.Date(),".jpeg"))
p_aitchison_group = aitchison_fun(phylo_mindsource_beta_filter_tranform,"PCoA", F, "LIBRA_Half",paste0("Unweighted UniFrac LIBRA Groups") ,paste0(path, "/unweighted_mindsource_libra_half_",Sys.Date(),".jpeg"))



p_weight_part = weighted_fun_part(phylo_mindsource_beta_filter_tranform,"PCoA", T, "Participant",paste0("Weighted UniFrac Participant") ,paste0(path, "/weighted_mindsource_Participant_",Sys.Date(),".jpeg"))
p_unweight_part = weighted_fun_part(phylo_mindsource_beta_filter_tranform,"PCoA", F, "Participant",paste0("Unweighted UniFrac Participant") ,paste0(path, "/unweighted_mindsource_Participant",Sys.Date(),".jpeg"))

# NEW: Aitchison ordination plot, colored by Participant (uses raw counts)
p_aitch_part = aitchison_fun_part(physeq_object_mindsource_final,
                                  "PCoA",
                                  "Participant",
                                  paste0("Aitchison Distance Participant") ,
                                  paste0(path, "/aitchison_mindsource_Participant_",Sys.Date(),".jpeg"))

coord_beta_w_part = p_weight_part + labs(tag="B") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
coord_beta_uw_part = p_unweight_part + labs(tag="A") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())

# NEW: tagged Aitchison panel, following the same pattern as the UniFrac panels above
coord_beta_aitch_part = p_aitch_part$plot + labs(tag="C") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())

ggsave(paste0(path, "supp1B_beta_w_part_figure_", Sys.Date(), ".jpeg"), coord_beta_w_part, height = 6, width = 8, units="in",device="jpeg")
ggsave(paste0(path, "supp1A_beta_uw_part_figure_", Sys.Date(), ".jpeg"), coord_beta_uw_part, height = 6, width = 8, units="in",device="jpeg")
ggsave(paste0(path, "supp1C_beta_aitch_part_figure_", Sys.Date(), ".jpeg"), coord_beta_aitch_part, height = 6, width = 8, units="in",device="jpeg")




#######
#####update from here
#############
library(phyloseq)
library(vegan)
library(microbiome)
library(ggplot2)
library(permaBrioche)  # devtools::install_github("biobakery/permaBrioche")

weighted_fun = function (data_set,da_method, weight_T_F, color_group, title_name,path_out){
  weight_ordination = ordinate(data_set, da_method, "unifrac", weighted=weight_T_F)
  p_weight <- plot_ordination(data_set,weight_ordination, type="samples", color=color_group, title=title_name)+ theme_bw()
  p_weight2 = p_weight + 
    stat_ellipse(type = "norm", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 16, face="bold"))+
    labs(color="Brain Health Group \n(LIBRA Score)") +
    #theme(plot.tag.position = c(.885, .98), plot.tag = element_text(size=7.25))+
    scale_colour_manual(
      values = c( "#A77903", 
                  "#1B3558"),
      labels= c("Lower Risk Group", "Higher Risk Group")
    )
  ggsave(path_out, p_weight2,device="jpeg", height = 6, width = 8, units = "in")
  return(p_weight2)
}

?ordinate
weighted_fun_part = function (data_set,da_method, weight_T_F, color_group, title_name,path_out){
  weight_ordination = ordinate(data_set, da_method, "unifrac", weighted=weight_T_F)
  p_weight <- plot_ordination(data_set,weight_ordination, type="samples", color=color_group, title=title_name)+ theme_bw()
  p_weight2 = p_weight + 
    stat_ellipse(type = "norm", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 12, face="bold"))+
    labs(color="Participant")
  #theme(plot.tag.position = c(.885, .98), plot.tag = element_text(size=7.25))+
  ggsave(path_out, p_weight2,device="jpeg", height = 6, width = 8, units = "in")
  return(p_weight2)
}

# ---------------------------------------------------------------------------
# NEW: Aitchison-distance versions of the ordination helper functions.
# Aitchison distance = Euclidean distance on CLR-transformed counts. It needs
# raw (or lightly filtered) counts, not already-normalized relative
# abundance, and a pseudocount to avoid log(0). A precomputed vegan `dist`
# object can be handed straight to phyloseq::ordinate() via `distance = `.
# ---------------------------------------------------------------------------
aitchison_fun = function (data_set, da_method, color_group, title_name, path_out){
  otu <- as(otu_table(data_set), "matrix")
  if (taxa_are_rows(data_set)) otu <- t(otu)          # vegan wants samples in rows
  aitch_dist <- vegdist(otu + 1, method = "aitchison") # pseudocount avoids log(0)
  
  aitch_ordination = ordinate(data_set, da_method, distance = aitch_dist)
  p_aitch <- plot_ordination(data_set, aitch_ordination, type="samples", color=color_group, title=title_name) + theme_bw()
  p_aitch2 = p_aitch + 
    stat_ellipse(type = "norm", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 16, face="bold"))+
    labs(color="Brain Health Group \n(LIBRA Score)") +
    scale_colour_manual(
      values = c( "#A77903", 
                  "#1B3558"),
      labels= c("Lower Risk Group", "Higher Risk Group")
    )
  ggsave(path_out, p_aitch2,device="jpeg", height = 6, width = 8, units = "in")
  return(list(plot = p_aitch2, dist = aitch_dist))
  
}

aitchison_fun_part = function (data_set, da_method, color_group, title_name, path_out){
  otu <- as(otu_table(data_set), "matrix")
  if (taxa_are_rows(data_set)) otu <- t(otu)
  aitch_dist <- vegdist(otu + 1, method = "robust.aitchison")
  
  aitch_ordination = ordinate(data_set, da_method, distance = aitch_dist)
  p_aitch <- plot_ordination(data_set, aitch_ordination, type="samples", color=color_group, title=title_name) + theme_bw()
  p_aitch2 = p_aitch + 
    stat_ellipse(type = "norm", linetype = 2) +
    theme_bw() +
    theme(text = element_text(size = 12, face="bold"))+
    labs(color="Participant")
  ggsave(path_out, p_aitch2,device="jpeg", height = 6, width = 8, units = "in")
  return(list(plot = p_aitch2, dist = aitch_dist))
}


ps.meta_mindsource = meta(physeq_object_mindsource_final)
#changes data to abundance counts and moves 0 to a low value
phylo_mindsource_beta_filter_tranform = transform_sample_counts(physeq_object_mindsource_final, function(x) x/sum(x))

ps_transform = transform_sample_counts(physeq_object_mindsource_final, function(x) 1E6 * x/sum(x))
ps.meta_mindsource = meta(physeq_object_mindsource_final)

# permaBrioche's PERMANOVA_repeat_measures needs a column that uniquely
# identifies each sample/row in the metadata (not the repeated-measures
# blocking variable). Adjust "SampleID" below to match your actual column
# name if it's already present in ps.meta_mindsource.
ps.meta_mindsource$SampleID <- rownames(ps.meta_mindsource)

weighted.unifrac.dist = UniFrac(ps_transform ,  weighted = T, normalized = TRUE, parallel = FALSE,fast = TRUE)
unweighted.unifrac.dist = UniFrac(ps_transform ,  weighted = F, normalized = TRUE, parallel = FALSE,fast = TRUE)

# ---------------------------------------------------------------------------
# Robust Aitchison distance (vegan). Unlike the plain "aitchison" method used
# for the ordination plots above, "robust.aitchison" handles zeros directly
# (CLR computed only over the non-zero parts of each sample), so no
# pseudocount is added here.
# ---------------------------------------------------------------------------
otu_all <- as(otu_table(physeq_object_mindsource_final), "matrix")
if (taxa_are_rows(physeq_object_mindsource_final)) otu_all <- t(otu_all)
robust.aitchison.dist <- vegdist(otu_all, method = "robust.aitchison")

# ---------------------------------------------------------------------------
# Repeated-measures PERMANOVA via permaBrioche, replacing adonis2().
# Participant is now handled by blocking_variable (permutations are
# restricted within/between subjects as appropriate) rather than being
# included as a term in the formula, since it defined the old error strata.
# ---------------------------------------------------------------------------
res_perm_weighted <- PERMANOVA_repeat_measures(
  formula           = weighted.unifrac.dist ~ Time_num*Total_LIBRA_Score*BMI_Current*age*as.factor(sex)*as.factor(Ever_Homeless),
  data              = ps.meta_mindsource,
  sample_id         = "Sample_ID",
  blocking_variable = "Participant",
  permutations      = 9999
)
res_perm_weighted

res_perm_unweighted <- PERMANOVA_repeat_measures(
  formula           = unweighted.unifrac.dist ~ Time_num*Total_LIBRA_Score*BMI_Current*age*as.factor(sex)*as.factor(Ever_Homeless),
  data              = ps.meta_mindsource,
  sample_id         = "SampleID",
  blocking_variable = "Participant",
  permutations      = 9999
)
res_perm_unweighted

res_perm_aitchison <- PERMANOVA_repeat_measures(
  formula           = robust.aitchison.dist ~ BMI_Current+Total_LIBRA_Score,#*BMI_Current*age*as.factor(sex)*as.factor(Ever_Homeless),
  data              = ps.meta_mindsource,
  sample_id         = "Sample_ID",
  blocking_variable = "Participant",
  permutations      = 9999
)
res_perm_aitchison



#Functions for beta, heatmap, and stack bar plot
p_weight_group = weighted_fun(phylo_mindsource_beta_filter_tranform,"PCoA", T, "LIBRA_Half",paste0("Weighted UniFrac LIBRA Groups") ,paste0(path, "/weighted_mindsource_libra_half_",Sys.Date(),".jpeg"))
p_unweight_group = weighted_fun(phylo_mindsource_beta_filter_tranform,"PCoA", F, "LIBRA_Half",paste0("Unweighted UniFrac LIBRA Groups") ,paste0(path, "/unweighted_mindsource_libra_half_",Sys.Date(),".jpeg"))

# NEW: Aitchison ordination plot, colored by LIBRA group (uses raw counts)
p_aitch_group = aitchison_fun(physeq_object_mindsource_final,"PCoA","LIBRA_Half",paste0("Aitchison Distance LIBRA Groups") ,paste0(path, "/aitchison_mindsource_libra_half_",Sys.Date(),".jpeg"))
p_aitch_group_plot = p_aitch_group$plot


p_weight_part = weighted_fun_part(phylo_mindsource_beta_filter_tranform,"PCoA", T, "Participant",paste0("Weighted UniFrac Participant") ,paste0(path, "/weighted_mindsource_Participant_",Sys.Date(),".jpeg"))
p_unweight_part = weighted_fun_part(phylo_mindsource_beta_filter_tranform,"PCoA", F, "Participant",paste0("Unweighted UniFrac Participant") ,paste0(path, "/unweighted_mindsource_Participant",Sys.Date(),".jpeg"))

# NEW: Aitchison ordination plot, colored by Participant (uses raw counts)
p_aitch_part = aitchison_fun_part(physeq_object_mindsource_final,"PCoA","Participant",paste0("Aitchison Distance Participant") ,paste0(path, "/aitchison_mindsource_Participant_",Sys.Date(),".jpeg"))

coord_beta_w_part = p_weight_part + labs(tag="B") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())
coord_beta_uw_part = p_unweight_part + labs(tag="A") +  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())

# NEW: tagged Aitchison panel, following the same pattern as the UniFrac panels above
coord_beta_aitch_part = p_aitch_part$plot + #labs(tag="C") +  
  theme(plot.tag = element_text(size = 14, face = "bold"),plot.title = element_blank())

ggsave(paste0(path, "supp1B_beta_w_part_figure_", Sys.Date(), ".jpeg"), coord_beta_w_part, height = 6, width = 8, units="in",device="jpeg")
ggsave(paste0(path, "supp1A_beta_uw_part_figure_", Sys.Date(), ".jpeg"), coord_beta_uw_part, height = 6, width = 8, units="in",device="jpeg")
ggsave(paste0(path, "supp1C_beta_aitch_part_figure_", Sys.Date(), ".jpeg"), coord_beta_aitch_part, height = 6, width = 8, units="in",device="jpeg")
ggsave(paste0(path, "supp1d_beta_rpca_part_figure_", Sys.Date(), ".tiff"), phylo_rpca_part, height = 6, width = 8, units="in",device="tiff")


##########Aitchison Distance Plot
aitchison_data_long = data.frame(layer_data(beta_long$Aitchison, 1)) %>% select(group, colour, x,y) %>% mutate(LIBRA_Half = case_when(colour=="#E31A1C" ~ "Lower",
                                                                                                                                      T~"Upper"))

aitchison_distance = ggplot(aitchison_data_long, aes(x=x, y=y,group=LIBRA_Half))+
  geom_line(aes(group=group),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = aitchison_data_long)+
  labs(tag = "A", x="Time (Months)", y="Aitchison Distance\n from Baseline", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c(  "#A77903", 
                 "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )


#'#0047AB' '#006666' '#002D72' '#FFB81C' "#1F3A5F"


lookup <- tibble(
  group = 1:17,
  Participant = c("P11","P14","P17","P19","P24","P25","P27","P28","P30",
                  "P31","P34","P42","P45","P48","P51","P52","P61")
)

aitchison_data_long2 <- aitchison_data_long %>% left_join(lookup, by = "group")
meta_and_alpha_libra_beta_long_filter = meta_and_alpha_libra_beta_long %>%
  select(Participant, Total_LIBRA_Score, age, sex, BMI_Current, Ever_Homeless) %>% distinct()
aitchison_merge = left_join(aitchison_data_long2, meta_and_alpha_libra_beta_long_filter, by="Participant")


#beta lmer
library(lme4)
library(lmerTest)

aitchison_merge$Ever_Homeless = factor(aitchison_merge$Ever_Homeless, levels=c("No", "Unknown", "Yes"))

aitchison_model <- lmer(y ~ Total_LIBRA_Score*x+ age + as.factor(sex) + BMI_Current + factor(Ever_Homeless) + (1 | Participant), data = aitchison_merge)
summary(aitchison_model)
(anova(aitchison_model))

aitchison_model_gamm <- gamm4(
  y ~ t2(x, Total_LIBRA_Score) + age + as.factor(sex) + 
    BMI_Current + as.factor(Ever_Homeless),
  random = ~ (1 | Participant),
  data = aitchison_merge
)

summary(aitchison_model_gamm$gam)   # smooth-term significance, deviance explained, etc.
summary(aitchison_model_gamm$mer)   # the underlying lmer fit, if you want the random-effect variance
plot(aitchison_model_gamm$gam, pages = 1)  # visualize the fitted smooth surface



sup_table_uni_dist_model = flextable(tidy(uni_model))
sup_table_uni_dist_model = set_caption(sup_table_uni_dist_model, caption = "Supplementary Table x. Mixed Effects Model Unweighted Unifrac Distance") %>% autofit()
save_as_docx(sup_table_uni_dist_model, path = paste0(path, "Supplementary_Table_uni_mod._", Sys.Date(), ".docx"))



sup_table_uni_dist_model = flextable(tidy(uni_model))
sup_table_uni_dist_model = set_caption(sup_table_uni_dist_model, caption = "Supplementary Table x. Mixed Effects Model Unweighted Unifrac Distance") %>% autofit()
save_as_docx(sup_table_uni_dist_model, path = paste0(path, "Supplementary_Table_uni_mod._", Sys.Date(), ".docx"))




uni_data_long = data.frame(layer_data(beta_long$UniFrac, 1)) %>% select(group, colour, x,y) %>% mutate(LIBRA_Half = case_when(colour=="#E31A1C" ~ "Lower",
                                                                                                                              T~"Upper"))

uni_data_long2 <- uni_data_long %>% left_join(lookup, by = "group")
meta_and_alpha_libra_beta_long_filter = meta_and_alpha_libra_beta_long %>%
  select(Participant, Total_LIBRA_Score, age, sex, BMI_Current, Ever_Homeless) %>% distinct()
uni_merge = left_join(uni_data_long2, meta_and_alpha_libra_beta_long_filter, by="Participant")


uni_distance = ggplot(uni_data_long, aes(x=x, y=y,group=LIBRA_Half))+
  geom_line(aes(group=group),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = uni_data_long)+
  labs(tag = "A", x="Time (Months)", y="Unweighted Unifrac\n Distance from Baseline", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c(  "#A77903", 
                 "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )


#'#0047AB' '#006666' '#002D72' '#FFB81C' "#1F3A5F"


#beta lmer
library(lme4)
library(lmerTest)

aitchison_model <- lmer(y ~ Total_LIBRA_Score*x+ age + as.factor(sex) + BMI_Current + factor(Ever_Homeless) + (1 | Participant), data = uni_merge)
summary(aitchison_model)
(anova(aitchison_model))

aitchison_model_gamm <- gamm4(
  y ~ t2(x, Total_LIBRA_Score) + age + as.factor(sex) + 
    BMI_Current + as.factor(Ever_Homeless),
  random = ~ (1 | Participant),
  data = uni_merge
)

summary(aitchison_model_gamm$gam)   # smooth-term significance, deviance explained, etc.
summary(aitchison_model_gamm$mer)   # the underlying lmer fit, if you want the random-effect variance
plot(aitchison_model_gamm$gam, pages = 1)  # visualize the fitted smooth surface



uni_model <- lmer(y ~ LIBRA_Half*x+ (1 | group), data = uni_data_long)
summary(uni_model)
(anova(uni_model))

sup_table_uni_dist_model = flextable(tidy(uni_model))
sup_table_uni_dist_model = set_caption(sup_table_uni_dist_model, caption = "Supplementary Table x. Mixed Effects Model Unweighted Unifrac Distance") %>% autofit()
save_as_docx(sup_table_uni_dist_model, path = paste0(path, "Supplementary_Table_uni_mod._", Sys.Date(), ".docx"))





wuni_data_long = data.frame(layer_data(beta_long$WUniFrac, 1)) %>% select(group, colour, x,y) %>% mutate(LIBRA_Half = case_when(colour=="#E31A1C" ~ "Lower",
                                                                                                                                T~"Upper"))

wuni_data_long2 <- wuni_data_long %>% left_join(lookup, by = "group")
meta_and_alpha_libra_beta_long_filter = meta_and_alpha_libra_beta_long %>%
  select(Participant, Total_LIBRA_Score, age, sex, BMI_Current, Ever_Homeless) %>% distinct()
wuni_merge = left_join(wuni_data_long2, meta_and_alpha_libra_beta_long_filter, by="Participant")

wuni_distance = ggplot(wuni_data_long, aes(x=x, y=y,group=LIBRA_Half))+
  geom_line(aes(group=group),color="lightgrey")+
  #geom_point()+
  #facet_wrap(~Participant)+
  geom_smooth(aes(color=LIBRA_Half),se=F, data = wuni_data_long)+
  labs(tag = "A", x="Time (Months)", y="Weighted Unifrac\n Distance from Baseline", color="Brain Health Group \n(LIBRA Score)") +
  theme(text = element_text(size = 16, face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  scale_x_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12)) +
  scale_colour_manual(
    values = c(  "#A77903", 
                 "#1B3558"),
    labels= c("Lower Risk Group", "Higher Risk Group")
  )



aitchison_model <- lmer(y ~ Total_LIBRA_Score*x+ age + as.factor(sex) + BMI_Current + factor(Ever_Homeless) + (1 | Participant), data = wuni_merge)
summary(aitchison_model)
(anova(aitchison_model))

aitchison_model_gamm <- gamm4(
  y ~ t2(x, Total_LIBRA_Score) + age + as.factor(sex) + 
    BMI_Current + as.factor(Ever_Homeless),
  random = ~ (1 | Participant),
  data = wuni_merge
)

summary(aitchison_model_gamm$gam)   # smooth-term significance, deviance explained, etc.
summary(aitchison_model_gamm$mer)   # the underlying lmer fit, if you want the random-effect variance
plot(aitchison_model_gamm$gam, pages = 1)  # visualize the fitted smooth surface

#beta lmer
library(lme4)
library(lmerTest)


wuni_model <- lmer(y ~ LIBRA_Half*x+ (1 | group), data = wuni_data_long)
summary(wuni_model)
(anova(wuni_model))

sup_table_wuni_dist_model = flextable(tidy(wuni_model))
sup_table_wuni_dist_model = set_caption(sup_table_wuni_dist_model, caption = "Supplementary Table x. Mixed Effects Model Weighted Unifrac Distance") %>% autofit()
save_as_docx(sup_table_wuni_dist_model, path = paste0(path, "Supplementary_Table_wuni_mod._", Sys.Date(), ".docx"))

generate_beta_volatility_test_long(
  data.obj = data.obj,
  dist.obj = NULL,
  subject.var = "Participant",
  time.var = "Time_num", 
  group.var = "LIBRA_Half",
  dist.name = c("UniFrac","WUniFrac","Robust.Aitchison"),
)



# --- Aitchison distance (not native, add it manually) ---
otu_mat <- data.obj$feature.tab          # taxa in rows, samples in columns
aitch <- vegdist(t(otu_mat), method = "robust.aitchison")   # pseudocount for zeros
# or: vegdist(t(otu_mat), method = "robust.aitchison")   # zero-robust, no pseudocount
dist.obj = mStat_calculate_beta_diversity(data.obj,
                                          dist.name = c("BC", "Jaccard", "UniFrac", "GUniFrac", "WUniFrac", "JS")
                                          
)
dist.obj$Aitchison <- as.matrix(aitch)

# --- Now run the volatility test across all three ---
generate_beta_volatility_test_long(
  data.obj = data.obj,
  dist.obj = dist.obj,
  subject.var = "Participant",
  time.var = "Time_num",
  group.var = "LIBRA_Half",
  dist.name = c("UniFrac", "WUniFrac", "Aitchison")
)
