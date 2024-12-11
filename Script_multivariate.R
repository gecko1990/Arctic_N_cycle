### Multivariate analysis  

# Load libraries:
library(ggplot2)
library(dplyr)
library("RColorBrewer")
library(data.table)
library(reshape2)
library(tidyverse)
library(zoo)
library(vegan)
library(ggrepel)

##############################

### Import genetic data ####
TPM_data <- read.csv("./TPM_data_13.arctic_bins.orftable", header = TRUE, sep="\t")

#Import contig info and taxonomy info
contigs_bin <- read.csv("./contigs_reference", header=FALSE, sep="\t")
gtdb_bins <- read.csv("./genome_name_GTDB.txt", header=TRUE, sep="\t")

#merge bins, contigs and bin taxonomy
gtdb_bins_contigs <- merge(gtdb_bins, contigs_bin, by.x="user_genome", by.y="V1", all=TRUE)
#Merge table
master_tab <- merge(TPM_data, gtdb_bins_contigs, by.x="Contig.ID", by.y="V3")

#Import contextual data
site_data <- read.table("./Sample_info.txt", header=TRUE, sep="\t")

#### Obtain table based on MAGs

master_tab$tax_MAG <- paste(master_tab$classification,master_tab$Genome_name_new,sep="user___")

### Tab for RNA

master_tab_RNA <- master_tab[,c(27,4,7,12,15,17,19,21)]
master_tab_RNA_summed <- master_tab_RNA %>%
  group_by(tax_MAG) %>%
  summarise(across(everything(), sum)) %>%
  ungroup()
master_tab_RNA_relative <- master_tab_RNA_summed %>%
  mutate(across(-tax_MAG, ~ log2(.), .names = "percent_{col}"))
master_tab_RNA_relative <- master_tab_RNA_relative[,c(1,9:15)]
master_tab_RNA_relative_matrix <- as.matrix(master_tab_RNA_relative[,-1])
row.names(master_tab_RNA_relative_matrix) <- master_tab_RNA_relative$tax_MAG
master_tab_RNA_relative_matrix_t <- t(master_tab_RNA_relative_matrix)
row.names(master_tab_RNA_relative_matrix_t) <- gsub("percent_", "", row.names(master_tab_RNA_relative_matrix_t))

### RDA

#Import environmental data and standarize them
env_data <- read.csv("./metadata_envs_RNA.txt", header = TRUE, row.names = 1, sep="\t", dec = ",")
env_data <- env_data[,c(5:9)]

t.env_data <- decostand(env_data, method = "standardize") # standarization lineal normal

#Do the RDA
RDA_RNA_MAG <- rda(master_tab_RNA_relative_matrix_t ~ ., data = t.env_data)
RDA_RNA_MAG_eigenvalues <- RDA_RNA_MAG$CCA$eig
RDA_RNA_MAG_eigenvalues_PC <- (RDA_RNA_MAG_eigenvalues / sum(RDA_RNA_MAG_eigenvalues)) * 100

# Extract scores for samples
sample_scores_RDA_RNA_MAG <- scores(RDA_RNA_MAG, display = "sites") %>% as.data.frame()
sample_scores_RDA_RNA_MAG$Variable <- "sites"
species_scores_RDA_RNA_MAG <- scores(RDA_RNA_MAG, display = "species") %>% as.data.frame()
species_scores_RDA_RNA_MAG$Variable <- "Species"
env_scores_RDA_RNA_MAG <- scores(RDA_RNA_MAG, display = "bp") %>% as.data.frame()
env_scores_RDA_RNA_MAG$Variable <- "Envs"

data_RDA_RNA_MAG <- rbind(sample_scores_RDA_RNA_MAG,species_scores_RDA_RNA_MAG, env_scores_RDA_RNA_MAG)
data_RDA_RNA_MAG$site <- rownames(data_RDA_RNA_MAG)
data_RDA_RNA_MAG$site <- gsub(pattern = "TPM.", "", data_RDA_RNA_MAG$site)

#Combine scores and contextual
data_RDA_RNA_MAG <- merge(data_RDA_RNA_MAG, site_data, by.x=c("site"), by.y=c("Sample_simple"), all.x=TRUE)
data_RDA_RNA_MAG$Season <- as.factor(data_RDA_RNA_MAG$Season)
data_RDA_RNA_MAG$Season <- factor(data_RDA_RNA_MAG$Season, levels=c("Winter","Transition","Summer"))

#Decomposes taxonomy
data_RDA_RNA_MAG$Tax <- data_RDA_RNA_MAG$site
data_RDA_RNA_MAG <- separate(data_RDA_RNA_MAG, Tax, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";")

#Import color code
color_code_tax <- read.csv("./color_code_15_classes_bin.txt", sep="\t", header = TRUE)
color_code_tax <- color_code_tax[,c(3,6)]
data_RDA_RNA_MAG_final <- merge(data_RDA_RNA_MAG, color_code_tax, by="Class", all.x=TRUE )

figure_col_tax <- c("Nitrososphaeria"="#E01853","Poseidoniia"="#e8517e","Actinobacteriota"="#F99392",
                    "Actinobacteriota"="#F99392","Bacteroidota"="#fdb462","Chloroflexota"="#CD672F",
                    "Gemmatimonadetes"="#b3cd2f","Latescibacterota"="#89CB6C","Marinisomatia"="#40A635",
                    "Myxococcota_A"="#FCD84A","Nitrospinota"="#E7E099","Planctomycetota"="#74add1",
                    "Planctomycetota"="#74add1","Alphaproteobacteria"="#1F74CD","Gammaproteobacteria"="#5D478B",
                    "SAR324"="#F6BFCA","Verrucomicrobiae"="#8B461D")


#Plot
RDA_RNA_MAG_plot <- ggplot(arrange(data_RDA_RNA_MAG_final, Domain, Phylum, Class, Taxa), aes(x=RDA1, y=RDA2)) +
  #arrows are also environmental
  geom_segment(data=subset(data_RDA_RNA_MAG_final, Variable =="Envs"), aes(x=0, y=0, xend= RDA1, yend= RDA2), size= 1, color="#99d8c9", arrow = arrow(length =unit (0.5, "cm"))) +
  #arrows are sites
  geom_segment(data=subset(data_RDA_RNA_MAG_final, Variable =="sites"), aes(x=0, y=0, xend= RDA1, yend= RDA2, color=Season), size= 1, arrow = arrow(length =unit (0.5, "cm"))) +
  geom_text_repel(aes(label = Date), max.overlaps = Inf,   # Ensure no label is hidden
                  box.padding = 0.5,    # Padding around labels
                  point.padding = 0.3,  # Padding between points and labels
                  segment.color = 'grey50', size=6) +
  scale_color_manual("Season", values= figure_col_general) +
  #Species are points
  geom_point(data=subset(arrange(data_RDA_RNA_MAG_final, Domain, Phylum, Class, Taxa), Variable =="Species"), aes(fill=(factor(Taxa, levels=unique(Taxa)))),  size= 3, shape=21) +
  scale_fill_manual("Taxa", values= figure_col_tax) +
  geom_text_repel(data=subset(data_RDA_RNA_MAG_final, Variable =="Envs"), aes(label = site), max.overlaps = Inf,   # Ensure no label is hidden
                  box.padding = 0.5,    # Padding around labels
                  point.padding = 0.3,  # Padding between points and labels
                  segment.color = 'grey50', size=6) +
  #geom_text_repel(data=subset(PC_log_vegan_data, Score =="sites"), aes(label = metadata$Depth..cm.)) +
  #geom_point(shape = 1,size = 5,colour = "black") +
  theme(panel.border= element_rect(color = "black", fill =NA)) +
  theme(panel.background = element_rect (fill= 'white', color= 'white')) +
  theme(axis.text = element_text(size=15),axis.title = element_text(size=15, face="bold"), legend.text = element_text(size=15)) 

RDA_RNA_MAG_plot
ggsave("../results/RDA_rna_MAG_2.pdf", plot=RDA_RNA_MAG_plot, dpi=300, scale=4)

### Heatmap and dendogram

library(pheatmap)

heatmap_matrix <- master_tab_RNA_relative_matrix_t
row.names(heatmap_matrix) <- gsub("TPM.|_RNA", "", row.names(heatmap_matrix))
colnames(heatmap_matrix) <- gsub(".*user___", "", colnames(heatmap_matrix))

sample_info <- subset(site_data, Molecule =="RNA")
sample_info <- data.frame(Season=sample_info[, 6], row.names = sample_info[, 3])

genome_info <- data.frame(Taxa=data_RDA_RNA_MAG_final[, 17], row.names = data_RDA_RNA_MAG_final[, 2])
row.names(genome_info) <- gsub(".*user___", "", row.names(genome_info))

#Create the pheatmap
heatmap_MAG <- pheatmap(heatmap_matrix, #show_colnames = FALSE,
         annotation_col= genome_info, annotation_row = sample_info,
         annotation_colors=list(Taxa=c("Nitrososphaeria"="#E01853","Poseidoniia"="#e8517e","Actinobacteriota"="#F99392",
                                      "Bacteroidota"="#fdb462","Chloroflexota"="#CD672F","Gemmatimonadetes"="#b3cd2f",
                                      "Latescibacterota"="#89CB6C","Marinisomatia"="#40A635","Myxococcota_A"="#FCD84A",
                                      "Nitrospinota"="#E7E099","Planctomycetota"="#74add1","Alphaproteobacteria"="#1F74CD","Gammaproteobacteria"="#5D478B",
                                      "SAR324"="#F6BFCA","Verrucomicrobiae"="#8B461D"),
                                Season=c("Winter"="#9FDCF2","Transition"="#FFE5AA","Summer"="#F4826C"))
         ,scale="none", cellwidth=4.5 #, cellheight = 20
         ,color = colorRampPalette(rev(brewer.pal(n = 7, name ="Spectral")))(100),
         cluster_rows= F, fontsize_col=5
         )

heatmap_MAG
ggsave("../results/heatmap_MAG.pdf", plot=heatmap_MAG, dpi=300, scale=2.75)