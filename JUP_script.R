######################
###EOGHAN JUP STUFF###
######################
### Load necessary packages
library("TCGAbiolinks")
library("SummarizedExperiment")
library("ggplot2")
library("biomaRt")
library("pheatmap")
library("RCurl")

# USE TCGAbiolinks package to download data from TCGA
# first identify what data to download
TCGAbiolinks:::getGDCprojects()$project_id
TCGAbiolinks:::getProjectSummary("TCGA-BRCA") # BREAST CANCER DATA ONLY

query <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "HTSeq - FPKM-UQ", legacy = FALSE) # last version of gene exp data (FPKM-UQ)

# Download data, will take few minutes. It will be saved to your working directory
setwd("C:/Users/.../Desktop/Eoghan") # change working directory to whatever u want (better create a new one)
GDCdownload(query)

# prepare data and save it (so next time you dont need to donwload again, just load it). Will take a while (even more than download)
# more RAM, faster
data <- GDCprepare(query, save = TRUE, save.filename = "probando_Eoghan.rda", summarizedExperiment = TRUE)

#################################################################################
# after first time you download, you can skip the previous steps and start HERE #
#################################################################################
# change path to wherever you saved it
load("C:/Users/.../Desktop/Eoghan/probando_Eoghan.rda") 

# extract needed data from SummarizedExperiment object
# in the downloaded data, theres "Primary solid Tumor", "Solid Tissue Normal" and "Metastatic" labelled samples
tcga_brca_all <- as.data.frame(assay(data))  # gene expression data
brca_subtype <- as.data.frame(data$subtype_BRCA_Subtype_PAM50) #BC subtype, taken from this paper: doi.org/10.1016/j.ccell.2018.03.014
tissue_type <- as.data.frame(data$definition) # tissue type of each sample

View(tcga_brca_all) # genes in rows, samples in columns
View(brca_subtype) # same order
View(tissue_type) # same order

# load IRE1_gene_sig 
breast_IRE1_sig <- read.csv(text = getURL("https://raw.githubusercontent.com/xaitorx/Analysis_JUP_TCGA_CCLE/master/data/breast_IRE1_sig.csv"))

# calculate IRE1 score for each sample
# mean of IRE1 positively correlated genes - mean of negatively correlated
tcga_brca_sig <- tcga_brca_all[as.character(breast_IRE1_sig$ensembl_gene_id),]
IRE1_score_pos <- log2(colMeans(na.omit(tcga_brca_sig[which(breast_IRE1_sig$sig == "TRUE"),])))
IRE1_score_neg <- log2(colMeans(na.omit(tcga_brca_sig[which(breast_IRE1_sig$sig == "FALSE"),])))
IRE1_score <- IRE1_score_pos - IRE1_score_neg

### Retrieve JUP expression data (they are annotated with ensembl ID)
JUP <- as.data.frame(t(log2(tcga_brca_all[c("ENSG00000173801"),]))) #if you want to check other gene, just change ENSEMBL ID

# create data frame including all info together
JUP_all <- data.frame(JUP = JUP$ENSG00000173801,
                      PAM_50 = brca_subtype$`data$subtype_BRCA_Subtype_PAM50`,
                      tissue_type = tissue_type,
                      IRE1_score = IRE1_score,
                      IRE1_score_scaled = scale(IRE1_score), # scaled around mean, its the same just look better for plotting. values are relative to mean of these set of samples
                      row.names = row.names(JUP))

View(JUP_all)

write.csv(JUP_all, "JUP_all.csv") # save this table if you want, is the data that is plotted in the following graphs
# some of these steps could be done together, but I think step by step is easier to see/understand

# now lets plot things
#############################################
##### JUP mRNA expression across subtypes ###
#############################################
# data only for tumor samples
JUP_all_tum <- JUP_all[which(tissue_type$`data$definition` == "Primary solid Tumor"),] # create index with position of tumor samples only

# boxplots with ggplot2
# you can customize the plot as you want, ggplot2 function accepts many arguments
p <- ggplot(JUP_all_tum, aes(JUP_all_tum$PAM_50, JUP_all_tum$JUP, fill = JUP_all_tum$PAM_50))
p + geom_boxplot(outlier.shape = NA, size =2, alpha = 0.5) + 
  geom_point(aes(), size = 2, position = position_jitterdodge()) + 
  labs(x="", y = "JUP mRNA expression (log2(FPKM-UQ))", element_text(face = "bold", angle = 0)) + 
  scale_fill_manual(values =  c("darkred", "grey40", "cyan", "darkblue", "khaki", "grey40")) + # you can change colors here if you want
  theme(panel.background = element_rect( colour = "grey50"), axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + 
  theme_bw() +
  theme(text = element_text(size=20))

# we can also compare JUP expression across "Primary solid Tumor", "Solid Tissue Normal" and "Metastatic" labelled samples
p <- ggplot(JUP_all, aes(JUP_all$data.definition, JUP_all$JUP, fill = JUP_all$data.definition))
p + geom_boxplot(outlier.shape = NA, size =2, alpha = 0.5) + 
  geom_point(aes(), size = 2, position = position_jitterdodge()) + 
  labs(x="", y = "JUP mRNA expression", element_text(face = "bold", angle = 0)) + 
  scale_fill_manual(values =  c("green", "red", "black")) + # you can change colors here if you want
  theme(panel.background = element_rect( colour = "grey50"), axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + 
  theme_bw() +
  theme(text = element_text(size=20))

### Now lets do IRE1 score by subtype for example
p <- ggplot(JUP_all_tum, aes(JUP_all_tum$PAM_50, JUP_all_tum$IRE1_score_scaled, fill = JUP_all_tum$PAM_50))
p + geom_boxplot(outlier.shape = NA, size =2, alpha = 0.5) + 
  geom_point(aes(), size = 2, position = position_jitterdodge()) + 
  labs(x="", y = "Predicted IRE1 activity", element_text(face = "bold", angle = 0)) + 
  scale_fill_manual(values =  c("darkred", "grey40", "cyan", "darkblue", "khaki", "grey40")) + # you can change colors here if you want
  theme(panel.background = element_rect( colour = "grey50"), axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + 
  theme_bw() +
  theme(text = element_text(size=20))


#####################################
### COMPARE JUP EXP VS IRE1 SCORE ###
#####################################
# only for tumor samples
p <- ggplot(JUP_all_tum, aes(JUP_all_tum$IRE1_score_scaled, JUP_all_tum$JUP))
p + geom_point(aes(colour = JUP_all_tum$PAM_50, size = 2, alpha = 0.5)) +
  labs(x = "IRE1 activity z score", y = "JUP mRNA expression (log2(FPKM-UQ))", element_text(face = "bold", angle = 0)) +
  geom_smooth(method="lm", color = "red") + theme(axis.line = element_line(size = 3, colour = "grey80")) +
  scale_colour_manual(values = c("darkred", "black", "cyan", "darkblue", "khaki",  "grey50")) +  
  guides(colour = guide_legend("BC subtypes", override.aes = list(size=5), title.theme = element_text(face = "bold", angle = 0)), size = "none") 

# calculate Pearsons R^2 and P-value
cor(JUP_all_tum$IRE1_score_scaled, JUP_all_tum$JUP)
cor.test(JUP_all_tum$IRE1_score_scaled, JUP_all_tum$JUP)$p.value

### Add these results to the plot
p + geom_point(aes(colour = JUP_all_tum$PAM_50, size = 2, alpha = 0.5)) +
  labs(x = "IRE1 activity z score", y = "JUP mRNA expression (log2(FPKM-UQ))", element_text(face = "bold", angle = 0)) +
  geom_smooth(method="lm", color = "red") + theme(axis.line = element_line(size = 3, colour = "grey80")) +
  scale_colour_manual(values = c("darkred", "black", "cyan", "darkblue", "khaki",  "grey50")) +  
  guides(colour = guide_legend("BC subtypes", override.aes = list(size=5), title.theme = element_text(face = "bold", angle = 0)), size = "none") +
  annotate('text', x = -2, y = 17.5 ,label = "R^{2}==-0.44 ",parse = TRUE, size=6) + # need to manually adjust where to put the annotation 
  annotate('text', x = -2, y = 17,label = "P-value==8.122e-56 ",parse = TRUE, size=6) + # same
  theme(text = element_text(size=20))
  

##############################################
### JUP levels in breast cancer cell lines ###
##############################################
# data from Cancer Cel Line Enciclopedia (CCLE), Broad Institute
# download from here: https://portals.broadinstitute.org/ccle/data
# 2 files: Cell_lines_annotations_20181226.txt ; CCLE_RNAseq_genes_rpkm_20180929.gct.gz
# save them to working directory
# Cell_lines_annotations_20181226 <- read.delim("C:/Users/.../Desktop/Eoghan/Cell_lines_annotations_20181226.txt")
# df <- read.table("C:/Users/.../Desktop/Eoghan/CCLE_RNAseq_genes_rpkm_20180929.gct", sep = '\t',header = TRUE, skip = 2)
# Cell_lines_breast <- subset(Cell_lines_annotations_20181226, Cell_lines_annotations_20181226$Disease == "breast_cancer")
# breast_cell_lines <- cbind(df[,1:2] ,df[, as.character(Cell_lines_breast$CCLE_ID)])
# molecular subtypes taken from  doi: 10.7150/jca.18457

breast_cell_lines <- read.csv(text = getURL("https://raw.githubusercontent.com/xaitorx/Analysis_JUP_TCGA_CCLE/master/data/breast_cell_lines.csv"))
breast_cell_lines_subtype <- read.csv(text = getURL("https://raw.githubusercontent.com/xaitorx/Analysis_JUP_TCGA_CCLE/master/data/breast_cell_lines_subtype.csv"))

# calculate IRE1 score for each cell line
breast_IRE1_sig$SYMBOL <- as.character(breast_IRE1_sig$SYMBOL)
breast_cell_lines$Description <- as.character(breast_cell_lines$Description)
breast_cell_lines_sig <- merge(breast_IRE1_sig, breast_cell_lines, by.x = 1, by.y = 2, sort = FALSE)

# stupid CCLE annotation, duplicated entries. also 5 genes missing (?)
ind <- duplicated(breast_cell_lines_sig$SYMBOL)
breast_cell_lines_sig <- breast_cell_lines_sig[!ind,]

# keep only numeric values to calculate means
breast_cell_lines_uuu <- breast_cell_lines_sig[,5:ncol(breast_cell_lines_sig)]
row.names(breast_cell_lines_uuu) <- breast_cell_lines_sig$SYMBOL

# calculate IRE1 activity for each cell line
IRE1_score_pos_cl <- log2(colMeans(na.omit(breast_cell_lines_uuu[which(breast_cell_lines_sig$sig == "TRUE"),])))
IRE1_score_neg_cl <- log2(colMeans(na.omit(breast_cell_lines_uuu[which(breast_cell_lines_sig$sig == "FALSE"),])))
IRE1_score_cl <- as.data.frame(IRE1_score_pos_cl - IRE1_score_neg_cl)

# get expression values for JUP
JUP_cell_lines <- as.data.frame(t(subset(breast_cell_lines, breast_cell_lines$Description == "JUP")[,3:38])) # here is annotated with gene symbol

# summarize all in one table
JUP_all_cell_lines <- data.frame(JUP = JUP_cell_lines$`44971`,
                                 IRE1_score = IRE1_score_cl$`IRE1_score_pos_cl - IRE1_score_neg_cl`,
                                 subtype = breast_cell_lines_subtype$subtype,
                                 row.names = row.names(JUP_cell_lines))

JUP_all_cell_lines <- JUP_all_cell_lines[order(JUP_all_cell_lines$IRE1_score),]

write.csv(JUP_all_cell_lines, "JUP_all_cell_lines.csv") # save this table if you want, is the data that is plotted in the next graph

####################
### DRAW HEATMAP ###
####################
#although if we are playing with only 2 variables probably sccaterplot is easier interpretation
# everything can be customized: colors, borders, dendrograms, labels, etc. 
# prepare the annotation for the columns
annotation_columnas <- as.data.frame(JUP_all_cell_lines[,3])
colnames(annotation_columnas) <- "subtype"
annotation_columnas$subtype <- as.character(annotation_columnas$subtype)
annotation_columnas$subtype[c(3,19,28)] <- "Nan"
row.names(annotation_columnas) <- row.names(JUP_all_cell_lines)
ann_colors = list( subtype = c( Basal = "darkred", Her2 = "black", Nan ="grey50", LumA = "cyan", LumB = "darkblue" )) # you can change colores here, maybe keep consistent with the other plots

# No dendrogram (arranged columns ugly), Everything orderes acording to IRE1 score, from higher to lower
my_palette <- colorRampPalette(c("navy","white","red"))(n = 299) # scale blue to red, or whatever colors
pheatmap(t(JUP_all_cell_lines[,1:2]), 
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         scale = "row", 
         color = my_palette, 
         border_color = "black",
         fontsize_col = 12, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col =  annotation_columnas,
         annotation_colors = ann_colors,
         cellwidth = 18,
         cellheight = 18,
         angle_col = 45)

