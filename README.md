# Analysis_JUP_TCGA_CCLE
 JUP is a candidate IRE1 RIDD substrate in Triple Negative Breast Cancer. 
 Explore levels of JUP mRNA across breast cancer subtypes, cell lines. 
 Explore correlation with IRE1 activity. 

# Levels of IRE1 activity vs levels of JUP mRNA transcript

Across 1100 patient tumor samples from TCGA-BRCA project
Tumor subtype indicated by color code.

![Slide1] (https://github.com/xaitorx/Analysis_JUP_TCGA_CCLE/blob/master/pics/scatterplot_IRE1_JUP.jpg)


# Schematic workflow
![Slide3](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide3.JPG)

# PCA: Exploring dataset
![Slide4](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide4.JPG)
Biological groups show decent separation. Repeats 4&5 only metabolomic data

# OPLS-DA: Identifying sources of variation between DMSO - MKC

![Slide5](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide5.JPG)

![Slide6](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide6.JPG)
First principal component t[1] vs first orthogonal component to[1]
Y-axis shows the within group variation, and X-axis shows the between group variation (DMSO vs MKC)

![Slide7](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide7.JPG)

# Rank all variables according to loading scores in component 1

![Slide8](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide8.JPG)

![Slide9](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide9.JPG)

![Slide10](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide10.JPG)

![Slide11](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide11.JPG)

# Pathway enrichment – top enriched terms - REACTOME
Including all gene, miRNA, chemical compounds annotated to REACTOME pathways. On the preranked variables.
![Slide12](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide12.JPG)

![Slide13](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide13.JPG)

![Slide14](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide14.JPG)
Plenty of terms enriched at the bottom of the ranking (DMSO phenotype), but few at the top (MKC phenotype)
Explanation? 
Why top of ranking is more heterogeneous?


# Summarization of top enriched terms - REACTOME
Pathway information is inherently redundant, as genes often participate in multiple pathways, and databases may organize pathways hierarchically by including general and specific pathways with many shared genes. Consequently, pathway enrichment analysis often highlights several versions of the same pathway. Collapsing redundant pathways into a single biological theme simplifies interpretation.

cutoff: 0.1 FDR
![Slide15](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide15.JPG)

*Post-analysis – Overlapping signatures*: 
Add more gene sets to an existing network. This is done by calculating the overlap between gene sets of the current network and other signature gene set files, spotting significantly overlapping signatures. Help give biological context.
![Slide16](https://github.com/xaitorx/OPLS-DA_Integration/blob/master/pics/Slide16.JPG)
