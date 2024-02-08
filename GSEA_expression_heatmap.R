# Chenck ggvenn at https://bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_Nov22/Bulk_RNAseq_Course_Base/Markdowns/11_Annotation_and_Visualisation.html
#for your tfg!!!

# Check metascape!!


## Subset by Pathways
library(ComplexHeatmap)
myGSEA.df_H[myGSEA.df_H$ID=="HALLMARK_G2M_CHECKPOINT", 11]
### extract values from that list, deleting the "/" symbol
geneset_list<-unlist(myGSEA.df_H[myGSEA.df_H$ID=="HALLMARK_G2M_CHECKPOINT", 11]) %>% strsplit("/")
geneset_list<-geneset_list$core_enrichment
g = intersect(geneset_list, rownames(diffGenes))
g_index<-rownames(diffGenes) %in% g
geneset_expr<-diffGenes[g,]
s = rowttests(diffGenes, factor(condition))[, "statistic"]
## Complex heatmap for pathways
Heatmap(geneset_expr, name = "z-score", top_annotation = HeatmapAnnotation(Condition = condition),
        column_title = "HALLMARK_G2M_CHECKPOINT",
        right_annotation = rowAnnotation("T Test" = anno_barplot(s[g_index], axis_param = list(direction = "reverse"), width = unit(2, "cm"))),
        column_split = condition)
