###########################################
###########################################
# Code to reproduce the fomputational analysis of this paper:
# "GLUT1 inhibition blocks growth of RB1-positive Triple Negative Breast Cancer" by Wu et. al. 2019
#
# Code written by: Wail Ba-alawi, email: wail.ba-alawi@uhnresearch.ca
#
###########################################
###########################################



# **** IMPORTANT NOTES ****
# 1. Once the project is downloaded to the user computer, the user needs to navigate to the main directory of the project "TNBC_BAY-876-master"
# 2. Set the working directory to "TNBC_BAY-876-master" or double-click on TNBC_GLUT1.Rproj and it will set it automatically
# 3. Make sure you have downloaded all the needed data files from the figshare url


library(Biobase)
library(PharmacoGx)
library(ggplot2)
library(GSA)
library(piano)
library(fgsea)
library(dplyr)
library(EnrichmentBrowser) # devtools::install_github("lgeistlinger/EnrichmentBrowser")
library(ggrepel)


source("R/plotting_functions.R")


load("data/metrics_for_sensitivity.RData")
#######################

load("data/SLC2A1_expression_across_TCGA_METABRIC_PMH.RData")
# Fig 1: a,b,c --- Fig S1: a,b,c

plot_gene_bySubtype(TCGA_data_SLC2A1,"TCGA","PAM50","SLC2A1")
plot_gene_bySubtype_specific(TCGA_data_SLC2A1,"TCGA","PAM50","Basal","SLC2A1")


plot_gene_bySubtype(Metabric_data_SLC2A1,"METABRIC","PAM50","SLC2A1")
plot_gene_bySubtype_specific(TCGA_data_SLC2A1,"TCGA","PAM50","Basal","SLC2A1")


plot_gene_bySubtype(PDX_data_SLC2A1,"PMH","PAM50","SLC2A1")
plot_gene_bySubtype_specific(PDX_data_SLC2A1,"PMH","PAM50","Basal","SLC2A1")

#######################


#######################
# Get IC50 values.
# Fig 1: d 

#######################

#######################
# Fig 3: a,b

# differential expression of most sensitive vs most resistant cell lines to Bay-876

mostSensitive <- c("HCC1806","HCC38","Hs-578-T")
mostResistant <- c("BT-549","MDA-MB-436","MDA-MB-468")


UHN_RNAseq_counts_info <- data.frame("samples"=c(mostSensitive,mostResistant),"condition"=c(rep("sensitive",length(mostSensitive)),rep("resistant",length(mostResistant))))
rownames(UHN_RNAseq_counts_info) <- UHN_RNAseq_counts_info$samples

load("data/PMCC_PSet.RData")

UHN_RNAseq_counts_sub <- summarizeMolecularProfiles(UHN,"rnaseq.counts",fill.missing = F,cell.lines = c(mostSensitive,mostResistant))
exprs(UHN_RNAseq_counts_sub) <- round(2^exprs(UHN_RNAseq_counts_sub)-1)

UHN_RNAseq_counts_sub$GROUP <- as.factor(ifelse(UHN_RNAseq_counts_info$condition=="sensitive",1,0))

deRes <- deAna(expr = UHN_RNAseq_counts_sub ,de.method = 'DESeq',stat.only = F,filter.by.expr=F)
deRes_final <- rowData(deRes)
deRes_final <- deRes_final[!is.na(deRes_final$ADJ.PVAL),]
rownames(deRes_final) <- deRes_final$gene_id


deRes_final2 <- as.data.frame(deRes_final)
mutateddf <- mutate(deRes_final2, sig=ifelse(deRes_final2$ADJ.PVAL<0.05, "FDR<0.05", "Not Sig")) #Will have different colors depending on significance
input <- cbind(gene=rownames(deRes_final2 ), mutateddf ) #convert the rownames to a column
input$FC <- -input$FC
mutateddf$FC <- -mutateddf$FC
volc = ggplot(mutateddf, aes(FC, -log10(PVAL))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance 
  scale_color_manual(values=c("red","black")) + 
  ggtitle("Fig 3a") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

volc 

tmp <- rowData(deRes)
tmp <- tmp[!is.na(tmp$ADJ.PVAL),]
genelevelstats <- tmp$DESeq2.STAT
names(genelevelstats) <- tmp$gene_id
genelevelstats <- genelevelstats[!is.na(genelevelstats)]
genelevelstats <- sort(genelevelstats,decreasing = T)
removeDuplicates <- function(IDs_source, IDs_target){
  
  ibx <- which(!is.na(IDs_target))
  idx <- which(!duplicated(IDs_target))
  
  return(intersect(ibx,idx))
  
}


load("data/genes_mappings.RData")

keepGenes <- removeDuplicates(names(genelevelstats),genes_mappings[names(genelevelstats),"Symbol"])

genelevelstats <- genelevelstats[keepGenes]
names(genelevelstats) <- genes_mappings[names(genelevelstats),"Symbol"]

genelevelstats <- sort(genelevelstats,decreasing = T)


gsc1 <- loadGSC("data/h.all.v6.1.symbols.gmt")

nPerm=1000000
gsea_out <- piano::runGSA(geneLevelStats=genelevelstats, geneSetStat="fgsea", gsc=gsc1,  nPerm=nPerm, ncpus=4, adjMethod="none", verbose=FALSE)

gseaResSummary <- piano::GSAsummaryTable(gsea_out)
gseares1 <- gseaResSummary[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE)
gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(nPerm+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
gseares1 <- gseares1[order(gseares1$pval,na.last = T),]
rownames(gseares1) <- gseares1$Name

geneSet <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
geneSet <- "HALLMARK_E2F_TARGETS"

geneSets <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_E2F_TARGETS")


geneSets_Associated_with_ResistanceToBAY876 <- gseares1[gseares1$`Stat (dist.dir)`<=0 & gseares1$FDR < 0.05,]

lapply(geneSets, function(geneSet){
  
  plotEnrichment(gsc1[[1]][[geneSet]],
                 genelevelstats) + labs(title=geneSet,  subtitle=paste("Estimate: " ,sprintf("%.3g", gseares1[geneSet,"Stat (dist.dir)"]),
                                                                       #", P-value: ",sprintf("%.1E",gseares1[geneSet,"pval"]),
                                                                       ", FDR: ",sprintf("%.1E",gseares1[geneSet,"FDR"])
                                                                       ,sep = "" )) +
    theme( panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,face = "bold"), plot.subtitle  = element_text(hjust = 0.5,face = "bold") 
           ,axis.title=element_text(face="bold"),axis.text=element_text(face="bold"),text = element_text(face = "bold",size = 20) )
  
  
})

#######################



#######################
# Fig 3: c

#######################


#######################
# Fig 3: d,e

# MD Anderson Proteomics associations

load("data/RPPA_MDAnderson.RData")
results <- t(apply(RPPA_MDAnderson, 2, function(x){
  if(sum(is.na(x))>2){
    return(c(0,1))
  }
  c <- cor.test(x,metrics_for_sensitivity[rownames(RPPA_MDAnderson),"IC50"])
  
  return(c(c$estimate,c$p.value))
}))

results <- cbind(results,p.adjust(results[,2]))
results <- results[order(results[,2],na.last = T,decreasing = F),]

results <- as.data.frame(results)
colnames(results) <- c("cor","PVAL","ADJ.PVAL")
mutateddf <- mutate(results, sig=ifelse(results$PVAL<0.05, "pval<0.05", "Not Sig")) #Will have different colors depending on significance
input <- cbind(gene=rownames(results ), mutateddf ) #convert the rownames to a column
volc = ggplot(input, aes(cor, -log10(PVAL))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance 
  scale_color_manual(values=c("black", "red")) + 
  ggtitle("MD Anderson BRCA RPPA [452 proteins - 12 TNBC cell lines]") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

volc + geom_text_repel(data=head(input, 8), aes(label=gene)) 



###
#PMCC Proteomics associations

load("data/RPPA_PMCC.RData")
commonCells <- intersect(rownames(metrics_for_sensitivity),rownames(RPPA_PMCC))


results <- t(apply(RPPA_PMCC[commonCells,], 2, function(x){
  if(sum(is.na(x))>2){
    return(c(0,1))
  }
  c <- cor.test(x,metrics_for_sensitivity[rownames(RPPA_PMCC[commonCells,]),"IC50"])
  
  return(c(c$estimate,c$p.value))
}))

results <- cbind(results,p.adjust(results[,2]))
results <- results[order(results[,2],na.last = T,decreasing = F),]

results <- as.data.frame(results)
colnames(results) <- c("cor","PVAL","ADJ.PVAL")
mutateddf <- mutate(results, sig=ifelse(results$PVAL<0.05, "pval<0.05", "Not Sig")) #Will have different colors depending on significance
input <- cbind(gene=rownames(results ), mutateddf,symbol=RPPA_PMCC_mapping[rownames(results )] ) #convert the rownames to a column
volc = ggplot(input, aes(cor, -log10(PVAL))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance 
  scale_color_manual(values=c("black", "red")) + 
  ggtitle("PMCC BRCA RPPA [218 proteins - 17 TNBC cell lines]") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

volc + geom_text_repel(data=head(input, 5), aes(label=gene)) 

##############################







