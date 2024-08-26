library(ggplot2)
library(ggsci)
library(readxl)
library(writexl)
library(Gmisc) 
library(dplyr)
library(patchwork)
library(cowplot)
library(ggVennDiagram)


# Define tau status
bf2.all.olink$tau_status <- ifelse(bf2.all.olink$tnic_cho_com_I_IV<1.36 ,0, 1)

## Define AT status based on CSF Ab ratio and tau-PET
bf2.all.olink$AT_status <- ifelse((bf2.all.olink$tau_status==1 & bf2.all.olink$Abnormal_CSF_Ab42_Ab40_Ratio==1),'A+T+',
                                  ifelse((bf2.all.olink$tau_status==0 & bf2.all.olink$Abnormal_CSF_Ab42_Ab40_Ratio==1),'A+T-',
                                         ifelse((bf2.all.olink$tau_status==0 & bf2.all.olink$Abnormal_CSF_Ab42_Ab40_Ratio==0),'A-T-','A-T+')))  

############ Group comparison analyses ###########

## Example comparing A-T- vs A+T- and then A+T- vs A+T+
tmp.df <- subset(bf2.all.olink,((AT_status=="A-T-" | AT_status=="A+T-") & !(study_cohort_baseline_variable=="BF2_E"))) ## Excluding cohort E
tmp.df <- subset(bf2.all.olink,((AT_status=="A+T-" | AT_status=="A+T+") & !(study_cohort_baseline_variable=="BF2_E"))) ## Excluding cohort E
tmp.df <- subset(bf2.all.olink,((AT_status=="A-T-" | AT_status=="A+T+") & !(study_cohort_baseline_variable=="BF2_E"))) ## Excluding cohort E

## select proteins 
id_lod <- which(list.lod<0.30)  ## 1331 proteins
list_to_keep <- list.all.olink2943[id_lod]

list.all.olink <- list_to_keep


olink.group.comp<-data.frame(matrix(nrow=1331,ncol=6))  ## keeping only proteins with 70% participants above LOD
colnames(olink.group.comp)=c('StdBeta_AT_comp_adjNPX','F_test_adjNPX','p_F_test_adjNPX',
                             'StdBeta_AT_comp_null','F_test_null','p_F_test_null')
olink.group.comp$olink_proteins<-list.all.olink

## Loop for group comparison
for (i in 1:length(list.all.olink)){
  svMisc::progress(i, max.value = length(list.all.olink))
  
  olink.protein <- list.all.olink[i]
  
  olink.name.tmp <- strsplit(olink.protein,split="__")[[1]][1]
  tmp <- tmp.df %>% select(contains(paste0(olink.name.tmp,"__")))
  
  names(tmp)<-c("current.olink","QC_Warning")
  
  # If values has WARN in QC, change for NA
  tmp$current.olink <- ifelse(tmp$QC_Warning=='WARN', NA, tmp$current.olink)
  
  ## Merge current data with tmp olink data
  tmp.tau<-cbind(tmp.df,tmp)
  
  m00 <- lm.beta(lm(current.olink ~  age + gender_baseline_variable + mean_NPX_zscore, data = tmp.tau))
  m11 <- lm.beta(lm(current.olink ~  AT_status + age + gender_baseline_variable + mean_NPX_zscore, data = tmp.tau))
  
  sum.lm11 <- summary(m11)
  f.test11 <- anova(m00,m11,test='F')
  
  m0null <- lm.beta(lm(current.olink ~  age + gender_baseline_variable, data = tmp.tau))
  m1null <- lm.beta(lm(current.olink ~  AT_status + age + gender_baseline_variable, data = tmp.tau))
  
  sum.lmnull <- summary(m1null)
  f.testnull <- anova(m0null,m1null,test='F')
  
  olink.group.comp[i,1]<-sum.lm11$coefficients[ 2, 2]  # Extract Std beta for olink protein
  olink.group.comp[i,2]<-f.test11$F[2]  # Save F-value
  olink.group.comp[i,3]<-f.test11$`Pr(>F)`[2]  # Save corresponding p-value
  
  olink.group.comp[i,4]<-sum.lmnull$coefficients[ 2, 2]  # Extract Std beta for olink protein
  olink.group.comp[i,5]<-f.testnull$F[2] 
  olink.group.comp[i,6]<-f.testnull$`Pr(>F)`[2]
}

olink.group.comp$p_fdr_adjNPX <- p.adjust(olink.group.comp$p_F_test_adjNPX,method="BH")
olink.group.comp$p_fdr_null <- p.adjust(olink.group.comp$p_F_test_null,method="BH")

### Volcano plots ####

library(ggrepel)
## Final plots 
df1$label <- ifelse(df1$p_fdr_adjNPX < 0.001 , df1$Gene_Symbol, "") 
df1$color <- ifelse((df1$p_fdr_adjNPX < 0.01 & df1$StdBeta_AT_comp_adjNPX < 0) , "down", 
                    ifelse((df1$p_fdr_adjNPX < 0.01 & df1$StdBeta_AT_comp_adjNPX > 0),"up","unchanged"))
p1 <- ggplot (data = df1, aes (x = StdBeta_AT_comp_adjNPX, y = -log10(p_fdr_adjNPX), col = color, label=label)) + geom_point () + 
  scale_color_manual(values = c("red","grey","blue"))+theme_minimal() + geom_text_repel(max.overlaps = Inf,size=4)
p1 <- p1 + geom_hline (yintercept=-log10(0.01), col="red", size=1) + xlab("Standardized beta")+ylab("log10(pFDR)") + ggtitle("A-T- vs A+T-")
p1


df2$label <- ifelse(df2$p_fdr_adjNPX < 0.001 , df2$Gene_Symbol, "")  
df2$color <- ifelse((df2$p_fdr_adjNPX < 0.01 & df2$StdBeta_AT_comp_adjNPX < 0) , "down", 
                    ifelse((df2$p_fdr_adjNPX < 0.01 & df2$StdBeta_AT_comp_adjNPX > 0),"up","unchanged"))
p2 <- ggplot (data = df2, aes (x = StdBeta_AT_comp_adjNPX, y = -log10(p_fdr_adjNPX), col = color, label=label)) + geom_point () + 
  scale_color_manual(values = c("red","grey","blue"))+theme_minimal() + geom_text_repel(max.overlaps = Inf,size=4)
p2 <- p2 + geom_hline (yintercept=-log10(0.01), col="red", size=1) + xlab("Standardized beta")+ylab("log10(pFDR)")+ ggtitle("A+T- vs A+T+")
p2


df3$label <- ifelse(df3$p_fdr_adjNPX < 0.000001 , df3$Gene_Symbol, "")  
df3$color <- ifelse((df3$p_fdr_adjNPX < 0.01 & df3$StdBeta_AT_comp_adjNPX < 0) , "down", 
                    ifelse((df3$p_fdr_adjNPX < 0.01 & df3$StdBeta_AT_comp_adjNPX > 0),"up","unchanged"))
p3 <- ggplot (data = df3, aes (x = -1*StdBeta_AT_comp_adjNPX, y = -log10(p_fdr_adjNPX), col = color, label=label)) + geom_point () + 
  scale_color_manual(values = c("red","grey","blue"))+theme_minimal() + geom_text_repel(max.overlaps = Inf,size=4)
p3 <- p3 + geom_hline (yintercept=-log10(0.01), col="red", size=1) + xlab("Standardized beta")+ylab("log10(pFDR)")+ ggtitle("A+T+ vs non-AD Ab-negative")
p3

df4$label <- ifelse(df4$p_fdr_adjNPX < 0.0001 , df4$Gene_Symbol, "")  
df4$color <- ifelse((df4$p_fdr_adjNPX < 0.01 & df4$StdBeta_AT_comp_adjNPX < 0) , "down", 
                    ifelse((df4$p_fdr_adjNPX < 0.01 & df4$StdBeta_AT_comp_adjNPX > 0),"up","unchanged"))
p4 <- ggplot (data = df4, aes (x = StdBeta_AT_comp_adjNPX, y = -log10(p_fdr_adjNPX), col = color, label=label)) + geom_point () + 
  scale_color_manual(values = c("red","grey","blue"))+theme_minimal() + geom_text_repel(max.overlaps = Inf,size=4)
p4 <- p4 + geom_hline (yintercept=-log10(0.01), col="red", size=1) + xlab("Standardized beta")+ylab("log10(pFDR)")+ ggtitle("A-T- vs non-AD Ab-negative")
p4

ggpubr::ggarrange(p1,p2,p3,p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

#### Venn diagrams ######
venn <- list(ATmin = ATmin$olink_proteins,
             AposTmin = AposTmin$olink_proteins,
             Pos_vs_nonAD = Pos_vs_nonAD$olink_proteins,
             Neg_vs_nonAD = Neg_vs_nonAD$olink_proteins
)

venn1 <- ggVennDiagram(venn[1:4], label_alpha = 0,
                       category.names = c("A-T- vs. A+T-","A+T- vs. A+T+","A+T+ vs. non-AD",
                                          "A-T- vs. non-AD"))+
  ggplot2::scale_fill_gradient(low="white",high = "red")+
  scale_color_brewer(palette = "Set1") ## Change bordder 


venn2 <- ggVennDiagram(venn[1:2], label_alpha = 0,
                       category.names = c("A-T- vs. A+T-","A+T- vs. A+T+"))+
  ggplot2::scale_fill_gradient(low="white",high = "red")+
  scale_color_brewer(palette = "Set1") ## Change bordder 