#!/usr/bin/Rscript

################################################################################
## Title:         epigenomics_preprocessing.R
## Description:   Imputation and normalization of epigenomics data
## Author:        Casper de Visser, Jenny van Dongen
## Date created:  2022-03-02
## Email:         casper.devisser@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################
memory.limit(size = 1500000)

#Read in data

args = commandArgs(trailingOnly=TRUE)

Epigenomics_MOESM1_path = args[1]
Epigenomics_MOESM4_path = args[2]
Epigenomics_MOESM5_path = args[3]
MethylationEPIC_path = args[4]
anno_epic_path = args[5]

output_dir = args[6]



###############################################################
######################## EPIC ANNOTATION ######################
###############################################################




options(stringsAsFactors=FALSE)
library(meffil)
meffil.list.featuresets()
# "450k"   "common" "epic"
epicanno <-meffil.get.features("epic")
commonanno <-meffil.get.features("common")

# cross-reactive probes
pidsley       <- read.csv(Epigenomics_MOESM1_path)  
#head(pidsley)
#table(pidsley$Total)
crossreactive <-  pidsley$X
#length(crossreactive)


# Probes overlapping genetic variants at targeted CpG sites
pidsley       <- read.csv(Epigenomics_MOESM4_path) 
#head(pidsley)
#dim(pidsley[which(pidsley$EUR_AF > 0.01),])
DNApolymorphismtarget <- pidsley[which(pidsley$EUR_AF > 0.01),"PROBE"] 


# probes with SNP in single base extension site
pidsley       <- read.csv(Epigenomics_MOESM5_path) 
#head(pidsley)
#dim(pidsley[which(pidsley$EUR_AF > 0.01),])
DNApolymorphismSBE <- pidsley[which(pidsley$EUR_AF > 0.01),"PROBE"] 



#length(intersect(DNApolymorphismSBE,DNApolymorphismtarget))
#length(unique(c(DNApolymorphismSBE,DNApolymorphismtarget,crossreactive))

# Illumina annotation - B4
a <- read.csv(MethylationEPIC_path) 
#dim(a)
#table(a$IlmnID==a$Name)
#a[which(!a$IlmnID==a$Name),c("IlmnID","Name")]
# remove control probes
a <- a[1:865918,]
#dim(a)      # 865918
#table(a$IlmnID==a$Name)
rownames(a) <- a$Name
 
load(anno_epic_path) 
#dim(epicanno)     # 866553 
#intersect(a$Name,epicanno$name)  # 865918
#head(epicanno)
rownames(epicanno) <- epicanno$name
epicanno <- epicanno[rownames(a),]

# check
### plot1 <- plot(a$MAPINFO,epicanno$position) #identical


# Phantom
#table(a$Phantom4_Enhancers)
#table(a$Phantom5_Enhancers)
a$Phantom4_Enhancer <- rep (0,nrow(a))
a$Phantom4_Enhancer[which(!a$Phantom4_Enhancers == "")]    <- 1 
a$Phantom5_Enhancer <- rep (0,nrow(a))
a$Phantom5_Enhancer[which(!a$Phantom5_Enhancers == "")]    <- 1
#table(a$Phantom4_Enhancer)
#table(a$Phantom5_Enhancer)
#table(a$Phantom4_Enhancer,a$Phantom5_Enhancer)      

#ENCODE
#table(a$Regulatory_Feature_Name)
#table(a$Regulatory_Feature_Group)  # not very useful
#table(a$DNase_Hypersensitivity_Evidence_Count)       #ENCODE
#length(which(!a$DNase_Hypersensitivity_NAME== "")) # 497025
#table(a$OpenChromatin_Evidence_Count)
#sum(table(a$OpenChromatin_Evidence_Count)) #119053
#table(a$TFBS_Evidence_Count)
#sum(table(a$TFBS_Evidence_Count))
#length(which(!a$TFBS_NAME== ""))  #134170 

#table(a$GencodeCompV12_Group)
#table(a$DMR)
#table(a$Random_Loci)
#table(a$X)

sel <- c("IlmnID", "Genome_Build", "CHR", "MAPINFO", "Relation_to_UCSC_CpG_Island","UCSC_RefGene_Group","Phantom4_Enhancer","Phantom4_Enhancers","Phantom5_Enhancer","Phantom5_Enhancers","Regulatory_Feature_Name", "Regulatory_Feature_Group","DNase_Hypersensitivity_Evidence_Count","DNase_Hypersensitivity_NAME","OpenChromatin_Evidence_Count","TFBS_Evidence_Count","TFBS_NAME","DMR")
a <- a[,sel]
#TODO: save this file as result?
### save(a, file="anno_epic_genomicfeatures.RData") 


save(epicanno, crossreactive, DNApolymorphismtarget, DNApolymorphismSBE, commonanno, file = output_dir)
