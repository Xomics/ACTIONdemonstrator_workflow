###############################################################
####################### PROBE FILTERING #######################
###############################################################
memory.limit(size = 1500000)

args = commandArgs(trailingOnly=TRUE)

epigenomics_values_path = args[1]
annotation_out_path = args[2]

output_dir = args[3]
#output_dir2 = args[4]

# Load data
betas <- read.csv(epigenomics_values_path, row.names=1) # beta_values
load(annotation_out_path)

#dim(betas) # 857793 

# Probes are excluded from all samples if they mapped to multiple locations in the genome, or if they overlapped with a single nucleotide polymorphism (SNP) or Insertion/Deletion (INDEL). Annotations of ambiguous mapping probes (based on an overlap of at least 47 bases per probe) and probes where genetic variants (SNPs or INDELS) with a minor allele frequency>0.01 in Europeans overlap with the targeted CpG or single base extension site (SBE) were obtained from Pidsley R, Zotenko E, Peters TJ, Lawrence MG, Risbridger GP, Molloy P, et al. Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling. Genome Biol. 2016;17. 

#length(intersect(rownames(betas),epicanno$name))       # note that the epic data contain some probes that are no longer in MethylationEPIC_v-1-0_B4
OL <- intersect(rownames(betas),epicanno$name)
rownames(epicanno) <- epicanno$name
epicanno <- epicanno[OL,]
betas <- betas[OL,]
#nrow(betas)


# 1. removal of cross-reactive probes and probes overlapping SNPs
rm <- c(crossreactive, DNApolymorphismtarget,DNApolymorphismSBE)
#length(which(rownames(betas) %in% rm))     # 52342
betas <- betas[which(!rownames(betas) %in% rm),]
#dim(betas)     #    805235 

# 2. removal of XY chromosomes
epicanno <- epicanno[rownames(betas),]
#table(rownames(epicanno)==rownames(betas))
#table(epicanno$chromosome)
rm <- which(epicanno$chromosome=="chrX" | epicanno$chromosome=="chrY")
epicanno <- epicanno[-rm,]
betas <- betas[rownames(epicanno),]
#dim(betas)   #  787711 

# 3. annotation of EPIC versus IL450k probes
#dim(commonanno)
#commonprobes <- intersect(rownames(betas), commonanno$name)
#EPICprobes   <- rownames(betas)[which(!rownames(betas) %in% commonprobes)]


# Save files
#save(betas, file = output_dir)
write.csv(betas, output_dir)
#save(EPICprobes, file = output_dir3)
