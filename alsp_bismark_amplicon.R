library(ggplot2)
library(ChAMP) #Just to load right stack
setwd("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/allele_sp/")
SNPa <- read.table("/home/ankitv/ref_av/Snpfiles/custom_commonSNPrefalt_GRCh38_Acomb_Snpfile.txt")
SNPb <- read.table("/home/ankitv/ref_av/Snpfiles/custom_commonSNPrefalt_GRCh38_Bcomb_Snpfile.txt")

S49_covA_ref <- read.table("49_S49_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S49_covA_ref) <- paste0(S49_covA_ref$V1, "_", S49_covA_ref$V2)
S49_covA_ref["Group"] <- "S49_covA_ref"
dim(S49_covA_ref)
colnames(S49_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S49_covA_ref["logcount"] <- log2(S49_covA_ref$M + S49_covA_ref$U)
S49_covA_alt <- read.table("49_S49_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S49_covA_alt) <- paste0(S49_covA_alt$V1, "_", S49_covA_alt$V2)
S49_covA_alt["Group"] <- "S49_covA_alt"
dim(S49_covA_alt)
colnames(S49_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S49_covA_alt["logcount"] <- log2(S49_covA_alt$M + S49_covA_alt$U)
library(ggplot2)
S49_covA <- rbind.data.frame(S49_covA_ref,S49_covA_alt)
ggplot(S49_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S49_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S49_covA.svg", width=15, height=8, units="cm", dpi=96)


S49_covB_ref <- read.table("49_S49_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S49_covB_ref) <- paste0(S49_covB_ref$V1, "_", S49_covB_ref$V2)
S49_covB_ref["Group"] <- "S49_covB_ref"
dim(S49_covB_ref)
colnames(S49_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S49_covB_ref["logcount"] <- log2(S49_covB_ref$M + S49_covB_ref$U)
S49_covB_alt <- read.table("49_S49_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S49_covB_alt) <- paste0(S49_covB_alt$V1, "_", S49_covB_alt$V2)
S49_covB_alt["Group"] <- "S49_covB_alt"
dim(S49_covB_alt)
colnames(S49_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S49_covB_alt["logcount"] <- log2(S49_covB_alt$M + S49_covB_alt$U)

S49_covB <- rbind.data.frame(S49_covB_ref,S49_covB_alt)
ggplot(S49_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("brown","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S49_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S49_covB.svg", width=15, height=8, units="cm", dpi=96)

S50_covA_ref <- read.table("50_S50_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S50_covA_ref) <- paste0(S50_covA_ref$V1, "_", S50_covA_ref$V2)
S50_covA_ref["Group"] <- "S50_covA_ref"
dim(S50_covA_ref)
colnames(S50_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S50_covA_ref["logcount"] <- log2(S50_covA_ref$M + S50_covA_ref$U)
S50_covA_alt <- read.table("50_S50_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S50_covA_alt) <- paste0(S50_covA_alt$V1, "_", S50_covA_alt$V2)
S50_covA_alt["Group"] <- "S50_covA_alt"
dim(S50_covA_alt)
colnames(S50_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S50_covA_alt["logcount"] <- log2(S50_covA_alt$M + S50_covA_alt$U)

S50_covA <- rbind.data.frame(S50_covA_ref,S50_covA_alt)
ggplot(S50_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S50_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S50_covA.svg", width=15, height=8, units="cm", dpi=96)


S50_covB_ref <- read.table("50_S50_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S50_covB_ref) <- paste0(S50_covB_ref$V1, "_", S50_covB_ref$V2)
S50_covB_ref["Group"] <- "S50_covB_ref"
dim(S50_covB_ref)
colnames(S50_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S50_covB_ref["logcount"] <- log2(S50_covB_ref$M + S50_covB_ref$U)
S50_covB_alt <- read.table("50_S50_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S50_covB_alt) <- paste0(S50_covB_alt$V1, "_", S50_covB_alt$V2)
S50_covB_alt["Group"] <- "S50_covB_alt"
dim(S50_covB_alt)
colnames(S50_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S50_covB_alt["logcount"] <- log2(S50_covB_alt$M + S50_covB_alt$U)

S50_covB <- rbind.data.frame(S50_covB_ref,S50_covB_alt)
ggplot(S50_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("brown","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S50_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S50_covB.svg", width=15, height=8, units="cm", dpi=96)

S51_covA_ref <- read.table("51_S51_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S51_covA_ref) <- paste0(S51_covA_ref$V1, "_", S51_covA_ref$V2)
S51_covA_ref["Group"] <- "S51_covA_ref"
dim(S51_covA_ref)
colnames(S51_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S51_covA_ref["logcount"] <- log2(S51_covA_ref$M + S51_covA_ref$U)
S51_covA_alt <- read.table("51_S51_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S51_covA_alt) <- paste0(S51_covA_alt$V1, "_", S51_covA_alt$V2)
S51_covA_alt["Group"] <- "S51_covA_alt"
dim(S51_covA_alt)
colnames(S51_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S51_covA_alt["logcount"] <- log2(S51_covA_alt$M + S51_covA_alt$U)

S51_covA <- rbind.data.frame(S51_covA_ref,S51_covA_alt)
ggplot(S51_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S51_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S51_covA.svg", width=15, height=8, units="cm", dpi=96)


S51_covB_ref <- read.table("51_S51_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S51_covB_ref) <- paste0(S51_covB_ref$V1, "_", S51_covB_ref$V2)
S51_covB_ref["Group"] <- "S51_covB_ref"
dim(S51_covB_ref)
colnames(S51_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S51_covB_ref["logcount"] <- log2(S51_covB_ref$M + S51_covB_ref$U)
S51_covB_alt <- read.table("51_S51_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S51_covB_alt) <- paste0(S51_covB_alt$V1, "_", S51_covB_alt$V2)
S51_covB_alt["Group"] <- "S51_covB_alt"
dim(S51_covB_alt)
colnames(S51_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S51_covB_alt["logcount"] <- log2(S51_covB_alt$M + S51_covB_alt$U)

S51_covB <- rbind.data.frame(S51_covB_ref,S51_covB_alt)
ggplot(S51_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("brown","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S51_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S51_covB.svg", width=15, height=8, units="cm", dpi=96)


S52_covA_ref <- read.table("52_S52_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S52_covA_ref) <- paste0(S52_covA_ref$V1, "_", S52_covA_ref$V2)
S52_covA_ref["Group"] <- "S52_covA_ref"
dim(S52_covA_ref)
colnames(S52_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S52_covA_ref["logcount"] <- log2(S52_covA_ref$M + S52_covA_ref$U)
S52_covA_alt <- read.table("52_S52_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S52_covA_alt) <- paste0(S52_covA_alt$V1, "_", S52_covA_alt$V2)
S52_covA_alt["Group"] <- "S52_covA_alt"
dim(S52_covA_alt)
colnames(S52_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S52_covA_alt["logcount"] <- log2(S52_covA_alt$M + S52_covA_alt$U)

S52_covA <- rbind.data.frame(S52_covA_ref,S52_covA_alt)
ggplot(S52_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S52_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S52_covA.svg", width=15, height=8, units="cm", dpi=96)


S52_covB_ref <- read.table("52_S52_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S52_covB_ref) <- paste0(S52_covB_ref$V1, "_", S52_covB_ref$V2)
S52_covB_ref["Group"] <- "S52_covB_ref"
dim(S52_covB_ref)
colnames(S52_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S52_covB_ref["logcount"] <- log2(S52_covB_ref$M + S52_covB_ref$U)
S52_covB_alt <- read.table("52_S52_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S52_covB_alt) <- paste0(S52_covB_alt$V1, "_", S52_covB_alt$V2)
S52_covB_alt["Group"] <- "S52_covB_alt"
dim(S52_covB_alt)
colnames(S52_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S52_covB_alt["logcount"] <- log2(S52_covB_alt$M + S52_covB_alt$U)

S52_covB <- rbind.data.frame(S52_covB_ref,S52_covB_alt)
ggplot(S52_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("brown","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S52_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S52_covB.svg", width=15, height=8, units="cm", dpi=96)

S53_covA_ref <- read.table("53_S53_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S53_covA_ref) <- paste0(S53_covA_ref$V1, "_", S53_covA_ref$V2)
S53_covA_ref["Group"] <- "S53_covA_ref"
dim(S53_covA_ref)
colnames(S53_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S53_covA_ref["logcount"] <- log2(S53_covA_ref$M + S53_covA_ref$U)
S53_covA_alt <- read.table("53_S53_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S53_covA_alt) <- paste0(S53_covA_alt$V1, "_", S53_covA_alt$V2)
S53_covA_alt["Group"] <- "S53_covA_alt"
dim(S53_covA_alt)
colnames(S53_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S53_covA_alt["logcount"] <- log2(S53_covA_alt$M + S53_covA_alt$U)

S53_covA <- rbind.data.frame(S53_covA_ref,S53_covA_alt)
ggplot(S53_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S53_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S53_covA.svg", width=15, height=8, units="cm", dpi=96)


S53_covB_ref <- read.table("53_S53_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S53_covB_ref) <- paste0(S53_covB_ref$V1, "_", S53_covB_ref$V2)
S53_covB_ref["Group"] <- "S53_covB_ref"
dim(S53_covB_ref)
colnames(S53_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S53_covB_ref["logcount"] <- log2(S53_covB_ref$M + S53_covB_ref$U)
S53_covB_alt <- read.table("53_S53_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S53_covB_alt) <- paste0(S53_covB_alt$V1, "_", S53_covB_alt$V2)
S53_covB_alt["Group"] <- "S53_covB_alt"
dim(S53_covB_alt)
colnames(S53_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S53_covB_alt["logcount"] <- log2(S53_covB_alt$M + S53_covB_alt$U)

S53_covB <- rbind.data.frame(S53_covB_ref,S53_covB_alt)
ggplot(S53_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("brown","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S53_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S53_covB.svg", width=15, height=8, units="cm", dpi=96)



S54_covA_ref <- read.table("54_S54_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S54_covA_ref) <- paste0(S54_covA_ref$V1, "_", S54_covA_ref$V2)
S54_covA_ref["Group"] <- "S54_covA_ref"
dim(S54_covA_ref)
colnames(S54_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S54_covA_ref["logcount"] <- log2(S54_covA_ref$M + S54_covA_ref$U)
S54_covA_alt <- read.table("54_S54_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S54_covA_alt) <- paste0(S54_covA_alt$V1, "_", S54_covA_alt$V2)
S54_covA_alt["Group"] <- "S54_covA_alt"
dim(S54_covA_alt)
colnames(S54_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S54_covA_alt["logcount"] <- log2(S54_covA_alt$M + S54_covA_alt$U)

S54_covA <- rbind.data.frame(S54_covA_ref,S54_covA_alt)
ggplot(S54_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("red","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S54_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S54_covA.svg", width=15, height=8, units="cm", dpi=96)


S54_covB_ref <- read.table("54_S54_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S54_covB_ref) <- paste0(S54_covB_ref$V1, "_", S54_covB_ref$V2)
S54_covB_ref["Group"] <- "S54_covB_ref"
dim(S54_covB_ref)
colnames(S54_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S54_covB_ref["logcount"] <- log2(S54_covB_ref$M + S54_covB_ref$U)
S54_covB_alt <- read.table("54_S54_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S54_covB_alt) <- paste0(S54_covB_alt$V1, "_", S54_covB_alt$V2)
S54_covB_alt["Group"] <- "S54_covB_alt"
dim(S54_covB_alt)
colnames(S54_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S54_covB_alt["logcount"] <- log2(S54_covB_alt$M + S54_covB_alt$U)

S54_covB <- rbind.data.frame(S54_covB_ref,S54_covB_alt)
ggplot(S54_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","red"))+ylim(0,100)+
  theme_bw()+ggtitle("S54_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S54_covB.svg", width=15, height=8, units="cm", dpi=96)

S55_covA_ref <- read.table("55_S55_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S55_covA_ref) <- paste0(S55_covA_ref$V1, "_", S55_covA_ref$V2)
S55_covA_ref["Group"] <- "S55_covA_ref"
dim(S55_covA_ref)
colnames(S55_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S55_covA_ref["logcount"] <- log2(S55_covA_ref$M + S55_covA_ref$U)
S55_covA_alt <- read.table("55_S55_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S55_covA_alt) <- paste0(S55_covA_alt$V1, "_", S55_covA_alt$V2)
S55_covA_alt["Group"] <- "S55_covA_alt"
dim(S55_covA_alt)
colnames(S55_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S55_covA_alt["logcount"] <- log2(S55_covA_alt$M + S55_covA_alt$U)

S55_covA <- rbind.data.frame(S55_covA_ref,S55_covA_alt)
S55_covA <- rbind.data.frame(S55_covA_ref)

ggplot(S55_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("red","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S55_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S55_covA.svg", width=15, height=8, units="cm", dpi=96)


S55_covB_ref <- read.table("55_S55_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S55_covB_ref) <- paste0(S55_covB_ref$V1, "_", S55_covB_ref$V2)
S55_covB_ref["Group"] <- "S55_covB_ref"
dim(S55_covB_ref)
colnames(S55_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S55_covB_ref["logcount"] <- log2(S55_covB_ref$M + S55_covB_ref$U)
S55_covB_alt <- read.table("55_S55_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S55_covB_alt) <- paste0(S55_covB_alt$V1, "_", S55_covB_alt$V2)
S55_covB_alt["Group"] <- "S55_covB_alt"
dim(S55_covB_alt)
colnames(S55_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S55_covB_alt["logcount"] <- log2(S55_covB_alt$M + S55_covB_alt$U)

S55_covB <- rbind.data.frame(S55_covB_ref,S55_covB_alt)
ggplot(S55_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","red"))+ylim(0,100)+
  theme_bw()+ggtitle("S55_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S55_covB.svg", width=15, height=8, units="cm", dpi=96)

S56_covA_ref <- read.table("56_S56_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S56_covA_ref) <- paste0(S56_covA_ref$V1, "_", S56_covA_ref$V2)
S56_covA_ref["Group"] <- "S56_covA_ref"
dim(S56_covA_ref)
colnames(S56_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S56_covA_ref["logcount"] <- log2(S56_covA_ref$M + S56_covA_ref$U)
S56_covA_alt <- read.table("56_S56_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S56_covA_alt) <- paste0(S56_covA_alt$V1, "_", S56_covA_alt$V2)
S56_covA_alt["Group"] <- "S56_covA_alt"
dim(S56_covA_alt)
colnames(S56_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S56_covA_alt["logcount"] <- log2(S56_covA_alt$M + S56_covA_alt$U)

S56_covA <- rbind.data.frame(S56_covA_ref,S56_covA_alt)
S56_covA <- rbind.data.frame(S56_covA_ref)

ggplot(S56_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("red","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S56_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S56_covA.svg", width=15, height=8, units="cm", dpi=96)


S56_covB_ref <- read.table("56_S56_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S56_covB_ref) <- paste0(S56_covB_ref$V1, "_", S56_covB_ref$V2)
S56_covB_ref["Group"] <- "S56_covB_ref"
dim(S56_covB_ref)
colnames(S56_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S56_covB_ref["logcount"] <- log2(S56_covB_ref$M + S56_covB_ref$U)
S56_covB_alt <- read.table("56_S56_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S56_covB_alt) <- paste0(S56_covB_alt$V1, "_", S56_covB_alt$V2)
S56_covB_alt["Group"] <- "S56_covB_alt"
dim(S56_covB_alt)
colnames(S56_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S56_covB_alt["logcount"] <- log2(S56_covB_alt$M + S56_covB_alt$U)

S56_covB <- rbind.data.frame(S56_covB_ref,S56_covB_alt)
ggplot(S56_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","red"))+ylim(0,100)+
  theme_bw()+ggtitle("S56_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S56_covB.svg", width=15, height=8, units="cm", dpi=96)


S57_covA_ref <- read.table("57_S57_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S57_covA_ref) <- paste0(S57_covA_ref$V1, "_", S57_covA_ref$V2)
S57_covA_ref["Group"] <- "S57_covA_ref"
dim(S57_covA_ref)
colnames(S57_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S57_covA_ref["logcount"] <- log2(S57_covA_ref$M + S57_covA_ref$U)
S57_covA_alt <- read.table("57_S57_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S57_covA_alt) <- paste0(S57_covA_alt$V1, "_", S57_covA_alt$V2)
S57_covA_alt["Group"] <- "S57_covA_alt"
dim(S57_covA_alt)
colnames(S57_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S57_covA_alt["logcount"] <- log2(S57_covA_alt$M + S57_covA_alt$U)

S57_covA <- rbind.data.frame(S57_covA_ref,S57_covA_alt)
S57_covA <- rbind.data.frame(S57_covA_ref)

ggplot(S57_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("red","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S57_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S57_covA.svg", width=15, height=8, units="cm", dpi=96)


S57_covB_ref <- read.table("57_S57_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S57_covB_ref) <- paste0(S57_covB_ref$V1, "_", S57_covB_ref$V2)
S57_covB_ref["Group"] <- "S57_covB_ref"
dim(S57_covB_ref)
colnames(S57_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S57_covB_ref["logcount"] <- log2(S57_covB_ref$M + S57_covB_ref$U)
S57_covB_alt <- read.table("57_S57_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S57_covB_alt) <- paste0(S57_covB_alt$V1, "_", S57_covB_alt$V2)
S57_covB_alt["Group"] <- "S57_covB_alt"
dim(S57_covB_alt)
colnames(S57_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S57_covB_alt["logcount"] <- log2(S57_covB_alt$M + S57_covB_alt$U)

S57_covB <- rbind.data.frame(S57_covB_ref,S57_covB_alt)
ggplot(S57_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","red"))+ylim(0,100)+
  theme_bw()+ggtitle("S57_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S57_covB.svg", width=15, height=8, units="cm", dpi=96)

S58_covA_ref <- read.table("58_S58_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov")
rownames(S58_covA_ref) <- paste0(S58_covA_ref$V1, "_", S58_covA_ref$V2)
S58_covA_ref["Group"] <- "S58_covA_ref"
dim(S58_covA_ref)
colnames(S58_covA_ref) <- c("chr","start","end","cov","M","U","Group")
S58_covA_ref["logcount"] <- log2(S58_covA_ref$M + S58_covA_ref$U)
S58_covA_alt <- read.table("58_S58_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov")
rownames(S58_covA_alt) <- paste0(S58_covA_alt$V1, "_", S58_covA_alt$V2)
S58_covA_alt["Group"] <- "S58_covA_alt"
dim(S58_covA_alt)
colnames(S58_covA_alt) <- c("chr","start","end","cov","M","U","Group")
S58_covA_alt["logcount"] <- log2(S58_covA_alt$M + S58_covA_alt$U)

S58_covA <- rbind.data.frame(S58_covA_ref,S58_covA_alt)
S58_covA <- rbind.data.frame(S58_covA_ref)

ggplot(S58_covA, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("red","green"))+ylim(0,100)+
  theme_bw()+ggtitle("S58_covA") + 
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S58_covA.svg", width=15, height=8, units="cm", dpi=96)


S58_covB_ref <- read.table("58_S58_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov")
rownames(S58_covB_ref) <- paste0(S58_covB_ref$V1, "_", S58_covB_ref$V2)
S58_covB_ref["Group"] <- "S58_covB_ref"
dim(S58_covB_ref)
colnames(S58_covB_ref) <- c("chr","start","end","cov","M","U","Group")
S58_covB_ref["logcount"] <- log2(S58_covB_ref$M + S58_covB_ref$U)
S58_covB_alt <- read.table("58_S58_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov")
rownames(S58_covB_alt) <- paste0(S58_covB_alt$V1, "_", S58_covB_alt$V2)
S58_covB_alt["Group"] <- "S58_covB_alt"
dim(S58_covB_alt)
colnames(S58_covB_alt) <- c("chr","start","end","cov","M","U","Group")
S58_covB_alt["logcount"] <- log2(S58_covB_alt$M + S58_covB_alt$U)

S58_covB <- rbind.data.frame(S58_covB_ref,S58_covB_alt)
ggplot(S58_covB, aes(x=as.numeric(start), y=cov))+
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group, size=logcount), shape=20)+
  scale_color_manual(values=c("blue","red"))+ylim(0,100)+
  theme_bw()+ggtitle("S58_covB")+
  scale_size(range = c(0,5), breaks = c(0,2,4,6,8,10), limits = c(0,NA))
ggsave("S58_covB.svg", width=15, height=8, units="cm", dpi=96)

#Prep for seqmonk
bedtools intersect -wa -wb -a myImpLoci.txt -b CG_hg38.txt | awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$11"\t"$6"."$7"."$8}' > myImpLoci_CGs_re.txt

#-------------------------------------------   END OF ANALYSIS ------------------------------------------#

#----------------
S54_covA <- read.table("54_S54_L001_R1_001_Acomb.coverage.txt")
rownames(S54_covA) <- paste0(S54_covA$V1, "_", S54_covA$V2)
S54_covA <- S54_covA[,c(5,6)]
colnames(S54_covA) <- c("T.M","T.U")
S54_covA["mePerc.T"] <- (S54_covA$T.M *100)/(S54_covA$T.M + S54_covA$T.U)
#S54_covA["mePerc.A"] <- (S54_covA$A.M *100)/(S54_covA$A.M + S54_covA$A.U)
S54_covB <- read.table("54_S54_L001_R1_001_Bcomb.coverage.txt")
rownames(S54_covB) <- paste0(S54_covB$V1, "_", S54_covB$V2)
S54_covB <- S54_covB[,c(5,6,11,12)]
colnames(S54_covB) <- c("T.M","T.U","C.M","C.U")
S54_covB["mePerc.T"] <- (S54_covB$T.M *100)/(S54_covB$T.M + S54_covB$T.U)
S54_covB["mePerc.C"] <- (S54_covB$C.M *100)/(S54_covB$C.M + S54_covB$C.U)

S55_covA <- read.table("55_S55_L001_R1_001_Acomb.coverage.txt")
rownames(S55_covA) <- paste0(S55_covA$V1, "_", S55_covA$V2)
S55_covA <- S55_covA[,c(5,6,11,12)]
colnames(S55_covA) <- c("T.M","T.U","A.M","A.U")
S55_covA["mePerc.T"] <- (S55_covA$T.M *100)/(S55_covA$T.M + S55_covA$T.U)
S55_covA["mePerc.A"] <- (S55_covA$A.M *100)/(S55_covA$A.M + S55_covA$A.U)
S55_covB <- read.table("55_S55_L001_R1_001_Bcomb.coverage.txt")
rownames(S55_covB) <- paste0(S55_covB$V1, "_", S55_covB$V2)
S55_covB <- S55_covB[,c(5,6,11,12)]
colnames(S55_covB) <- c("T.M","T.U","C.M","C.U")
S55_covB["mePerc.T"] <- (S55_covB$T.M *100)/(S55_covB$T.M + S55_covB$T.U)
S55_covB["mePerc.C"] <- (S55_covB$C.M *100)/(S55_covB$C.M + S55_covB$C.U)

S56_covA <- read.table("56_S56_L001_R1_001_Acomb.coverage.txt")
rownames(S56_covA) <- paste0(S56_covA$V1, "_", S56_covA$V2)
S56_covA <- S56_covA[,c(5,6,11,12)]
colnames(S56_covA) <- c("T.M","T.U","A.M","A.U")
S56_covA["mePerc.T"] <- (S56_covA$T.M *100)/(S56_covA$T.M + S56_covA$T.U)
S56_covA["mePerc.A"] <- (S56_covA$A.M *100)/(S56_covA$A.M + S56_covA$A.U)
S56_covB <- read.table("56_S56_L001_R1_001_Bcomb.coverage.txt")
rownames(S56_covB) <- paste0(S56_covB$V1, "_", S56_covB$V2)
S56_covB <- S56_covB[,c(5,6,11,12)]
colnames(S56_covB) <- c("T.M","T.U","C.M","C.U")
S56_covB["mePerc.T"] <- (S56_covB$T.M *100)/(S56_covB$T.M + S56_covB$T.U)
S56_covB["mePerc.C"] <- (S56_covB$C.M *100)/(S56_covB$C.M + S56_covB$C.U)

S57_covA <- read.table("57_S57_L001_R1_001_Acomb.coverage.txt")
rownames(S57_covA) <- paste0(S57_covA$V1, "_", S57_covA$V2)
S57_covA <- S57_covA[,c(5,6)]
#colnames(S57_covA) <- c("T.M","T.U","A.M","A.U")
S57_covA["mePerc.T"] <- (S57_covA$T.M *100)/(S57_covA$T.M + S57_covA$T.U)
S57_covA["mePerc.A"] <- (S57_covA$A.M *100)/(S57_covA$A.M + S57_covA$A.U)
S57_covB <- read.table("57_S57_L001_R1_001_Bcomb.coverage.txt")
rownames(S57_covB) <- paste0(S57_covB$V1, "_", S57_covB$V2)
S57_covB <- S57_covB[,c(5,6,11,12)]
colnames(S57_covB) <- c("T.M","T.U","C.M","C.U")
S57_covB["mePerc.T"] <- (S57_covB$T.M *100)/(S57_covB$T.M + S57_covB$T.U)
S57_covB["mePerc.C"] <- (S57_covB$C.M *100)/(S57_covB$C.M + S57_covB$C.U)

S58_covA <- read.table("58_S58_L001_R1_001_Acomb.coverage.txt")
rownames(S58_covA) <- paste0(S58_covA$V1, "_", S58_covA$V2)
S58_covA <- S58_covA[,c(5,6,11,12)]
colnames(S58_covA) <- c("T.M","T.U","A.M","A.U")
S58_covA["mePerc.T"] <- (S58_covA$T.M *100)/(S58_covA$T.M + S58_covA$T.U)
S58_covA["mePerc.A"] <- (S58_covA$A.M *100)/(S58_covA$A.M + S58_covA$A.U)
S58_covB <- read.table("58_S58_L001_R1_001_Bcomb.coverage.txt")
rownames(S58_covB) <- paste0(S58_covB$V1, "_", S58_covB$V2)
S58_covB <- S58_covB[,c(5,6,11,12)]
colnames(S58_covB) <- c("T.M","T.U","C.M","C.U")
S58_covB["mePerc.T"] <- (S58_covB$T.M *100)/(S58_covB$T.M + S58_covB$T.U)
S58_covB["mePerc.C"] <- (S58_covB$C.M *100)/(S58_covB$C.M + S58_covB$C.U)


#Biq Analyzer based: Allele specific
#S49
samtools fasta -1 S49_A_g1_R1.fa -2 S49_A_g1_R2.fa 49_S49_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S49_A_g1_R2.fa > S49_A_g1_R2_rc.fa

samtools fasta -1 S49_A_g2_R1.fa -2 S49_A_g2_R2.fa 49_S49_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.sort.bam
#Since R2 is reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S49_A_g2_R2.fa > S49_A_g2_R2_rc.fa

cat S49_A_g1_R1.fa S49_A_g1_R2_rc.fa > S49_A_g1_merge.fa
cat S49_A_g2_R1.fa S49_A_g2_R2_rc.fa > S49_A_g2_merge.fa

#S50
samtools fasta -1 S50_A_g1_R1.fa -2 S50_A_g1_R2.fa 50_S50_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S50_A_g1_R2.fa > S50_A_g1_R2_rc.fa

samtools fasta -1 S50_A_g2_R1.fa -2 S50_A_g2_R2.fa 50_S50_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.sort.bam
#Since R2 is reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S50_A_g2_R2.fa > S50_A_g2_R2_rc.fa

#S51
samtools fasta -1 S51_A_g1_R1.fa -2 S51_A_g1_R2.fa 51_S51_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S51_A_g1_R2.fa > S51_A_g1_R2_rc.fa

samtools fasta -1 S51_A_g2_R1.fa -2 S51_A_g2_R2.fa 51_S51_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.sort.bam
#Since R2 is reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S51_A_g2_R2.fa > S51_A_g2_R2_rc.fa



#S52
samtools fasta -1 S52_A_g1_R1.fa -2 S52_A_g1_R2.fa 52_S52_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S52_A_g1_R2.fa > S52_A_g1_R2_rc.fa

samtools fasta -1 S52_A_g2_R1.fa -2 S52_A_g2_R2.fa 52_S52_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.sort.bam
#Since R2 is reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S52_A_g2_R2.fa > S52_A_g2_R2_rc.fa


#S53
samtools fasta -1 S53_A_g1_R1.fa -2 S53_A_g1_R2.fa 53_S53_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome1.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S53_A_g1_R2.fa > S53_A_g1_R2_rc.fa

samtools fasta -1 S53_A_g2_R1.fa -2 S53_A_g2_R2.fa 53_S53_L001_R1_001_val_1_bismark.sortedByReadname_Acomb.genome2.sort.bam
#Since R2 is reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S53_A_g2_R2.fa > S53_A_g2_R2_rc.fa

#image convert
convert -resize 1024x1024 patternmap_Bisulfite.svg patternmap_Bisulfite.png
#S54
samtools fasta -1 S54_B_g1_R1.fa -2 S54_B_g1_R2.fa 54_S54_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S54_B_g1_R1.fa > S54_B_g1_R1_rc.fa

samtools fasta -1 S54_B_g2_R1.fa -2 S54_B_g2_R2.fa 54_S54_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S54_B_g2_R1.fa > S54_B_g2_R1_rc.fa

#S55
samtools fasta -1 S55_B_g1_R1.fa -2 S55_B_g1_R2.fa 55_S55_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S55_B_g1_R1.fa > S55_B_g1_R1_rc.fa

samtools fasta -1 S55_B_g2_R1.fa -2 S55_B_g2_R2.fa 55_S55_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S55_B_g2_R1.fa > S55_B_g2_R1_rc.fa

#S56
samtools fasta -1 S56_B_g1_R1.fa -2 S56_B_g1_R2.fa 56_S56_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S56_B_g1_R1.fa > S56_B_g1_R1_rc.fa

samtools fasta -1 S56_B_g2_R1.fa -2 S56_B_g2_R2.fa 56_S56_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S56_B_g2_R1.fa > S56_B_g2_R1_rc.fa


#S57
samtools fasta -1 S57_B_g1_R1.fa -2 S57_B_g1_R2.fa 57_S57_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S57_B_g1_R1.fa > S57_B_g1_R1_rc.fa

samtools fasta -1 S57_B_g2_R1.fa -2 S57_B_g2_R2.fa 57_S57_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S57_B_g2_R1.fa > S57_B_g2_R1_rc.fa


#S58
samtools fasta -1 S58_B_g1_R1.fa -2 S58_B_g1_R2.fa 58_S58_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome1.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S58_B_g1_R1.fa > S58_B_g1_R1_rc.fa

samtools fasta -1 S58_B_g2_R1.fa -2 S58_B_g2_R2.fa 58_S58_L001_R1_001_val_1_bismark.sortedByReadname_Bcomb.genome2.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S58_B_g2_R1.fa > S58_B_g2_R1_rc.fa
head S58_B_g2_R1_rc.fa -n 8000 > S58_B_g2_R1_rc_top.fa
head S58_B_g1_R1_rc.fa -n 8000 > S58_B_g1_R1_rc_top.fa
#Convert also reference Nnat
/home/ankitv/tools_av/seqtk/seqtk seq -r ref_Nnat.fa > ref_Nnat_rc.fa


#Biq Analyzer based: Bulk
#S49
samtools fasta -1 S49_R1.fa -2 S49_R2.fa 49_S49_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S49_R2.fa > S49_R2_rc.fa


#S50
samtools fasta -1 S50_R1.fa -2 S50_R2.fa 50_S50_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S50_R2.fa > S50_R2_rc.fa

#S51
samtools fasta -1 S51_R1.fa -2 S51_R2.fa 51_S51_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S51_R2.fa > S51_R2_rc.fa


#S52
samtools fasta -1 S52_R1.fa -2 S52_R2.fa 52_S52_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S52_R2.fa > S52_R2_rc.fa


#S53
samtools fasta -1 S53_R1.fa -2 S53_R2.fa 53_S53_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R2 is in reverse complimentary it need to be bring to forward (R2 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S53_R2.fa > S53_R2_rc.fa

#image convert
convert -resize 1024x1024 patternmap_Bisulfite.svg patternmap_Bisulfite.png
#S54
samtools fasta -1 S54_R1.fa -2 S54_R2.fa 54_S54_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S54_R1.fa > S54_R1_rc.fa
head S54_R1_rc.fa -n 8000 > S54_R1_rc_top.fa

#S55
samtools fasta -1 S55_R1.fa -2 S55_R2.fa 55_S55_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S55_R1.fa > S55_R1_rc.fa

#S56
samtools fasta -1 S56_R1.fa -2 S56_R2.fa 56_S56_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S56_R1.fa > S56_R1_rc.fa
head S56_R1_rc.fa -n 8000 > S56_R1_rc_top.fa


#S57
samtools fasta -1 S57_R1.fa -2 S57_R2.fa 57_S57_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S57_R1.fa > S57_R1_rc.fa


#S58
samtools fasta -1 S58_R1.fa -2 S58_R2.fa 58_S58_L001_R1_001_val_1_bismark_bt2_pe.sort.bam
#Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)
/home/ankitv/tools_av/seqtk/seqtk seq -r S58_R1.fa > S58_R1_rc.fa

head S58_R1_rc.fa -n 10000 > S58_R1_rc_top.fa
#Convert also reference Nnat
/home/ankitv/tools_av/seqtk/seqtk seq -r ref_Nnat.fa > ref_Nnat_rc.fa

#Import results
#S49_R2
sort -k1,1 -u S49_R2_pattern.tsv > S49_R2_patternu.tsv
S49_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S49_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S49_R2_pattern) <- paste0(S49_R2_pattern$ID,"_",S49_R2_pattern$Sample)
dim(S49_R2_pattern)
head(S49_R2_pattern)
S49_R2_pattern <- subset(S49_R2_pattern, select=-c(ID,Reference,Sample))
S49_R2_pattern <-  data.frame(S49_R2_pattern)
S49_R2_patternst <- stack((S49_R2_pattern))
countS49_R2_patternst <- count(S49_R2_patternst,"values")
countS49_R2_patternst <- data.frame(t(countS49_R2_patternst))
colnames(countS49_R2_patternst) <- c("U", "M", "A")
countS49_R2_patternst <- countS49_R2_patternst[-1,]
countS49_R2_patternst["sample"] <- "S49_R2"
countS49_R2_patternst <- data.frame(countS49_R2_patternst)
countS49_R2_patternst$U <- as.numeric(as.character(countS49_R2_patternst$U))
countS49_R2_patternst$M <- as.numeric(as.character(countS49_R2_patternst$M))
countS49_R2_patternst["methCpGsprop"] <- countS49_R2_patternst$M / (countS49_R2_patternst$U + countS49_R2_patternst$M)


#S50_R2
sort -k1,1 -u S50_R2_pattern.tsv > S50_R2_patternu.tsv
S50_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S50_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S50_R2_pattern) <- paste0(S50_R2_pattern$ID,"_",S50_R2_pattern$Sample)
dim(S50_R2_pattern)
head(S50_R2_pattern)
S50_R2_pattern <- subset(S50_R2_pattern, select=-c(ID,Reference,Sample))
S50_R2_pattern <-  data.frame(S50_R2_pattern)
S50_R2_patternst <- stack((S50_R2_pattern))
countS50_R2_patternst <- count(S50_R2_patternst,"values")
countS50_R2_patternst <- data.frame(t(countS50_R2_patternst))
colnames(countS50_R2_patternst) <- c("U", "M", "A")
countS50_R2_patternst <- countS50_R2_patternst[-1,]
countS50_R2_patternst["sample"] <- "S50_R2"
countS50_R2_patternst <- data.frame(countS50_R2_patternst)
countS50_R2_patternst$U <- as.numeric(as.character(countS50_R2_patternst$U))
countS50_R2_patternst$M <- as.numeric(as.character(countS50_R2_patternst$M))
countS50_R2_patternst["methCpGsprop"] <- countS50_R2_patternst$M / (countS50_R2_patternst$U + countS50_R2_patternst$M)

#S51_R2
sort -k1,1 -u S51_R2_pattern.tsv > S51_R2_patternu.tsv
S51_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S51_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S51_R2_pattern) <- paste0(S51_R2_pattern$ID,"_",S51_R2_pattern$Sample)
dim(S51_R2_pattern)
head(S51_R2_pattern)
S51_R2_pattern <- subset(S51_R2_pattern, select=-c(ID,Reference,Sample))
S51_R2_pattern <-  data.frame(S51_R2_pattern)
S51_R2_patternst <- stack((S51_R2_pattern))
countS51_R2_patternst <- count(S51_R2_patternst,"values")
countS51_R2_patternst <- data.frame(t(countS51_R2_patternst))
colnames(countS51_R2_patternst) <- c("U", "M", "A")
countS51_R2_patternst <- countS51_R2_patternst[-1,]
countS51_R2_patternst["sample"] <- "S51_R2"
countS51_R2_patternst <- data.frame(countS51_R2_patternst)
countS51_R2_patternst$U <- as.numeric(as.character(countS51_R2_patternst$U))
countS51_R2_patternst$M <- as.numeric(as.character(countS51_R2_patternst$M))
countS51_R2_patternst["methCpGsprop"] <- countS51_R2_patternst$M / (countS51_R2_patternst$U + countS51_R2_patternst$M)


#S52_R2
sort -k1,1 -u S52_R2_pattern.tsv > S52_R2_patternu.tsv
S52_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S52_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S52_R2_pattern) <- paste0(S52_R2_pattern$ID,"_",S52_R2_pattern$Sample)
dim(S52_R2_pattern)
head(S52_R2_pattern)
S52_R2_pattern <- subset(S52_R2_pattern, select=-c(ID,Reference,Sample))
S52_R2_pattern <-  data.frame(S52_R2_pattern)
S52_R2_patternst <- stack((S52_R2_pattern))
countS52_R2_patternst <- count(S52_R2_patternst,"values")
countS52_R2_patternst <- data.frame(t(countS52_R2_patternst))
colnames(countS52_R2_patternst) <- c("U", "M", "A")
countS52_R2_patternst <- countS52_R2_patternst[-1,]
countS52_R2_patternst["sample"] <- "S52_R2"
countS52_R2_patternst <- data.frame(countS52_R2_patternst)
countS52_R2_patternst$U <- as.numeric(as.character(countS52_R2_patternst$U))
countS52_R2_patternst$M <- as.numeric(as.character(countS52_R2_patternst$M))
countS52_R2_patternst["methCpGsprop"] <- countS52_R2_patternst$M / (countS52_R2_patternst$U + countS52_R2_patternst$M)

#S53_R2
sort -k1,1 -u S53_R2_pattern.tsv > S53_R2_patternu.tsv
S53_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S53_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S53_R2_pattern) <- paste0(S53_R2_pattern$ID,"_",S53_R2_pattern$Sample)
dim(S53_R2_pattern)
head(S53_R2_pattern)
S53_R2_pattern <- subset(S53_R2_pattern, select=-c(ID,Reference,Sample))
S53_R2_pattern <-  data.frame(S53_R2_pattern)
S53_R2_patternst <- stack((S53_R2_pattern))
countS53_R2_patternst <- count(S53_R2_patternst,"values")
countS53_R2_patternst <- data.frame(t(countS53_R2_patternst))
colnames(countS53_R2_patternst) <- c("U", "M", "A")
countS53_R2_patternst <- countS53_R2_patternst[-1,]
countS53_R2_patternst["sample"] <- "S53_R2"
countS53_R2_patternst <- data.frame(countS53_R2_patternst)
countS53_R2_patternst$U <- as.numeric(as.character(countS53_R2_patternst$U))
countS53_R2_patternst$M <- as.numeric(as.character(countS53_R2_patternst$M))
countS53_R2_patternst["methCpGsprop"] <- countS53_R2_patternst$M / (countS53_R2_patternst$U + countS53_R2_patternst$M)

#S49_genome1_R2
sort -k1,1 -u S49_genome1_R2_pattern.tsv > S49_genome1_R2_patternu.tsv
S49_genome1_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S49_genome1_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S49_genome1_R2_pattern) <- paste0(S49_genome1_R2_pattern$ID,"_",S49_genome1_R2_pattern$Sample)
dim(S49_genome1_R2_pattern)
head(S49_genome1_R2_pattern)
S49_genome1_R2_pattern <- subset(S49_genome1_R2_pattern, select=-c(ID,Reference,Sample))
S49_genome1_R2_pattern <-  data.frame(S49_genome1_R2_pattern)
S49_genome1_R2_patternst <- stack((S49_genome1_R2_pattern))
countS49_genome1_R2_patternst <- count(S49_genome1_R2_patternst,"values")
countS49_genome1_R2_patternst <- data.frame(t(countS49_genome1_R2_patternst))
colnames(countS49_genome1_R2_patternst) <- c("U", "M", "A")
countS49_genome1_R2_patternst <- countS49_genome1_R2_patternst[-1,]
countS49_genome1_R2_patternst["sample"] <- "S49_genome1_R2"
countS49_genome1_R2_patternst <- data.frame(countS49_genome1_R2_patternst)
countS49_genome1_R2_patternst$U <- as.numeric(as.character(countS49_genome1_R2_patternst$U))
countS49_genome1_R2_patternst$M <- as.numeric(as.character(countS49_genome1_R2_patternst$M))
countS49_genome1_R2_patternst["methCpGsprop"] <- countS49_genome1_R2_patternst$M / (countS49_genome1_R2_patternst$U + countS49_genome1_R2_patternst$M)


#S50_genome1_R2
sort -k1,1 -u S50_genome1_R2_pattern.tsv > S50_genome1_R2_patternu.tsv
S50_genome1_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S50_genome1_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S50_genome1_R2_pattern) <- paste0(S50_genome1_R2_pattern$ID,"_",S50_genome1_R2_pattern$Sample)
dim(S50_genome1_R2_pattern)
head(S50_genome1_R2_pattern)
S50_genome1_R2_pattern <- subset(S50_genome1_R2_pattern, select=-c(ID,Reference,Sample))
S50_genome1_R2_pattern <-  data.frame(S50_genome1_R2_pattern)
S50_genome1_R2_patternst <- stack((S50_genome1_R2_pattern))
countS50_genome1_R2_patternst <- count(S50_genome1_R2_patternst,"values")
countS50_genome1_R2_patternst <- data.frame(t(countS50_genome1_R2_patternst))
colnames(countS50_genome1_R2_patternst) <- c("U", "M", "A")
countS50_genome1_R2_patternst <- countS50_genome1_R2_patternst[-1,]
countS50_genome1_R2_patternst["sample"] <- "S50_genome1_R2"
countS50_genome1_R2_patternst <- data.frame(countS50_genome1_R2_patternst)
countS50_genome1_R2_patternst$U <- as.numeric(as.character(countS50_genome1_R2_patternst$U))
countS50_genome1_R2_patternst$M <- as.numeric(as.character(countS50_genome1_R2_patternst$M))
countS50_genome1_R2_patternst["methCpGsprop"] <- countS50_genome1_R2_patternst$M / (countS50_genome1_R2_patternst$U + countS50_genome1_R2_patternst$M)

#S51_genome1_R2
sort -k1,1 -u S51_genome1_R2_pattern.tsv > S51_genome1_R2_patternu.tsv
S51_genome1_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S51_genome1_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S51_genome1_R2_pattern) <- paste0(S51_genome1_R2_pattern$ID,"_",S51_genome1_R2_pattern$Sample)
dim(S51_genome1_R2_pattern)
head(S51_genome1_R2_pattern)
S51_genome1_R2_pattern <- subset(S51_genome1_R2_pattern, select=-c(ID,Reference,Sample))
S51_genome1_R2_pattern <-  data.frame(S51_genome1_R2_pattern)
S51_genome1_R2_patternst <- stack((S51_genome1_R2_pattern))
countS51_genome1_R2_patternst <- count(S51_genome1_R2_patternst,"values")
countS51_genome1_R2_patternst <- data.frame(t(countS51_genome1_R2_patternst))
colnames(countS51_genome1_R2_patternst) <- c("U", "M", "A")
countS51_genome1_R2_patternst <- countS51_genome1_R2_patternst[-1,]
countS51_genome1_R2_patternst["sample"] <- "S51_genome1_R2"
countS51_genome1_R2_patternst <- data.frame(countS51_genome1_R2_patternst)
countS51_genome1_R2_patternst$U <- as.numeric(as.character(countS51_genome1_R2_patternst$U))
countS51_genome1_R2_patternst$M <- as.numeric(as.character(countS51_genome1_R2_patternst$M))
countS51_genome1_R2_patternst["methCpGsprop"] <- countS51_genome1_R2_patternst$M / (countS51_genome1_R2_patternst$U + countS51_genome1_R2_patternst$M)


#S52_genome1_R2
sort -k1,1 -u S52_genome1_R2_pattern.tsv > S52_genome1_R2_patternu.tsv
S52_genome1_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S52_genome1_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S52_genome1_R2_pattern) <- paste0(S52_genome1_R2_pattern$ID,"_",S52_genome1_R2_pattern$Sample)
dim(S52_genome1_R2_pattern)
head(S52_genome1_R2_pattern)
S52_genome1_R2_pattern <- subset(S52_genome1_R2_pattern, select=-c(ID,Reference,Sample))
S52_genome1_R2_pattern <-  data.frame(S52_genome1_R2_pattern)
S52_genome1_R2_patternst <- stack((S52_genome1_R2_pattern))
countS52_genome1_R2_patternst <- count(S52_genome1_R2_patternst,"values")
countS52_genome1_R2_patternst <- data.frame(t(countS52_genome1_R2_patternst))
colnames(countS52_genome1_R2_patternst) <- c("U", "M", "A")
countS52_genome1_R2_patternst <- countS52_genome1_R2_patternst[-1,]
countS52_genome1_R2_patternst["sample"] <- "S52_genome1_R2"
countS52_genome1_R2_patternst <- data.frame(countS52_genome1_R2_patternst)
countS52_genome1_R2_patternst$U <- as.numeric(as.character(countS52_genome1_R2_patternst$U))
countS52_genome1_R2_patternst$M <- as.numeric(as.character(countS52_genome1_R2_patternst$M))
countS52_genome1_R2_patternst["methCpGsprop"] <- countS52_genome1_R2_patternst$M / (countS52_genome1_R2_patternst$U + countS52_genome1_R2_patternst$M)

#S53_genome1_R2
sort -k1,1 -u S53_genome1_R2_pattern.tsv > S53_genome1_R2_patternu.tsv
S53_genome1_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S53_genome1_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S53_genome1_R2_pattern) <- paste0(S53_genome1_R2_pattern$ID,"_",S53_genome1_R2_pattern$Sample)
dim(S53_genome1_R2_pattern)
head(S53_genome1_R2_pattern)
S53_genome1_R2_pattern <- subset(S53_genome1_R2_pattern, select=-c(ID,Reference,Sample))
S53_genome1_R2_pattern <-  data.frame(S53_genome1_R2_pattern)
S53_genome1_R2_patternst <- stack((S53_genome1_R2_pattern))
countS53_genome1_R2_patternst <- count(S53_genome1_R2_patternst,"values")
countS53_genome1_R2_patternst <- data.frame(t(countS53_genome1_R2_patternst))
colnames(countS53_genome1_R2_patternst) <- c("U", "M", "A")
countS53_genome1_R2_patternst <- countS53_genome1_R2_patternst[-1,]
countS53_genome1_R2_patternst["sample"] <- "S53_genome1_R2"
countS53_genome1_R2_patternst <- data.frame(countS53_genome1_R2_patternst)
countS53_genome1_R2_patternst$U <- as.numeric(as.character(countS53_genome1_R2_patternst$U))
countS53_genome1_R2_patternst$M <- as.numeric(as.character(countS53_genome1_R2_patternst$M))
countS53_genome1_R2_patternst["methCpGsprop"] <- countS53_genome1_R2_patternst$M / (countS53_genome1_R2_patternst$U + countS53_genome1_R2_patternst$M)


#S49_genome2_R2
sort -k1,1 -u S49_genome2_R2_pattern.tsv > S49_genome2_R2_patternu.tsv
S49_genome2_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S49_genome2_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S49_genome2_R2_pattern) <- paste0(S49_genome2_R2_pattern$ID,"_",S49_genome2_R2_pattern$Sample)
dim(S49_genome2_R2_pattern)
head(S49_genome2_R2_pattern)
S49_genome2_R2_pattern <- subset(S49_genome2_R2_pattern, select=-c(ID,Reference,Sample))
S49_genome2_R2_pattern <-  data.frame(S49_genome2_R2_pattern)
S49_genome2_R2_patternst <- stack((S49_genome2_R2_pattern))
countS49_genome2_R2_patternst <- count(S49_genome2_R2_patternst,"values")
countS49_genome2_R2_patternst <- data.frame(t(countS49_genome2_R2_patternst))
colnames(countS49_genome2_R2_patternst) <- c("U", "M", "A")
countS49_genome2_R2_patternst <- countS49_genome2_R2_patternst[-1,]
countS49_genome2_R2_patternst["sample"] <- "S49_genome2_R2"
countS49_genome2_R2_patternst <- data.frame(countS49_genome2_R2_patternst)
countS49_genome2_R2_patternst$U <- as.numeric(as.character(countS49_genome2_R2_patternst$U))
countS49_genome2_R2_patternst$M <- as.numeric(as.character(countS49_genome2_R2_patternst$M))
countS49_genome2_R2_patternst["methCpGsprop"] <- countS49_genome2_R2_patternst$M / (countS49_genome2_R2_patternst$U + countS49_genome2_R2_patternst$M)


#S50_genome2_R2
sort -k1,1 -u S50_genome2_R2_pattern.tsv > S50_genome2_R2_patternu.tsv
S50_genome2_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S50_genome2_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S50_genome2_R2_pattern) <- paste0(S50_genome2_R2_pattern$ID,"_",S50_genome2_R2_pattern$Sample)
dim(S50_genome2_R2_pattern)
head(S50_genome2_R2_pattern)
S50_genome2_R2_pattern <- subset(S50_genome2_R2_pattern, select=-c(ID,Reference,Sample))
S50_genome2_R2_pattern <-  data.frame(S50_genome2_R2_pattern)
S50_genome2_R2_patternst <- stack((S50_genome2_R2_pattern))
countS50_genome2_R2_patternst <- count(S50_genome2_R2_patternst,"values")
countS50_genome2_R2_patternst <- data.frame(t(countS50_genome2_R2_patternst))
colnames(countS50_genome2_R2_patternst) <- c("U", "M", "A")
countS50_genome2_R2_patternst <- countS50_genome2_R2_patternst[-1,]
countS50_genome2_R2_patternst["sample"] <- "S50_genome2_R2"
countS50_genome2_R2_patternst <- data.frame(countS50_genome2_R2_patternst)
countS50_genome2_R2_patternst$U <- as.numeric(as.character(countS50_genome2_R2_patternst$U))
countS50_genome2_R2_patternst$M <- as.numeric(as.character(countS50_genome2_R2_patternst$M))
countS50_genome2_R2_patternst["methCpGsprop"] <- countS50_genome2_R2_patternst$M / (countS50_genome2_R2_patternst$U + countS50_genome2_R2_patternst$M)

#S51_genome2_R2
sort -k1,1 -u S51_genome2_R2_pattern.tsv > S51_genome2_R2_patternu.tsv
S51_genome2_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S51_genome2_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S51_genome2_R2_pattern) <- paste0(S51_genome2_R2_pattern$ID,"_",S51_genome2_R2_pattern$Sample)
dim(S51_genome2_R2_pattern)
head(S51_genome2_R2_pattern)
S51_genome2_R2_pattern <- subset(S51_genome2_R2_pattern, select=-c(ID,Reference,Sample))
S51_genome2_R2_pattern <-  data.frame(S51_genome2_R2_pattern)
S51_genome2_R2_patternst <- stack((S51_genome2_R2_pattern))
countS51_genome2_R2_patternst <- count(S51_genome2_R2_patternst,"values")
countS51_genome2_R2_patternst <- data.frame(t(countS51_genome2_R2_patternst))
colnames(countS51_genome2_R2_patternst) <- c("U", "M", "A")
countS51_genome2_R2_patternst <- countS51_genome2_R2_patternst[-1,]
countS51_genome2_R2_patternst["sample"] <- "S51_genome2_R2"
countS51_genome2_R2_patternst <- data.frame(countS51_genome2_R2_patternst)
countS51_genome2_R2_patternst$U <- as.numeric(as.character(countS51_genome2_R2_patternst$U))
countS51_genome2_R2_patternst$M <- as.numeric(as.character(countS51_genome2_R2_patternst$M))
countS51_genome2_R2_patternst["methCpGsprop"] <- countS51_genome2_R2_patternst$M / (countS51_genome2_R2_patternst$U + countS51_genome2_R2_patternst$M)


#S52_genome2_R2
sort -k1,1 -u S52_genome2_R2_pattern.tsv > S52_genome2_R2_patternu.tsv
S52_genome2_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S52_genome2_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S52_genome2_R2_pattern) <- paste0(S52_genome2_R2_pattern$ID,"_",S52_genome2_R2_pattern$Sample)
dim(S52_genome2_R2_pattern)
head(S52_genome2_R2_pattern)
S52_genome2_R2_pattern <- subset(S52_genome2_R2_pattern, select=-c(ID,Reference,Sample))
S52_genome2_R2_pattern <-  data.frame(S52_genome2_R2_pattern)
S52_genome2_R2_patternst <- stack((S52_genome2_R2_pattern))
countS52_genome2_R2_patternst <- count(S52_genome2_R2_patternst,"values")
countS52_genome2_R2_patternst <- data.frame(t(countS52_genome2_R2_patternst))
colnames(countS52_genome2_R2_patternst) <- c("U", "M", "A")
countS52_genome2_R2_patternst <- countS52_genome2_R2_patternst[-1,]
countS52_genome2_R2_patternst["sample"] <- "S52_genome2_R2"
countS52_genome2_R2_patternst <- data.frame(countS52_genome2_R2_patternst)
countS52_genome2_R2_patternst$U <- as.numeric(as.character(countS52_genome2_R2_patternst$U))
countS52_genome2_R2_patternst$M <- as.numeric(as.character(countS52_genome2_R2_patternst$M))
countS52_genome2_R2_patternst["methCpGsprop"] <- countS52_genome2_R2_patternst$M / (countS52_genome2_R2_patternst$U + countS52_genome2_R2_patternst$M)

#S53_genome2_R2
sort -k1,1 -u S53_genome2_R2_pattern.tsv > S53_genome2_R2_patternu.tsv
S53_genome2_R2_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S53_genome2_R2_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S53_genome2_R2_pattern) <- paste0(S53_genome2_R2_pattern$ID,"_",S53_genome2_R2_pattern$Sample)
dim(S53_genome2_R2_pattern)
head(S53_genome2_R2_pattern)
S53_genome2_R2_pattern <- subset(S53_genome2_R2_pattern, select=-c(ID,Reference,Sample))
S53_genome2_R2_pattern <-  data.frame(S53_genome2_R2_pattern)
S53_genome2_R2_patternst <- stack((S53_genome2_R2_pattern))
countS53_genome2_R2_patternst <- count(S53_genome2_R2_patternst,"values")
countS53_genome2_R2_patternst <- data.frame(t(countS53_genome2_R2_patternst))
colnames(countS53_genome2_R2_patternst) <- c("U", "M", "A")
countS53_genome2_R2_patternst <- countS53_genome2_R2_patternst[-1,]
countS53_genome2_R2_patternst["sample"] <- "S53_genome2_R2"
countS53_genome2_R2_patternst <- data.frame(countS53_genome2_R2_patternst)
countS53_genome2_R2_patternst$U <- as.numeric(as.character(countS53_genome2_R2_patternst$U))
countS53_genome2_R2_patternst$M <- as.numeric(as.character(countS53_genome2_R2_patternst$M))
countS53_genome2_R2_patternst["methCpGsprop"] <- countS53_genome2_R2_patternst$M / (countS53_genome2_R2_patternst$U + countS53_genome2_R2_patternst$M)

count_S49_to_S53_R2_patternst <- rbind.data.frame(countS49_R2_patternst,
                                                  countS50_R2_patternst,
                                                  countS51_R2_patternst,
                                                  countS52_R2_patternst,
                                                  countS53_R2_patternst,
                                                  countS49_genome1_R2_patternst,
                                                  countS50_genome1_R2_patternst,
                                                  countS51_genome1_R2_patternst,
                                                  countS52_genome1_R2_patternst,
                                                  countS53_genome1_R2_patternst,
                                                  countS49_genome2_R2_patternst,
                                                  countS50_genome2_R2_patternst,
                                                  countS51_genome2_R2_patternst,
                                                  countS52_genome2_R2_patternst,
                                                  countS53_genome2_R2_patternst)
count_S49_to_S53_R2_patternst["row"] <- rep(c("S49","S50","S51","S52","S53"),3)
count_S49_to_S53_R2_patternst <- data.frame(count_S49_to_S53_R2_patternst)
ggplot(count_S49_to_S53_R2_patternst, aes(x=row, y=methCpGsprop,color=sample))+
  geom_point(aes(fill=row), position=position_dodge(width=0.2))+ ylim(0,1)+theme_classic()+
  scale_color_manual(values = c(rep(c("green","blue","black"),5)))

ggsave("count_S49_to_S53_R2_patternst.svg", width=10*1.25, height=6*1.25, units="cm", dpi=96)

sort -k1,1 -u S54_R1_pattern.tsv > S54_R1_patternu.tsv
sort -k1,1 -u S55_R1_pattern.tsv > S55_R1_patternu.tsv
sort -k1,1 -u S56_R1_pattern.tsv > S56_R1_patternu.tsv
sort -k1,1 -u S57_R1_pattern.tsv > S57_R1_patternu.tsv
sort -k1,1 -u S58_R1_pattern.tsv > S58_R1_patternu.tsv
sort -k1,1 -u S54_genome1_R1_pattern.tsv > S54_genome1_R1_patternu.tsv
sort -k1,1 -u S55_genome1_R1_pattern.tsv > S55_genome1_R1_patternu.tsv
sort -k1,1 -u S56_genome1_R1_pattern.tsv > S56_genome1_R1_patternu.tsv
sort -k1,1 -u S57_genome1_R1_pattern.tsv > S57_genome1_R1_patternu.tsv
sort -k1,1 -u S58_genome1_R1_pattern.tsv > S58_genome1_R1_patternu.tsv
sort -k1,1 -u S54_genome2_R1_pattern.tsv > S54_genome2_R1_patternu.tsv
sort -k1,1 -u S55_genome2_R1_pattern.tsv > S55_genome2_R1_patternu.tsv
sort -k1,1 -u S56_genome2_R1_pattern.tsv > S56_genome2_R1_patternu.tsv
sort -k1,1 -u S57_genome2_R1_pattern.tsv > S57_genome2_R1_patternu.tsv
sort -k1,1 -u S58_genome2_R1_pattern.tsv > S58_genome2_R1_patternu.tsv

#S54_R1
sort -k1,1 -u S54_R1_pattern.tsv > S54_R1_patternu.tsv
S54_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S54_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S54_R1_pattern) <- paste0(S54_R1_pattern$ID,"_",S54_R1_pattern$Sample)
dim(S54_R1_pattern)
head(S54_R1_pattern)
S54_R1_pattern <- subset(S54_R1_pattern, select=-c(ID,Reference,Sample))
S54_R1_pattern <-  data.frame(S54_R1_pattern)
S54_R1_patternst <- stack((S54_R1_pattern))
countS54_R1_patternst <- count(S54_R1_patternst,"values")
countS54_R1_patternst <- data.frame(t(countS54_R1_patternst))
colnames(countS54_R1_patternst) <- c("U", "M", "A")
countS54_R1_patternst <- countS54_R1_patternst[-1,]
countS54_R1_patternst["sample"] <- "S54_R1"
countS54_R1_patternst <- data.frame(countS54_R1_patternst)
countS54_R1_patternst$U <- as.numeric(as.character(countS54_R1_patternst$U))
countS54_R1_patternst$M <- as.numeric(as.character(countS54_R1_patternst$M))
countS54_R1_patternst["methCpGsprop"] <- countS54_R1_patternst$M / (countS54_R1_patternst$U + countS54_R1_patternst$M)


#S55_R1
sort -k1,1 -u S55_R1_pattern.tsv > S55_R1_patternu.tsv
S55_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S55_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S55_R1_pattern) <- paste0(S55_R1_pattern$ID,"_",S55_R1_pattern$Sample)
dim(S55_R1_pattern)
head(S55_R1_pattern)
S55_R1_pattern <- subset(S55_R1_pattern, select=-c(ID,Reference,Sample))
S55_R1_pattern <-  data.frame(S55_R1_pattern)
S55_R1_patternst <- stack((S55_R1_pattern))
countS55_R1_patternst <- count(S55_R1_patternst,"values")
countS55_R1_patternst <- data.frame(t(countS55_R1_patternst))
colnames(countS55_R1_patternst) <- c("U", "M", "A")
countS55_R1_patternst <- countS55_R1_patternst[-1,]
countS55_R1_patternst["sample"] <- "S55_R1"
countS55_R1_patternst <- data.frame(countS55_R1_patternst)
countS55_R1_patternst$U <- as.numeric(as.character(countS55_R1_patternst$U))
countS55_R1_patternst$M <- as.numeric(as.character(countS55_R1_patternst$M))
countS55_R1_patternst["methCpGsprop"] <- countS55_R1_patternst$M / (countS55_R1_patternst$U + countS55_R1_patternst$M)

#S56_R1
sort -k1,1 -u S56_R1_pattern.tsv > S56_R1_patternu.tsv
S56_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S56_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S56_R1_pattern) <- paste0(S56_R1_pattern$ID,"_",S56_R1_pattern$Sample)
dim(S56_R1_pattern)
head(S56_R1_pattern)
S56_R1_pattern <- subset(S56_R1_pattern, select=-c(ID,Reference,Sample))
S56_R1_pattern <-  data.frame(S56_R1_pattern)
S56_R1_patternst <- stack((S56_R1_pattern))
countS56_R1_patternst <- count(S56_R1_patternst,"values")
countS56_R1_patternst <- data.frame(t(countS56_R1_patternst))
colnames(countS56_R1_patternst) <- c("U", "M", "A")
countS56_R1_patternst <- countS56_R1_patternst[-1,]
countS56_R1_patternst["sample"] <- "S56_R1"
countS56_R1_patternst <- data.frame(countS56_R1_patternst)
countS56_R1_patternst$U <- as.numeric(as.character(countS56_R1_patternst$U))
countS56_R1_patternst$M <- as.numeric(as.character(countS56_R1_patternst$M))
countS56_R1_patternst["methCpGsprop"] <- countS56_R1_patternst$M / (countS56_R1_patternst$U + countS56_R1_patternst$M)


#S57_R1
sort -k1,1 -u S57_R1_pattern.tsv > S57_R1_patternu.tsv
S57_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S57_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S57_R1_pattern) <- paste0(S57_R1_pattern$ID,"_",S57_R1_pattern$Sample)
dim(S57_R1_pattern)
head(S57_R1_pattern)
S57_R1_pattern <- subset(S57_R1_pattern, select=-c(ID,Reference,Sample))
S57_R1_pattern <-  data.frame(S57_R1_pattern)
S57_R1_patternst <- stack((S57_R1_pattern))
countS57_R1_patternst <- count(S57_R1_patternst,"values")
countS57_R1_patternst <- data.frame(t(countS57_R1_patternst))
colnames(countS57_R1_patternst) <- c("U", "M", "A")
countS57_R1_patternst <- countS57_R1_patternst[-1,]
countS57_R1_patternst["sample"] <- "S57_R1"
countS57_R1_patternst <- data.frame(countS57_R1_patternst)
countS57_R1_patternst$U <- as.numeric(as.character(countS57_R1_patternst$U))
countS57_R1_patternst$M <- as.numeric(as.character(countS57_R1_patternst$M))
countS57_R1_patternst["methCpGsprop"] <- countS57_R1_patternst$M / (countS57_R1_patternst$U + countS57_R1_patternst$M)

#S58_R1
sort -k1,1 -u S58_R1_pattern.tsv > S58_R1_patternu.tsv
S58_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S58_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S58_R1_pattern) <- paste0(S58_R1_pattern$ID,"_",S58_R1_pattern$Sample)
dim(S58_R1_pattern)
head(S58_R1_pattern)
S58_R1_pattern <- subset(S58_R1_pattern, select=-c(ID,Reference,Sample))
S58_R1_pattern <-  data.frame(S58_R1_pattern)
S58_R1_patternst <- stack((S58_R1_pattern))
countS58_R1_patternst <- count(S58_R1_patternst,"values")
countS58_R1_patternst <- data.frame(t(countS58_R1_patternst))
colnames(countS58_R1_patternst) <- c("U", "M", "A")
countS58_R1_patternst <- countS58_R1_patternst[-1,]
countS58_R1_patternst["sample"] <- "S58_R1"
countS58_R1_patternst <- data.frame(countS58_R1_patternst)
countS58_R1_patternst$U <- as.numeric(as.character(countS58_R1_patternst$U))
countS58_R1_patternst$M <- as.numeric(as.character(countS58_R1_patternst$M))
countS58_R1_patternst["methCpGsprop"] <- countS58_R1_patternst$M / (countS58_R1_patternst$U + countS58_R1_patternst$M)


#S54_genome1_R1
sort -k1,1 -u S54_genome1_R1_pattern.tsv > S54_genome1_R1_patternu.tsv
S54_genome1_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S54_genome1_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S54_genome1_R1_pattern) <- paste0(S54_genome1_R1_pattern$ID,"_",S54_genome1_R1_pattern$Sample)
dim(S54_genome1_R1_pattern)
head(S54_genome1_R1_pattern)
S54_genome1_R1_pattern <- subset(S54_genome1_R1_pattern, select=-c(ID,Reference,Sample))
S54_genome1_R1_pattern <-  data.frame(S54_genome1_R1_pattern)
S54_genome1_R1_patternst <- stack((S54_genome1_R1_pattern))
countS54_genome1_R1_patternst <- count(S54_genome1_R1_patternst,"values")
countS54_genome1_R1_patternst <- data.frame(t(countS54_genome1_R1_patternst))
colnames(countS54_genome1_R1_patternst) <- c("U", "M", "A")
countS54_genome1_R1_patternst <- countS54_genome1_R1_patternst[-1,]
countS54_genome1_R1_patternst["sample"] <- "S54_genome1_R1"
countS54_genome1_R1_patternst <- data.frame(countS54_genome1_R1_patternst)
countS54_genome1_R1_patternst$U <- as.numeric(as.character(countS54_genome1_R1_patternst$U))
countS54_genome1_R1_patternst$M <- as.numeric(as.character(countS54_genome1_R1_patternst$M))
countS54_genome1_R1_patternst["methCpGsprop"] <- countS54_genome1_R1_patternst$M / (countS54_genome1_R1_patternst$U + countS54_genome1_R1_patternst$M)


#S55_genome1_R1
sort -k1,1 -u S55_genome1_R1_pattern.tsv > S55_genome1_R1_patternu.tsv
S55_genome1_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S55_genome1_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S55_genome1_R1_pattern) <- paste0(S55_genome1_R1_pattern$ID,"_",S55_genome1_R1_pattern$Sample)
dim(S55_genome1_R1_pattern)
head(S55_genome1_R1_pattern)
S55_genome1_R1_pattern <- subset(S55_genome1_R1_pattern, select=-c(ID,Reference,Sample))
S55_genome1_R1_pattern <-  data.frame(S55_genome1_R1_pattern)
S55_genome1_R1_patternst <- stack((S55_genome1_R1_pattern))
countS55_genome1_R1_patternst <- count(S55_genome1_R1_patternst,"values")
countS55_genome1_R1_patternst <- data.frame(t(countS55_genome1_R1_patternst))
colnames(countS55_genome1_R1_patternst) <- c("U", "M", "A")
countS55_genome1_R1_patternst <- countS55_genome1_R1_patternst[-1,]
countS55_genome1_R1_patternst["sample"] <- "S55_genome1_R1"
countS55_genome1_R1_patternst <- data.frame(countS55_genome1_R1_patternst)
countS55_genome1_R1_patternst$U <- as.numeric(as.character(countS55_genome1_R1_patternst$U))
countS55_genome1_R1_patternst$M <- as.numeric(as.character(countS55_genome1_R1_patternst$M))
countS55_genome1_R1_patternst["methCpGsprop"] <- countS55_genome1_R1_patternst$M / (countS55_genome1_R1_patternst$U + countS55_genome1_R1_patternst$M)

#S56_genome1_R1
sort -k1,1 -u S56_genome1_R1_pattern.tsv > S56_genome1_R1_patternu.tsv
S56_genome1_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S56_genome1_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S56_genome1_R1_pattern) <- paste0(S56_genome1_R1_pattern$ID,"_",S56_genome1_R1_pattern$Sample)
dim(S56_genome1_R1_pattern)
head(S56_genome1_R1_pattern)
S56_genome1_R1_pattern <- subset(S56_genome1_R1_pattern, select=-c(ID,Reference,Sample))
S56_genome1_R1_pattern <-  data.frame(S56_genome1_R1_pattern)
S56_genome1_R1_patternst <- stack((S56_genome1_R1_pattern))
countS56_genome1_R1_patternst <- count(S56_genome1_R1_patternst,"values")
countS56_genome1_R1_patternst <- data.frame(t(countS56_genome1_R1_patternst))
colnames(countS56_genome1_R1_patternst) <- c("U", "M", "A")
countS56_genome1_R1_patternst <- countS56_genome1_R1_patternst[-1,]
countS56_genome1_R1_patternst["sample"] <- "S56_genome1_R1"
countS56_genome1_R1_patternst <- data.frame(countS56_genome1_R1_patternst)
countS56_genome1_R1_patternst$U <- as.numeric(as.character(countS56_genome1_R1_patternst$U))
countS56_genome1_R1_patternst$M <- as.numeric(as.character(countS56_genome1_R1_patternst$M))
countS56_genome1_R1_patternst["methCpGsprop"] <- countS56_genome1_R1_patternst$M / (countS56_genome1_R1_patternst$U + countS56_genome1_R1_patternst$M)


#S57_genome1_R1
sort -k1,1 -u S57_genome1_R1_pattern.tsv > S57_genome1_R1_patternu.tsv
S57_genome1_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S57_genome1_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S57_genome1_R1_pattern) <- paste0(S57_genome1_R1_pattern$ID,"_",S57_genome1_R1_pattern$Sample)
dim(S57_genome1_R1_pattern)
head(S57_genome1_R1_pattern)
S57_genome1_R1_pattern <- subset(S57_genome1_R1_pattern, select=-c(ID,Reference,Sample))
S57_genome1_R1_pattern <-  data.frame(S57_genome1_R1_pattern)
S57_genome1_R1_patternst <- stack((S57_genome1_R1_pattern))
countS57_genome1_R1_patternst <- count(S57_genome1_R1_patternst,"values")
countS57_genome1_R1_patternst <- data.frame(t(countS57_genome1_R1_patternst))
colnames(countS57_genome1_R1_patternst) <- c("U", "M", "A")
countS57_genome1_R1_patternst <- countS57_genome1_R1_patternst[-1,]
countS57_genome1_R1_patternst["sample"] <- "S57_genome1_R1"
countS57_genome1_R1_patternst <- data.frame(countS57_genome1_R1_patternst)
countS57_genome1_R1_patternst$U <- as.numeric(as.character(countS57_genome1_R1_patternst$U))
countS57_genome1_R1_patternst$M <- as.numeric(as.character(countS57_genome1_R1_patternst$M))
countS57_genome1_R1_patternst["methCpGsprop"] <- countS57_genome1_R1_patternst$M / (countS57_genome1_R1_patternst$U + countS57_genome1_R1_patternst$M)

#S58_genome1_R1
sort -k1,1 -u S58_genome1_R1_pattern.tsv > S58_genome1_R1_patternu.tsv
S58_genome1_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S58_genome1_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S58_genome1_R1_pattern) <- paste0(S58_genome1_R1_pattern$ID,"_",S58_genome1_R1_pattern$Sample)
dim(S58_genome1_R1_pattern)
head(S58_genome1_R1_pattern)
S58_genome1_R1_pattern <- subset(S58_genome1_R1_pattern, select=-c(ID,Reference,Sample))
S58_genome1_R1_pattern <-  data.frame(S58_genome1_R1_pattern)
S58_genome1_R1_patternst <- stack((S58_genome1_R1_pattern))
countS58_genome1_R1_patternst <- count(S58_genome1_R1_patternst,"values")
countS58_genome1_R1_patternst <- data.frame(t(countS58_genome1_R1_patternst))
colnames(countS58_genome1_R1_patternst) <- c("U", "M", "A")
countS58_genome1_R1_patternst <- countS58_genome1_R1_patternst[-1,]
countS58_genome1_R1_patternst["sample"] <- "S58_genome1_R1"
countS58_genome1_R1_patternst <- data.frame(countS58_genome1_R1_patternst)
countS58_genome1_R1_patternst$U <- as.numeric(as.character(countS58_genome1_R1_patternst$U))
countS58_genome1_R1_patternst$M <- as.numeric(as.character(countS58_genome1_R1_patternst$M))
countS58_genome1_R1_patternst["methCpGsprop"] <- countS58_genome1_R1_patternst$M / (countS58_genome1_R1_patternst$U + countS58_genome1_R1_patternst$M)


#S54_genome2_R1
sort -k1,1 -u S54_genome2_R1_pattern.tsv > S54_genome2_R1_patternu.tsv
S54_genome2_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S54_genome2_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S54_genome2_R1_pattern) <- paste0(S54_genome2_R1_pattern$ID,"_",S54_genome2_R1_pattern$Sample)
dim(S54_genome2_R1_pattern)
head(S54_genome2_R1_pattern)
S54_genome2_R1_pattern <- subset(S54_genome2_R1_pattern, select=-c(ID,Reference,Sample))
S54_genome2_R1_pattern <-  data.frame(S54_genome2_R1_pattern)
S54_genome2_R1_patternst <- stack((S54_genome2_R1_pattern))
countS54_genome2_R1_patternst <- count(S54_genome2_R1_patternst,"values")
countS54_genome2_R1_patternst <- data.frame(t(countS54_genome2_R1_patternst))
colnames(countS54_genome2_R1_patternst) <- c("U", "M", "A")
countS54_genome2_R1_patternst <- countS54_genome2_R1_patternst[-1,]
countS54_genome2_R1_patternst["sample"] <- "S54_genome2_R1"
countS54_genome2_R1_patternst <- data.frame(countS54_genome2_R1_patternst)
countS54_genome2_R1_patternst$U <- as.numeric(as.character(countS54_genome2_R1_patternst$U))
countS54_genome2_R1_patternst$M <- as.numeric(as.character(countS54_genome2_R1_patternst$M))
countS54_genome2_R1_patternst["methCpGsprop"] <- countS54_genome2_R1_patternst$M / (countS54_genome2_R1_patternst$U + countS54_genome2_R1_patternst$M)


#S55_genome2_R1
sort -k1,1 -u S55_genome2_R1_pattern.tsv > S55_genome2_R1_patternu.tsv
S55_genome2_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S55_genome2_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S55_genome2_R1_pattern) <- paste0(S55_genome2_R1_pattern$ID,"_",S55_genome2_R1_pattern$Sample)
dim(S55_genome2_R1_pattern)
head(S55_genome2_R1_pattern)
S55_genome2_R1_pattern <- subset(S55_genome2_R1_pattern, select=-c(ID,Reference,Sample))
S55_genome2_R1_pattern <-  data.frame(S55_genome2_R1_pattern)
S55_genome2_R1_patternst <- stack((S55_genome2_R1_pattern))
countS55_genome2_R1_patternst <- count(S55_genome2_R1_patternst,"values")
countS55_genome2_R1_patternst <- data.frame(t(countS55_genome2_R1_patternst))
colnames(countS55_genome2_R1_patternst) <- c("U", "M", "A")
countS55_genome2_R1_patternst <- countS55_genome2_R1_patternst[-1,]
countS55_genome2_R1_patternst["sample"] <- "S55_genome2_R1"
countS55_genome2_R1_patternst <- data.frame(countS55_genome2_R1_patternst)
countS55_genome2_R1_patternst$U <- as.numeric(as.character(countS55_genome2_R1_patternst$U))
countS55_genome2_R1_patternst$M <- as.numeric(as.character(countS55_genome2_R1_patternst$M))
countS55_genome2_R1_patternst["methCpGsprop"] <- countS55_genome2_R1_patternst$M / (countS55_genome2_R1_patternst$U + countS55_genome2_R1_patternst$M)

#S56_genome2_R1
sort -k1,1 -u S56_genome2_R1_pattern.tsv > S56_genome2_R1_patternu.tsv
S56_genome2_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S56_genome2_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S56_genome2_R1_pattern) <- paste0(S56_genome2_R1_pattern$ID,"_",S56_genome2_R1_pattern$Sample)
dim(S56_genome2_R1_pattern)
head(S56_genome2_R1_pattern)
S56_genome2_R1_pattern <- subset(S56_genome2_R1_pattern, select=-c(ID,Reference,Sample))
S56_genome2_R1_pattern <-  data.frame(S56_genome2_R1_pattern)
S56_genome2_R1_patternst <- stack((S56_genome2_R1_pattern))
countS56_genome2_R1_patternst <- count(S56_genome2_R1_patternst,"values")
countS56_genome2_R1_patternst <- data.frame(t(countS56_genome2_R1_patternst))
colnames(countS56_genome2_R1_patternst) <- c("U", "M", "A")
countS56_genome2_R1_patternst <- countS56_genome2_R1_patternst[-1,]
countS56_genome2_R1_patternst["sample"] <- "S56_genome2_R1"
countS56_genome2_R1_patternst <- data.frame(countS56_genome2_R1_patternst)
countS56_genome2_R1_patternst$U <- as.numeric(as.character(countS56_genome2_R1_patternst$U))
countS56_genome2_R1_patternst$M <- as.numeric(as.character(countS56_genome2_R1_patternst$M))
countS56_genome2_R1_patternst["methCpGsprop"] <- countS56_genome2_R1_patternst$M / (countS56_genome2_R1_patternst$U + countS56_genome2_R1_patternst$M)


#S57_genome2_R1
sort -k1,1 -u S57_genome2_R1_pattern.tsv > S57_genome2_R1_patternu.tsv
S57_genome2_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S57_genome2_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S57_genome2_R1_pattern) <- paste0(S57_genome2_R1_pattern$ID,"_",S57_genome2_R1_pattern$Sample)
dim(S57_genome2_R1_pattern)
head(S57_genome2_R1_pattern)
S57_genome2_R1_pattern <- subset(S57_genome2_R1_pattern, select=-c(ID,Reference,Sample))
S57_genome2_R1_pattern <-  data.frame(S57_genome2_R1_pattern)
S57_genome2_R1_patternst <- stack((S57_genome2_R1_pattern))
countS57_genome2_R1_patternst <- count(S57_genome2_R1_patternst,"values")
countS57_genome2_R1_patternst <- data.frame(t(countS57_genome2_R1_patternst))
colnames(countS57_genome2_R1_patternst) <- c("U", "M", "A")
countS57_genome2_R1_patternst <- countS57_genome2_R1_patternst[-1,]
countS57_genome2_R1_patternst["sample"] <- "S57_genome2_R1"
countS57_genome2_R1_patternst <- data.frame(countS57_genome2_R1_patternst)
countS57_genome2_R1_patternst$U <- as.numeric(as.character(countS57_genome2_R1_patternst$U))
countS57_genome2_R1_patternst$M <- as.numeric(as.character(countS57_genome2_R1_patternst$M))
countS57_genome2_R1_patternst["methCpGsprop"] <- countS57_genome2_R1_patternst$M / (countS57_genome2_R1_patternst$U + countS57_genome2_R1_patternst$M)

#S58_genome2_R1
sort -k1,1 -u S58_genome2_R1_pattern.tsv > S58_genome2_R1_patternu.tsv
S58_genome2_R1_pattern <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/directalign/S58_genome2_R1_patternu.tsv", sep = "\t", header = T,stringsAsFactors= FALSE)
#rownames(S58_genome2_R1_pattern) <- paste0(S58_genome2_R1_pattern$ID,"_",S58_genome2_R1_pattern$Sample)
dim(S58_genome2_R1_pattern)
head(S58_genome2_R1_pattern)
S58_genome2_R1_pattern <- subset(S58_genome2_R1_pattern, select=-c(ID,Reference,Sample))
S58_genome2_R1_pattern <-  data.frame(S58_genome2_R1_pattern)
S58_genome2_R1_patternst <- stack((S58_genome2_R1_pattern))
countS58_genome2_R1_patternst <- count(S58_genome2_R1_patternst,"values")
countS58_genome2_R1_patternst <- data.frame(t(countS58_genome2_R1_patternst))
colnames(countS58_genome2_R1_patternst) <- c("U", "M", "A")
countS58_genome2_R1_patternst <- countS58_genome2_R1_patternst[-1,]
countS58_genome2_R1_patternst["sample"] <- "S58_genome2_R1"
countS58_genome2_R1_patternst <- data.frame(countS58_genome2_R1_patternst)
countS58_genome2_R1_patternst$U <- as.numeric(as.character(countS58_genome2_R1_patternst$U))
countS58_genome2_R1_patternst$M <- as.numeric(as.character(countS58_genome2_R1_patternst$M))
countS58_genome2_R1_patternst["methCpGsprop"] <- countS58_genome2_R1_patternst$M / (countS58_genome2_R1_patternst$U + countS58_genome2_R1_patternst$M)


count_S54_to_S58_R1_patternst <- rbind.data.frame(countS54_R1_patternst,
                                                  countS55_R1_patternst,
                                                  countS56_R1_patternst,
                                                  countS57_R1_patternst,
                                                  countS58_R1_patternst,
                                                  countS54_genome1_R1_patternst,
                                                  countS55_genome1_R1_patternst,
                                                  countS56_genome1_R1_patternst,
                                                  countS57_genome1_R1_patternst,
                                                  countS58_genome1_R1_patternst,
                                                  countS54_genome2_R1_patternst,
                                                  countS55_genome2_R1_patternst,
                                                  countS56_genome2_R1_patternst,
                                                  countS57_genome2_R1_patternst,
                                                  countS58_genome2_R1_patternst)
count_S54_to_S58_R1_patternst["row"] <- rep(c("S54","S55","S56","S57","S58"),3)
count_S54_to_S58_R1_patternst <- data.frame(count_S54_to_S58_R1_patternst)
ggplot(count_S54_to_S58_R1_patternst, aes(x=row, y=methCpGsprop,color=sample))+
  geom_point(aes(fill=row), position=position_dodge(width=0.2))+ ylim(0,1)+theme_classic()+
  scale_color_manual(values = c(rep(c("red","blue","black"),5)))

ggsave("count_S54_to_S58_R1_patternst.svg", width=10*1.25, height=6*1.25, units="cm", dpi=96)


#Get specific reads from Fastq
/home/ankitv/tools_av/bbmap/bbduk.sh in=49_S49_L001_R2_001.fastq out=unmatch_49_S49_L001_R2_001.fastq outm=match_49_S49_L001_R2_001.fastq literal=AGTTATTTTTATAGTGGAGAGAGATGGCGTTTAAGTGCAAATTTGTTAGTAGTTTTTT k=58




#New analysis amplicon Dec 2021
#Python script for bsseq_paired_nnatset.py
#Convert to fasta bsseq_paired_nnatset_tofasta.py (later conversion step was added to bsseq_paired_nnatset.py also to make a complete pipeline from one script)
#Prepare BiQ Analyzer sheet
#Follow these steps:

#Use this way since Nnat was amplified on reverse strand I used minus strand as forward strand
#Prep for set1
#ls -1 *.gz | awk -F'_' '{print $2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/"$1"_R1.fa"}' | tail -n 5 > minus_SseriesBiQsheetR1.tsv
#ls -1 *.gz | awk -F'_' '{print $2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/"$1"_R2.fa"}' | tail -n 5 > minus_SseriesBiQsheetR2.tsv

#Prep for set2
#minus, since Nnat is amplified in reverse strand and also BiQ manual suggest to use the strand which was originally amplified by PCR
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R1.fa"}' > minus_BiQsheetreplicate1.tsv
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R2.fa"}' > minus_BiQsheetreplicate2.tsv

#Note:***************************** VERY IMPORTANT: ADD HEADER to EACH TSV Files otherwise it will not take the first file ****************************

#Get names ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1" "$1"_L001_R1_001"" "$1"_L001_R2_001"}'
#Import BiQ Analyzer Samplesheet to BiQ HiMOD
#Run BiQ Analyzer (Double click icon --> Create new project --> Click directory where you want result --> Load from data structure table (.tsv) --> press |> run)
#Export -> All results in one TSV -> table_results.tsv
#Make another file which contain information about group sample_color.txt eg. N15_S7	Downs	N15_S7%Downs
#Run Python on reference fasta to get CG positions, determine_CpG_position_reverse.py (reverse because Nnat was reverse)

#------------------------minus--------------------------#
#Read 1 is  in  reverse as checked in IGV and fasta and BiQ alignmnet both for read and ref, 
#minus ref will be ok to use, assignment to CG position in reverse coordinataes will be Ok

#minus_R1
#Use python to compare results
python tabstsv_split.py > set2_meth_pattern_bypython_R1.txt
#Remove first column
grep "Sample" set2_meth_pattern_bypython_R1.txt -v > set2_meth_pattern_bypython_R1_re.txt
#---
set2_meth_pattern_bypython_minus_R1 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_Nseries/set2_meth_pattern_bypython_R1_re.txt", header = F, stringsAsFactors = F)
head(set2_meth_pattern_bypython_minus_R1)
rownames(set2_meth_pattern_bypython_minus_R1) <- paste0(set2_meth_pattern_bypython_minus_R1$V1, "%",
                                                        set2_meth_pattern_bypython_minus_R1$V2)
set2_meth_pattern_bypython_minus_R1 <- set2_meth_pattern_bypython_minus_R1[,-1]
head(set2_meth_pattern_bypython_minus_R1)
dim(set2_meth_pattern_bypython_minus_R1)
#Check of reads covered 
#Total reads after python based CpG extraction wc -l set2_meth_pattern_bypython_R1_re.txt = 549173 
#Total reads in N series R1 (divide by 2 as fasta first lines are >header) wc -l /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/*_R1.fa = 1098586/2 = 549173
#So all the reads for each sample are covered
set2_meth_pattern_bypython_minus_R1[set2_meth_pattern_bypython_minus_R1 == "x"] <- NA

#Remove NAs
set2_tab_ref_subNnat_minus_R1_Nseries_py <- na.omit(set2_meth_pattern_bypython_minus_R1)
head(set2_tab_ref_subNnat_minus_R1_Nseries_py)
dim(set2_tab_ref_subNnat_minus_R1_Nseries_py) #Same number as R step
summary(set2_tab_ref_subNnat_minus_R1_Nseries_py)
set2_tab_ref_subNnat_minus_R1_Nseries_py <- set2_tab_ref_subNnat_minus_R1_Nseries_py[,c(1,3:14)]
colnames(set2_tab_ref_subNnat_minus_R1_Nseries_py) <- c("Sample", "CG1r1","CG2r1","CG3r1","CG4r1","CG5r1","CG6r1","CG7r1","CG8r1","CG9r1","CG10r1","CG11r1","CG12r1")
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG1r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG1r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG2r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG2r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG3r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG3r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG4r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG4r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG5r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG5r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG6r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG6r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG7r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG7r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG8r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG8r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG9r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG9r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG10r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG10r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG11r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG11r1)
set2_tab_ref_subNnat_minus_R1_Nseries_py$CG12r1 <- as.numeric(set2_tab_ref_subNnat_minus_R1_Nseries_py$CG12r1)
head(set2_tab_ref_subNnat_minus_R1_Nseries_py)
dim(set2_tab_ref_subNnat_minus_R1_Nseries_py)
set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate = aggregate(set2_tab_ref_subNnat_minus_R1_Nseries_py[,c(2:13)],by=list(set2_tab_ref_subNnat_minus_R1_Nseries_py$Sample), mean)
head(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate)

#check indiv
#Test mean manually the aggregation process

for (i in unique(set2_tab_ref_subNnat_minus_R1_Nseries_py$Sample)){
   temporayset <- set2_tab_ref_subNnat_minus_R1_Nseries_py[which(set2_tab_ref_subNnat_minus_R1_Nseries_py$Sample == i),]
   print(i)
   print(colMeans(temporayset[,2:13]))
}
#copy the rstudio outport, adjust by \n  \t and other if required
#Import
test_R1_mean_Nseries_manual <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_Nseries/test_R1_mean_Nseries_manual.txt", header = F, stringsAsFactors = F)
colnames(test_R1_mean_Nseries_manual) <- colnames(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate)
rownames(test_R1_mean_Nseries_manual) <- test_R1_mean_Nseries_manual$Sample
#Check
head(test_R1_mean_Nseries_manual[,2:13])
head(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate[,2:13])
tail(test_R1_mean_Nseries_manual[,2:13])
tail(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate[,2:13])
colMeans(test_R1_mean_Nseries_manual[,2:13])
colMeans(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate[,2:13])
#So as also checked manually aggregate function calculates the mean methylation for each samples correctly
rownames(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate) <- set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate$Group.1
pheatmap::pheatmap(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate[,2:13])
colnames(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate) <- c("Sample", "CG1r1","CG2r1","CG3r1","CG4r1","CG5r1","CG6r1","CG7r1","CG8r1","CG9r1","CG10r1","CG11r1","CG12r1")
newsample_R1_color <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_Nseries/sample_color.txt")
colnames(newsample_R1_color) <- c("Sample","Condition","Combine")
set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col <- merge(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate, newsample_R1_color, by="Sample")
rownames(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col) <- set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col$Combine
set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col <- set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col[,c(2:13)]
head(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col)
#Remove undetermined
set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col <- set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col[c(-48),]
sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col <- stack(as.matrix(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col))
head(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col)
write.table(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col,"/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_Nseries/set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col.txt", sep = "\t", quote = F, append = F)
sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col <- data.frame(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col)
sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col <- splitstackshape::cSplit(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col,"row","%")
sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col <- data.frame(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col)
#python3 determine_CpG_position_reverse.py enter-- /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_nospace.fa
ref_Nnat_CGsitesR1 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_Nseries/CGSites.txt")
colnames(ref_Nnat_CGsitesR1) <- c("col", "sites")
sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_cg <- merge(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col, ref_Nnat_CGsitesR1, by="col")
dim(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_cg) #12 * 47 = 564

ggplot(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_cg, aes(x=as.numeric(sites), y=value, group=row_1))+
  geom_line(aes(color=row_1), size = 1)+
  geom_vline(xintercept = unique(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_cg$sites), colour = "lightblue", linetype="solid")+
  #geom_point(aes(color=row_1), shape=20)+
  theme_classic()+ ylim(0,1)+
  scale_color_manual(values=c(rep("red",29),rep("#00CC00",18),"blue"))
ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_Nseries/sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_cg.svg", width=25*1.25, height=9*1.25, units="cm", dpi=96)

#Read 2 is also in  reverse as checked in IGV and fasta and BiQ alignmnet both for read and ref, 
#minus ref will be ok to use, assignment to CG position in reverse coordinataes will be Ok
#minus_R2
#Use python to compare results
#read 2 covers 13 CpGs and first 7 need to be removed because it overlap last 7 CpGs of read 1 see image here /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/image_comparison_r1_r2_removal.jpg
python tabstsv_split_R2.py > set2_meth_pattern_bypython_R2.txt
#Remove first column
grep "Sample" set2_meth_pattern_bypython_R2.txt -v > set2_meth_pattern_bypython_R2_re.txt

set2_meth_pattern_bypython_minus_R2 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2_Nseries/set2_meth_pattern_bypython_R2_re.txt", header = F, stringsAsFactors = F)
head(set2_meth_pattern_bypython_minus_R2)
dim(set2_meth_pattern_bypython_minus_R2)
rownames(set2_meth_pattern_bypython_minus_R2) <- paste0(set2_meth_pattern_bypython_minus_R2$V1, "%",
                                                        set2_meth_pattern_bypython_minus_R2$V2)
set2_meth_pattern_bypython_minus_R2 <- set2_meth_pattern_bypython_minus_R2[,-1]
head(set2_meth_pattern_bypython_minus_R2)
set2_meth_pattern_bypython_minus_R2[set2_meth_pattern_bypython_minus_R2 == "x"] <- NA

#Remove NAs
set2_tab_ref_subNnat_minus_R2_Nseries_py <- na.omit(set2_meth_pattern_bypython_minus_R2)
head(set2_tab_ref_subNnat_minus_R2_Nseries_py)
dim(set2_tab_ref_subNnat_minus_R2_Nseries_py) #Same number as R step
summary(set2_tab_ref_subNnat_minus_R2_Nseries_py)
set2_tab_ref_subNnat_minus_R2_Nseries_py <- set2_tab_ref_subNnat_minus_R2_Nseries_py[,c(1,3:15)]
colnames(set2_tab_ref_subNnat_minus_R2_Nseries_py) <- c("Sample", "CG1r2","CG2r2","CG3r2","CG4r2","CG5r2","CG6r2","CG7r2","CG8r2","CG9r2","CG10r2","CG11r2","CG12r2","CG13r2")
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG1r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG1r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG2r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG2r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG3r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG3r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG4r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG4r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG5r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG5r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG6r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG6r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG7r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG7r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG8r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG8r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG9r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG9r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG10r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG10r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG11r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG11r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG12r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG12r2)
set2_tab_ref_subNnat_minus_R2_Nseries_py$CG13r2 <- as.numeric(set2_tab_ref_subNnat_minus_R2_Nseries_py$CG13r2)

dim(set2_tab_ref_subNnat_minus_R2_Nseries_py)
set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate = aggregate(set2_tab_ref_subNnat_minus_R2_Nseries_py[,c(2:14)],by=list(set2_tab_ref_subNnat_minus_R2_Nseries_py$Sample), mean)
head(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate)
rownames(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate) <- set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate$Group.1
pheatmap::pheatmap(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate[,2:14])
colnames(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate) <- c("Sample","CG1r2","CG2r2","CG3r2","CG4r2","CG5r2","CG6r2","CG7r2","CG8r2","CG9r2","CG10r2","CG11r2","CG12r2","CG13r2")
newsample_R2_color <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2_Nseries/sample_color.txt")
colnames(newsample_R2_color) <- c("Sample","Condition","Combine")
set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col <- merge(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate, newsample_R2_color, by="Sample")
rownames(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col) <- set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col$Combine
set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col <- set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col[,c(2:14)]
head(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col)
#Remove undetermined
set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col <- set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col[c(-48),]
sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col <- stack(as.matrix(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col))
head(sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col)
write.table(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col,"/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2_Nseries/set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col.txt", sep = "\t", quote = F, append = F)
sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col <- data.frame(sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col)
sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col <- splitstackshape::cSplit(sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col,"row","%")
sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col <- data.frame(sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col)
#python3 determine_CpG_position_reverse.py  enter-- /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2_nospace.fa
#tail -6 CGSites.txt > CGSites_filtered.txt
ref_Nnat_CGsitesR2 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2_Nseries/CGSites_filtered.txt")
colnames(ref_Nnat_CGsitesR2) <- c("col", "sites")
sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_cg <- merge(sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col, ref_Nnat_CGsitesR2, by="col")
dim(sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_cg)
ggplot(sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_cg, aes(x=as.numeric(sites), y=value, group=row_1))+
  geom_line(aes(color=row_2))+
  geom_vline(xintercept = unique(sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_cg$sites), colour = "lightblue", linetype="solid")+
  geom_point(aes(color=row_2), shape=20)+
  theme_classic()+ ylim(0,1)+
  scale_color_manual(values=c("#00CC00","red","blue"))
ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2_Nseries/sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_cg_merge.svg", width=25*1.25, height=9*1.25, units="cm", dpi=96)

dim(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_cg)
dim(sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_cg)
sset2_tab_ref_subNnat_minusR1nR2_Nseries_py_aggregate_cg <- rbind.data.frame(sset2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_cg,sset2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_cg)
dim(sset2_tab_ref_subNnat_minusR1nR2_Nseries_py_aggregate_cg)
head(sset2_tab_ref_subNnat_minusR1nR2_Nseries_py_aggregate_cg)
ggplot(sset2_tab_ref_subNnat_minusR1nR2_Nseries_py_aggregate_cg, aes(x=as.numeric(sites), y=value, group=row_1))+
  geom_line(aes(color=row_1), size = 1)+
  geom_vline(xintercept = unique(sset2_tab_ref_subNnat_minusR1nR2_Nseries_py_aggregate_cg$sites), colour = "lightblue", linetype="solid")+
  #geom_point(aes(color=row_1), shape=20)+
  theme_classic()+ ylim(0,1)+
  scale_color_manual(values=c(rep("red",29),rep("#00CC00",18),"blue"))
ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/sset2_tab_ref_subNnat_minusR1nR2_Nseries_py_aggregate_cg.svg", width=25*1.25, height=9*1.25, units="cm", dpi=96)

#Do the same for older data of five samples:remove S55, S56 and S58 as they are repeated, use minus reference
#S series
#Read 1 is  in  reverse as checked in IGV and fasta and BiQ alignmnet both for read and ref, 
#minus ref will be ok to use, assignment to CG position in reverse coordinataes will be Ok

#minus_R1
#Use python to compare results
cd /media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R1_Sseries
python tabstsv_split.py > set1_meth_pattern_bypython_R1.txt
#Remove first column
grep "Sample" set1_meth_pattern_bypython_R1.txt -v > set1_meth_pattern_bypython_R1_re.txt
#---
set1_meth_pattern_bypython_minus_R1 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R1_Sseries/set1_meth_pattern_bypython_R1_re.txt", header = F, stringsAsFactors = F)
head(set1_meth_pattern_bypython_minus_R1)
dim(set1_meth_pattern_bypython_minus_R1)
rownames(set1_meth_pattern_bypython_minus_R1) <- paste0(set1_meth_pattern_bypython_minus_R1$V1, "%",
                                                        set1_meth_pattern_bypython_minus_R1$V2)
set1_meth_pattern_bypython_minus_R1 <- set1_meth_pattern_bypython_minus_R1[,-1]
head(set1_meth_pattern_bypython_minus_R1)
set1_meth_pattern_bypython_minus_R1[set1_meth_pattern_bypython_minus_R1 == "x"] <- NA

#Remove NAs
set1_tab_ref_subNnat_minus_R1_Sseries_py <- na.omit(set1_meth_pattern_bypython_minus_R1)
head(set1_tab_ref_subNnat_minus_R1_Sseries_py)
dim(set1_tab_ref_subNnat_minus_R1_Sseries_py) #Same number as R step
summary(set1_tab_ref_subNnat_minus_R1_Sseries_py)
set1_tab_ref_subNnat_minus_R1_Sseries_py <- set1_tab_ref_subNnat_minus_R1_Sseries_py[,c(1,3:14)]

colnames(set1_tab_ref_subNnat_minus_R1_Sseries_py) <- c("Sample", "CG1r1","CG2r1","CG3r1","CG4r1","CG5r1","CG6r1","CG7r1","CG8r1","CG9r1","CG10r1","CG11r1","CG12r1")
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG1r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG1r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG2r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG2r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG3r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG3r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG4r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG4r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG5r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG5r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG6r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG6r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG7r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG7r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG8r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG8r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG9r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG9r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG10r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG10r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG11r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG11r1)
set1_tab_ref_subNnat_minus_R1_Sseries_py$CG12r1 <- as.numeric(set1_tab_ref_subNnat_minus_R1_Sseries_py$CG12r1)



set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate = aggregate(set1_tab_ref_subNnat_minus_R1_Sseries_py[,c(2:13)],by=list(set1_tab_ref_subNnat_minus_R1_Sseries_py$Sample), mean)
head(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate)
rownames(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate) <- set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate$Group.1
write.table(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate[,2:13], "/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate.txt", col.names =T, sep ="\t", quote = F, append  =F,row.names = T)

pheatmap::pheatmap(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate[,2:13])
colnames(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate) <- c("Sample", "CG1r1","CG2r1","CG3r1","CG4r1","CG5r1","CG6r1","CG7r1","CG8r1","CG9r1","CG10r1","CG11r1","CG12r1")
set1sample_R1_color <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R1_Sseries/sample_color.txt")
colnames(set1sample_R1_color) <- c("Sample","Condition","Combine")
set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col <- merge(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate, set1sample_R1_color, by="Sample")
rownames(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col) <- set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col$Combine
set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col <- set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col[,c(2:13)]
head(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col)
sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col <- stack(as.matrix(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col))
head(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col)
write.table(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col,"/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R1_Sseries/set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col.txt", sep = "\t", quote = F, append = F)
sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col <- data.frame(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col)
sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col <- splitstackshape::cSplit(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col,"row","%")
sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col <- data.frame(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col)
#python3 determine_CpG_position_reverse.py enter-- /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_nospace.fa
set1ref_Nnat_CGsitesR1 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R1_Sseries/CGSites.txt")
colnames(set1ref_Nnat_CGsitesR1) <- c("col", "sites")
sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_cg <- merge(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col, set1ref_Nnat_CGsitesR1, by="col")
dim(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_cg)

ggplot(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_cg, aes(x=as.numeric(sites), y=value, group=row_1))+
  geom_line(aes(color=row_1), size = 1)+
  geom_vline(xintercept = unique(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_cg$sites), colour = "lightblue", linetype="solid")+
  #geom_point(aes(color=row_1), shape=20)+
  theme_classic()+ ylim(0,1)+
  scale_color_manual(values=c(rep("#00CC00",1),rep("red",1)))
ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R1_Sseries/sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_cg.svg", width=25*1.25, height=4*1.25, units="cm", dpi=96)

#Read 2 is also in  reverse as checked in IGV and fasta and BiQ alignmnet both for read and ref, 
#minus ref will be ok to use, assignment to CG position in reverse coordinataes will be Ok
#minus_R2
#Use python to compare results
python tabstsv_split_R2.py > set1_meth_pattern_bypython_R2.txt
#Remove first column
grep "Sample" set1_meth_pattern_bypython_R2.txt -v > set1_meth_pattern_bypython_R2_re.txt

set1_meth_pattern_bypython_minus_R2 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R2_Sseries/set1_meth_pattern_bypython_R2_re.txt", header = F, stringsAsFactors = F)
head(set1_meth_pattern_bypython_minus_R2)
rownames(set1_meth_pattern_bypython_minus_R2) <- paste0(set1_meth_pattern_bypython_minus_R2$V1, "%",
                                                        set1_meth_pattern_bypython_minus_R2$V2)
set1_meth_pattern_bypython_minus_R2 <- set1_meth_pattern_bypython_minus_R2[,-1]
head(set1_meth_pattern_bypython_minus_R2)
dim(set1_meth_pattern_bypython_minus_R2)
set1_meth_pattern_bypython_minus_R2[set1_meth_pattern_bypython_minus_R2 == "x"] <- NA

#Remove NAs
set1_tab_ref_subNnat_minus_R2_Sseries_py <- na.omit(set1_meth_pattern_bypython_minus_R2)
head(set1_tab_ref_subNnat_minus_R2_Sseries_py)
dim(set1_tab_ref_subNnat_minus_R2_Sseries_py) #Same number as R step
summary(set1_tab_ref_subNnat_minus_R2_Sseries_py)
set1_tab_ref_subNnat_minus_R2_Sseries_py <- set1_tab_ref_subNnat_minus_R2_Sseries_py[,c(1,3:15)]

colnames(set1_tab_ref_subNnat_minus_R2_Sseries_py) <- c("Sample", "CG1r2","CG2r2","CG3r2","CG4r2","CG5r2","CG6r2","CG7r2","CG8r2","CG9r2","CG10r2","CG11r2","CG12r2","CG13r2")
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG1r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG1r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG2r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG2r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG3r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG3r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG4r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG4r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG5r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG5r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG6r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG6r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG7r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG7r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG8r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG8r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG9r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG9r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG10r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG10r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG11r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG11r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG12r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG12r2)
set1_tab_ref_subNnat_minus_R2_Sseries_py$CG13r2 <- as.numeric(set1_tab_ref_subNnat_minus_R2_Sseries_py$CG13r2)


dim(set1_tab_ref_subNnat_minus_R2_Sseries_py)
set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate = aggregate(set1_tab_ref_subNnat_minus_R2_Sseries_py[,c(2:14)],by=list(set1_tab_ref_subNnat_minus_R2_Sseries_py$Sample), mean)
head(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate)
rownames(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate) <- set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate$Group.1
pheatmap::pheatmap(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate[,2:14])

write.table(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate[,2:14], "/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate.txt", col.names =T, sep ="\t", quote = F, append  =F,row.names = T)
colnames(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate) <- c("Sample", "CG1r2","CG2r2","CG3r2","CG4r2","CG5r2","CG6r2","CG7r2","CG8r2","CG9r2","CG10r2","CG11r2","CG12r2","CG13r2")
set1sample_R2_color <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R2_Sseries/sample_color.txt")
colnames(set1sample_R2_color) <- c("Sample","Condition","Combine")
set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col <- merge(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate, set1sample_R2_color, by="Sample")
rownames(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col) <- set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col$Combine
set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col <- set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col[,c(2:14)]
head(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col)
sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col <- stack(as.matrix(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col))
head(sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col)
write.table(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col,"/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R2_Sseries/set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col.txt", sep = "\t", quote = F, append = F)
sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col <- data.frame(sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col)
sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col <- splitstackshape::cSplit(sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col,"row","%")
sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col <- data.frame(sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col)
#python3 determine_CpG_position_reverse.py  enter-- /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2_nospace.fa
#tail -6 CGSites.txt > CGSites_filtered.txt
set1ref_Nnat_CGsitesR2 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R2_Sseries/CGSites_filtered.txt")
colnames(set1ref_Nnat_CGsitesR2) <- c("col", "sites")
sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_cg <- merge(sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col, set1ref_Nnat_CGsitesR2, by="col")
dim(sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_cg)
ggplot(sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_cg, aes(x=as.numeric(sites), y=value, group=row_1))+
  geom_line(aes(color=row_2))+
  geom_vline(xintercept = unique(sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_cg$sites), colour = "lightblue", linetype="solid")+
  geom_point(aes(color=row_2), shape=20)+
  theme_classic()+ ylim(0,1)+
  scale_color_manual(values=c("#00CC00","red","blue"))
ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/minus_ref_subNnat_R2_Sseries/sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_cg_merge.svg", width=25*1.25, height=9*1.25, units="cm", dpi=96)

dim(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_cg)
dim(sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_cg)
sset1_tab_ref_subNnat_minusR1nR2_Sseries_py_aggregate_cg <- rbind.data.frame(sset1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_cg,sset1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_cg)
dim(sset1_tab_ref_subNnat_minusR1nR2_Sseries_py_aggregate_cg)
head(sset1_tab_ref_subNnat_minusR1nR2_Sseries_py_aggregate_cg)
ggplot(sset1_tab_ref_subNnat_minusR1nR2_Sseries_py_aggregate_cg, aes(x=as.numeric(sites), y=value, group=row_1))+
  geom_line(aes(color=row_1), size = 1)+
  geom_vline(xintercept = unique(sset1_tab_ref_subNnat_minusR1nR2_Sseries_py_aggregate_cg$sites), colour = "lightblue", linetype="solid")+
  #geom_point(aes(color=row_1), shape=20)+
  theme_classic()+ ylim(0,1)+
  scale_color_manual(values=c(rep("#00CC00",1),rep("red",1)))
ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/sset1_tab_ref_subNnat_minusR1nR2_Sseries_py_aggregate_cg.svg", width=25*1.25, height=9*1.25, units="cm", dpi=96)

paste set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate.txt set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate.txt  > set1_tab_ref_subNnat_minus_R1nR2_Sseries_py_aggregate.txt

#merge set1 and set2

set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg <- rbind.data.frame(sset1_tab_ref_subNnat_minusR1nR2_Sseries_py_aggregate_cg,sset2_tab_ref_subNnat_minusR1nR2_Nseries_py_aggregate_cg)
dim(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg)
head(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg)

#Line chart
ggplot(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg, aes(x=as.numeric(sites), y=value, group=row_1))+
  geom_line(aes(color=row_1), size = 0.5, linetype="solid")+
  geom_vline(xintercept = unique(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg$sites), colour = "lightblue", linetype="solid", size = 0.1)+
  #geom_point(aes(color=row_1), shape=20)+
  theme_classic()+ ylim(0,1)+
  scale_color_manual(values=c(rep("#00CC00",1),rep("red",1),rep("red",29),rep("#00CC00",18),"blue"))
ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/Linechart_set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg.svg", width=25*1.25, height=9*1.25, units="cm", dpi=96)

dim(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg)
dim(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg[which(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg$row_2 == "Downs"),])
dim(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg[which(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg$row_2 == "Control"),])
#Jitter chart
ggplot(set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg, aes(x=row_2, y=value,color=row_2)) + 
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0.3), size = 0.5)+
  theme_classic()+ ylim(0,1)+
  scale_color_manual(values=c("#00CC00","red"))
ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/Jitterchart_set1nset2_tab_ref_subNnat_minusR1nR2_SNseries_py_aggregate_cg.svg", width=9*1.25, height=11*1.25, units="cm", dpi=96)


#Heatmap aggregate
head(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col)
head(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col)
head(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col)
head(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col)

set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col_merge <- set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col
set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col_merge["Sample"] <- rownames(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col_merge)
set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col_merge <- set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col
set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col_merge["Sample"] <- rownames(set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col_merge)

set1_tab_ref_subNnat_minus_R1nR2_Sseries_py_aggregate_col_merge <- merge(set1_tab_ref_subNnat_minus_R1_Sseries_py_aggregate_col_merge, set1_tab_ref_subNnat_minus_R2_Sseries_py_aggregate_col_merge, by ="Sample")

set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col_merge <- set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col
set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col_merge["Sample"] <- rownames(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col_merge)
set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col_merge <- set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col
set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col_merge["Sample"] <- rownames(set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col_merge)

set2_tab_ref_subNnat_minus_R1nR2_Nseries_py_aggregate_col_merge <- merge(set2_tab_ref_subNnat_minus_R1_Nseries_py_aggregate_col_merge, set2_tab_ref_subNnat_minus_R2_Nseries_py_aggregate_col_merge, by ="Sample")

set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge <- rbind.data.frame(set1_tab_ref_subNnat_minus_R1nR2_Sseries_py_aggregate_col_merge,
                                                                                     set2_tab_ref_subNnat_minus_R1nR2_Nseries_py_aggregate_col_merge)


dim(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge)
head(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge)
summary(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge[,c(2:13,21:26)])
rownames(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge) <- set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge$Sample
pheatmap::pheatmap(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge[c(2:49,1),c(2:13,21:26)],
                   color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), 
                   breaks = breaksListd,fontsize = 5,
                   clustering_distance_cols = "euclidean",
                   cluster_rows = F,cluster_cols = T,
                   clustering_method = "ward.D", 
                   border_color = "white",
                   show_rownames = T, show_colnames = T,
                   cellheight = 5,
                   cellwidth=5,legend = T,
                   main= "Average Meth CpGs")
#save as set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge.svg


tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge <- data.frame(t(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge[c(2:49,1),c(2:13,21:26)]))
dim(tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge[,31:49])
tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge1 <- as.matrix(tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge)
head(tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge1[,31:49])
tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge["Allcontrol"] <- rowMedians(tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge1[,31:49])
tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge["AllDownS"] <- rowMedians(tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge1[,1:30])
head(tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge)
dim(tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge)
subtset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge <- tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge[,1:50] - tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge[,50]
head(subtset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge)

breaksListm = seq(-0.5, 0.5, by = 0.001)
pheatmap::pheatmap(t(subtset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge),
                   color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListm)), 
                   breaks = breaksListm,fontsize = 5,
                   clustering_distance_cols = "euclidean",
                   cluster_rows = F,cluster_cols = T,
                   clustering_method = "ward.D", 
                   border_color = "white",
                   show_rownames = T, show_colnames = T,
                   cellheight = 5,
                   cellwidth=5,legend = T,
                   main= "Average Meth CpGs")
#save as subtset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge.svg

#PCA Plot
set1nset2_Nnat_minus_R1nR2_SN_agg_PCA = data.frame(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge[,c(2:13,21:26)])
#write.table(tset1nset2_Nnat_minus_R1nR2_SN_agg_PCA , "tset1nset2_Nnat_minus_R1nR2_SN_agg_PCAdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(set1nset2_Nnat_minus_R1nR2_SN_agg_PCA)
head(set1nset2_Nnat_minus_R1nR2_SN_agg_PCA)

set1nset2_Nnat_minus_R1nR2_SN_agg_PCA["Color"] <-  rownames(set1nset2_Nnat_minus_R1nR2_SN_agg_PCA)
dim(set1nset2_Nnat_minus_R1nR2_SN_agg_PCA)
set1nset2_Nnat_minus_R1nR2_SN_agg_PCA <- splitstackshape::cSplit(set1nset2_Nnat_minus_R1nR2_SN_agg_PCA, "Color", "%")
head(set1nset2_Nnat_minus_R1nR2_SN_agg_PCA)
set1nset2_Nnat_minus_R1nR2_SN_agg_PCA <- data.frame(set1nset2_Nnat_minus_R1nR2_SN_agg_PCA)
set1n2R1n2SNaggdfx <-set1nset2_Nnat_minus_R1nR2_SN_agg_PCA[c(1:(length(set1nset2_Nnat_minus_R1nR2_SN_agg_PCA)-2))]
set1n2R1n2SNaggC <-prcomp(set1n2R1n2SNaggdfx, center = TRUE, scale. = TRUE)
set1n2R1n2SNaggCi<-data.frame(set1n2R1n2SNaggC$x,Color=set1nset2_Nnat_minus_R1nR2_SN_agg_PCA$Color_2)
percentageset1n2R1n2SN <- round(set1n2R1n2SNaggC$sdev^2 / sum(set1n2R1n2SNaggC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentageset1n2R1n2SN <- paste( colnames(set1n2R1n2SNaggCi), "(", paste( as.character(percentageset1n2R1n2SN), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pset1n2R1n2SN <-ggplot(set1n2R1n2SNaggCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentageset1n2R1n2SN[1]) + ylab(percentageset1n2R1n2SN[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("#00CC00","red"))+
  scale_shape_manual(values=c(21,21))
pset1n2R1n2SN <- pset1n2R1n2SN +theme_classic()
pset1n2R1n2SN + xlim(-11,11)+ ylim(-5,5)
ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/PCA_set1nset2_Nnat_minus_R1nR2_SN_agg_PCA_CntrlIndiv.svg", width=9*1.25, height=8*1.25, units="cm", dpi=96)



#Scatter plot
ggplot(tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge, aes(x=Allcontrol, y=AllDownS)) +
  geom_point(size=2, shape=23)+theme_classic()+
  geom_smooth(method=lm, se=FALSE, linetype="dashed")

#save as Scatter_tset1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge.png


#Heatmap individual CpGs
head(set2_tab_ref_subNnat_minus_R1_Nseries_py)
head(set2_tab_ref_subNnat_minus_R2_Nseries_py)
head(set1_tab_ref_subNnat_minus_R1_Sseries_py)
head(set1_tab_ref_subNnat_minus_R2_Sseries_py)

library(gridExtra)
set1nset2_tab_ref_subNnat_minusR1_SNseries_py <- rbind.data.frame(set1_tab_ref_subNnat_minus_R1_Sseries_py, set2_tab_ref_subNnat_minus_R1_Nseries_py)
dim(set1nset2_tab_ref_subNnat_minusR1_SNseries_py)
head(set1nset2_tab_ref_subNnat_minusR1_SNseries_py)
unique(set1nset2_tab_ref_subNnat_minusR1_SNseries_py$Sample)[c(4,6:22)]
breaksListd = seq(0, 1, by = 0.001)


#Set1 and set2 R1 Cases 30
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minusR1_SNseries_py$Sample)[c(4,6:34)]){
  samplecodes <- j
  print(samplecodes)
  tempframe <- set1nset2_tab_ref_subNnat_minusR1_SNseries_py[which(set1nset2_tab_ref_subNnat_minusR1_SNseries_py$Sample == samplecodes),]
  #print(tempframe)
  pheatmapoutput <- pheatmap::pheatmap(tempframe[,2:13],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- pheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=10, nrow=3)
}

#save as set1nset2_tab_ref_subNnat_minusR1_SNseries_py_cases1_30.pdf


#Set1 and set2 R1 Control 19
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minusR1_SNseries_py$Sample)[c(1,35:52)]){
  samplecodes <- j
  print(samplecodes)
  tempframe <- set1nset2_tab_ref_subNnat_minusR1_SNseries_py[which(set1nset2_tab_ref_subNnat_minusR1_SNseries_py$Sample == samplecodes),]
  #print(tempframe)
  pheatmapoutput <- pheatmap::pheatmap(tempframe[,2:13],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- pheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=10, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minusR1_SNseries_py_control1_19.pdf


set1nset2_tab_ref_subNnat_minusR2_SNseries_py <- rbind.data.frame(set1_tab_ref_subNnat_minus_R2_Sseries_py, set2_tab_ref_subNnat_minus_R2_Nseries_py)
dim(set1nset2_tab_ref_subNnat_minusR2_SNseries_py)
head(set1nset2_tab_ref_subNnat_minusR2_SNseries_py)
unique(set1nset2_tab_ref_subNnat_minusR2_SNseries_py$Sample)[c(4,6:34)]
breaksListd = seq(0, 1, by = 0.001)

#Set1 and set2 R2 Cases 1-30
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minusR2_SNseries_py$Sample)[c(4,6:34)]){
  samplecodes <- j
  print(samplecodes)
  tempframe <- set1nset2_tab_ref_subNnat_minusR2_SNseries_py[which(set1nset2_tab_ref_subNnat_minusR2_SNseries_py$Sample == samplecodes),]
  #print(tempframe)
  pheatmapoutput <- pheatmap::pheatmap(tempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- pheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=10, nrow=3)
}

#save as set1nset2_tab_ref_subNnat_minusR2_SNseries_py_cases1_30.pdf


#Set1 and set2 R2 Control 1-19
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minusR2_SNseries_py$Sample)[c(1,35:52)]){
  samplecodes <- j
  print(samplecodes)
  tempframe <- set1nset2_tab_ref_subNnat_minusR2_SNseries_py[which(set1nset2_tab_ref_subNnat_minusR2_SNseries_py$Sample == samplecodes),]
  #print(tempframe)
  pheatmapoutput <- pheatmap::pheatmap(tempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- pheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=10, nrow=2)
}


#save as set1nset2_tab_ref_subNnat_minusR2_SNseries_py_control1_19.pdf

#convert pdf to jpeg https://pdftoimage.com/


#Variants
python3 Het_SNP_detector.py  > Sample_SNPs_37520583.txt

study_participants <- read.delim("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/variants/study_participants.txt", header = F, stringsAsFactors = F)
study_participants <- data.frame(study_participants)
colnames(study_participants) <- c("sample","Sample1","condition","altID","age","gender")
Sample_SNPs_37520583 <- read.delim("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/variants/Sample_SNPs_37520583.txt", header = F, stringsAsFactors = F)
Sample_SNPs_37520583 <- data.frame(Sample_SNPs_37520583)
colnames(Sample_SNPs_37520583) <- c("sample","base1","base2","proportion","zygosity")
study_participants_info <- merge(study_participants, Sample_SNPs_37520583, by ="sample", all.x=T,all.y=T)
head(study_participants_info)
dim(study_participants_info)
writexl::write_xlsx(study_participants_info, "/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/variants/study_participants_info.xlsx")

study_participants_info_het <- study_participants_info[which(study_participants_info$zygosity == "Heterozygous"),]
dim(study_participants_info_het)
head(study_participants_info_het)


#Correlation of methylation with Age and gender
head(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge)
set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re <- data.frame(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge[c(2:49,1),c(2:13,21:26)])
dim(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re)
set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re1 <- as.matrix(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re)
head(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re1)
dim(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re1)
set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re["Medians"] <- rowMedians(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re1[,1:18])
head(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re)
dim(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re)
set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re["sample"] <- rownames(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re)
set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re <- splitstackshape::cSplit(set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re, "sample", "%")
study_participants_info_medianmeth <- merge(study_participants_info, set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re, by.x ="Sample1",by.y ="sample_1", all.x=F,all.y=T)
head(study_participants_info_medianmeth)
dim(study_participants_info_medianmeth)
writexl::write_xlsx(study_participants_info_medianmeth, "/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/study_participants_info_medianmeth.xlsx")
plot(study_participants_info_medianmeth$age, study_participants_info_medianmeth$Medians)
write.table(study_participants_info_medianmeth, "/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/study_participants_info_medianmeth.txt", quote = F, sep = F, append = F)

study_participants_info_medianmeth <- read.delim("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/study_participants_info_medianmeth.txt")
#Scatter plot
ggplot(study_participants_info_medianmeth, aes(x=age, y=Medians, color = condition)) +
  geom_point(size=2, shape=23)+theme_classic()+
  scale_color_manual(values = c("#00CC00","red"))+
  geom_smooth(method=lm, se=FALSE, linetype="dashed")+ylim(0,1)

ggsave("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/Scatter_set1nset2_tab_ref_subNnat_minus_R1nR2_SNseries_py_aggregate_col_merge_re.png", width=9*1.25, height=8*1.25, units="cm", dpi=96)

#####------------------------Allele Specific ---------------------######
#minus_as1B_R2

#Read 2 is also in  reverse as checked in IGV and fasta and BiQ alignmnet both for read and ref, 
#minus_as1B ref will be ok to use, assignment to CG position in reverse coordinataes will be Ok
#minus_as1B_R2
#set1 and set2 analysed in one go
#Use python to compare results
python tabstsv_split_R2.py > set1nset2_meth_pattern_as1B_bypython_R2.txt
#Remove first column
grep "Sample" set1nset2_meth_pattern_as1B_bypython_R2.txt -v > set1nset2_meth_pattern_as1B_bypython_R2_re.txt

set1nset2_meth_pattern_bypython_minus_as1B_R2 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as1B_R2/set1nset2_meth_pattern_as1B_bypython_R2_re.txt", header = F, stringsAsFactors = F)
head(set1nset2_meth_pattern_bypython_minus_as1B_R2)
dim(set1nset2_meth_pattern_bypython_minus_as1B_R2)
rownames(set1nset2_meth_pattern_bypython_minus_as1B_R2) <- paste0(set1nset2_meth_pattern_bypython_minus_as1B_R2$V1, "%",
                                                                  set1nset2_meth_pattern_bypython_minus_as1B_R2$V2)
set1nset2_meth_pattern_bypython_minus_as1B_R2 <- set1nset2_meth_pattern_bypython_minus_as1B_R2[,-1]
head(set1nset2_meth_pattern_bypython_minus_as1B_R2)
set1nset2_meth_pattern_bypython_minus_as1B_R2[set1nset2_meth_pattern_bypython_minus_as1B_R2 == "x"] <- NA

#Remove NAs
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py <- na.omit(set1nset2_meth_pattern_bypython_minus_as1B_R2)
head(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py)
dim(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py) #Same number as R step
summary(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py <- set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py[,c(1,3:15)]
colnames(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py) <- c("Sample", "CG1r2","CG2r2","CG3r2","CG4r2","CG5r2","CG6r2","CG7r2","CG8r2","CG9r2","CG10r2","CG11r2","CG12r2","CG13r2")
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG1r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG1r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG2r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG2r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG3r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG3r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG4r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG4r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG5r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG5r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG6r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG6r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG7r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG7r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG8r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG8r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG9r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG9r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG10r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG10r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG11r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG11r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG12r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG12r2)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG13r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$CG13r2)

dim(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py)
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate = aggregate(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py[,c(2:14)],by=list(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py$Sample), mean)
head(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate)
rownames(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate) <- set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate$Group.1
pheatmap::pheatmap(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate[,2:14])
colnames(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate) <- c("Sample","CG1r2","CG2r2","CG3r2","CG4r2","CG5r2","CG6r2","CG7r2","CG8r2","CG9r2","CG10r2","CG11r2","CG12r2","CG13r2")

#Take only heterozygosity samples
study_participants_info_het_as1B <- data.frame(paste0(study_participants_info_het$Sample1,"_as1B"))
colnames(study_participants_info_het_as1B) <- "Sample"
as1Bnewsample_R2_color <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/sample_color_as1B.txt")
colnames(as1Bnewsample_R2_color) <- c("Sample","Condition","Combine")

as1Bnewsample_R2_color_het <- merge(as1Bnewsample_R2_color, study_participants_info_het_as1B, by="Sample")
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col <- merge(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate, as1Bnewsample_R2_color_het, by="Sample")
rownames(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col) <- set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col$Combine
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col <- set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col[,c(2:14)]
head(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col)
#Remove undetermined
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col <- set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col[c(-48),]
sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col <- stack(as.matrix(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col))
head(sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col)
write.table(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col,"/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as1B_R2/set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col.txt", sep = "\t", quote = F, append = F)
sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col <- data.frame(sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col)
sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col <- splitstackshape::cSplit(sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col,"row","%")
sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col <- data.frame(sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col)
#python3 determine_CpG_position_reverse.py  enter-- /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2_nospace.fa
as1Bref_Nnat_CGsitesR2 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as1B_R2/CGSites.txt")
colnames(as1Bref_Nnat_CGsitesR2) <- c("col", "sites")
sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_cg <- merge(sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_col, as1Bref_Nnat_CGsitesR2, by="col")
dim(sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_cg)



#Set1 and set2 R2 Cases 1-13
set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het <- merge(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py, as1Bnewsample_R2_color_het, by="Sample", all.y=T)
dim(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het)
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het$Sample)[c(1:13)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het_cases1_13.pdf

#Set1 and set2 R2 Control 1-11

plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het$Sample)[c(14:24)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het_control1_11.pdf
N18_S8_as1B <- set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het$Sample == "N18_S8_as1B"),]
pheatmap::pheatmap(N18_S8_as1B[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= "N18_S8")

#minus_as2B_R2

#Read 2 is also in  reverse as checked in IGV and fasta and BiQ alignmnet both for read and ref, 
#minus_as2B ref will be ok to use, assignment to CG position in reverse coordinataes will be Ok
#minus_as2B_R2
#Use python to compare results
python tabstsv_split_R2.py > set1nset2_meth_pattern_as2B_bypython_R2.txt
#Remove first column
grep "Sample" set1nset2_meth_pattern_as2B_bypython_R2.txt -v > set1nset2_meth_pattern_as2B_bypython_R2_re.txt

set1nset2_meth_pattern_bypython_minus_as2B_R2 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as2B_R2/set1nset2_meth_pattern_as2B_bypython_R2_re.txt", header = F, stringsAsFactors = F)
head(set1nset2_meth_pattern_bypython_minus_as2B_R2)
dim(set1nset2_meth_pattern_bypython_minus_as2B_R2)
rownames(set1nset2_meth_pattern_bypython_minus_as2B_R2) <- paste0(set1nset2_meth_pattern_bypython_minus_as2B_R2$V1, "%",
                                                                  set1nset2_meth_pattern_bypython_minus_as2B_R2$V2)
set1nset2_meth_pattern_bypython_minus_as2B_R2 <- set1nset2_meth_pattern_bypython_minus_as2B_R2[,-1]
head(set1nset2_meth_pattern_bypython_minus_as2B_R2)
set1nset2_meth_pattern_bypython_minus_as2B_R2[set1nset2_meth_pattern_bypython_minus_as2B_R2 == "x"] <- NA

#Remove NAs
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py <- na.omit(set1nset2_meth_pattern_bypython_minus_as2B_R2)
head(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py)
dim(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py) #Same number as R step
summary(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py <- set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py[,c(1,3:15)]
colnames(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py) <- c("Sample", "CG1r2","CG2r2","CG3r2","CG4r2","CG5r2","CG6r2","CG7r2","CG8r2","CG9r2","CG10r2","CG11r2","CG12r2","CG13r2")
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG1r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG1r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG2r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG2r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG3r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG3r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG4r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG4r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG5r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG5r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG6r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG6r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG7r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG7r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG8r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG8r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG9r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG9r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG10r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG10r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG11r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG11r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG12r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG12r2)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG13r2 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$CG13r2)

dim(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py)
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate = aggregate(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py[,c(2:14)],by=list(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py$Sample), mean)
head(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate)
rownames(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate) <- set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate$Group.1
pheatmap::pheatmap(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate[,2:14])
colnames(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate) <- c("Sample","CG1r2","CG2r2","CG3r2","CG4r2","CG5r2","CG6r2","CG7r2","CG8r2","CG9r2","CG10r2","CG11r2","CG12r2","CG13r2")

#Take only heterozygosity samples
study_participants_info_het_as2B <- data.frame(paste0(study_participants_info_het$Sample1,"_as2B"))
colnames(study_participants_info_het_as2B) <- "Sample"
as2Bnewsample_R2_color <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/sample_color_as2B.txt")
colnames(as2Bnewsample_R2_color) <- c("Sample","Condition","Combine")

as2Bnewsample_R2_color_het <- merge(as2Bnewsample_R2_color, study_participants_info_het_as2B, by="Sample")
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col <- merge(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate, as2Bnewsample_R2_color_het, by="Sample")
rownames(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col) <- set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col$Combine
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col <- set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col[,c(2:14)]
head(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col)
#Remove undetermined
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col <- set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col[c(-48),]
sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col <- stack(as.matrix(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col))
head(sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col)
write.table(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col,"/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as2B_R2/set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col.txt", sep = "\t", quote = F, append = F)
sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col <- data.frame(sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col)
sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col <- splitstackshape::cSplit(sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col,"row","%")
sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col <- data.frame(sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col)
#python3 determine_CpG_position_reverse.py  enter-- /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2_nospace.fa
as2Bref_Nnat_CGsitesR2 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as2B_R2/CGSites.txt")
colnames(as2Bref_Nnat_CGsitesR2) <- c("col", "sites")
sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_cg <- merge(sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_col, as2Bref_Nnat_CGsitesR2, by="col")
dim(sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_cg)



#Set1 and set2 R2 Cases 1-13
set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het <- merge(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py, as2Bnewsample_R2_color_het, by="Sample", all.y=T)
dim(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het)
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het$Sample)[c(1:13)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het_cases1_13.pdf

#Set1 and set2 R2 Control 1-11

plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het$Sample)[c(14:24)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het_control1_11.pdf


#convert pdf to jpeg https://pdftoimage.com/

head(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het,2)
head(set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het,2)
#convert pdf to jpeg https://pdftoimage.com/

set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het <- rbind.data.frame(set1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_het, 
                 set1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_het)


dim(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het)
set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het <- data.frame(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het)
set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het <- set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het[order(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het$Sample),]
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het$Sample)[c(1:26)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=8, nrow=4)
}

#save as set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het_cases1_13.pdf

#Set1 and set2 R2 Control 1-11

plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het$Sample)[c(14:24)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R2_NnSseries_py_het_control1_11.pdf

#**************-------------------*****************#
#Read 1 is also in  reverse as checked in IGV and fasta and BiQ alignmnet both for read and ref, 
#minus_as1B_R1

#Read 2 is also in  reverse as checked in IGV and fasta and BiQ alignmnet both for read and ref, 
#minus_as1B ref will be ok to use, assignment to CG position in reverse coordinataes will be Ok
#minus_as1B_R1
#set1 and set2 analysed in one go
#Use python to compare results
python tabstsv_split_R1.py > set1nset2_meth_pattern_as1B_bypython_R1.txt
#Remove first column
grep "Sample" set1nset2_meth_pattern_as1B_bypython_R1.txt -v > set1nset2_meth_pattern_as1B_bypython_R1_re.txt

set1nset2_meth_pattern_bypython_minus_as1B_R1 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as1B_R1/set1nset2_meth_pattern_as1B_bypython_R1_re.txt", header = F, stringsAsFactors = F)
head(set1nset2_meth_pattern_bypython_minus_as1B_R1)
dim(set1nset2_meth_pattern_bypython_minus_as1B_R1)
rownames(set1nset2_meth_pattern_bypython_minus_as1B_R1) <- paste0(set1nset2_meth_pattern_bypython_minus_as1B_R1$V1, "%",
                                                                  set1nset2_meth_pattern_bypython_minus_as1B_R1$V2)
set1nset2_meth_pattern_bypython_minus_as1B_R1 <- set1nset2_meth_pattern_bypython_minus_as1B_R1[,-1]
head(set1nset2_meth_pattern_bypython_minus_as1B_R1)
set1nset2_meth_pattern_bypython_minus_as1B_R1[set1nset2_meth_pattern_bypython_minus_as1B_R1 == "x"] <- NA

#Remove NAs
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py <- na.omit(set1nset2_meth_pattern_bypython_minus_as1B_R1)
head(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py)
dim(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py) #Same number as R step
summary(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py <- set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py[,c(1,3:14)]
colnames(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py) <- c("Sample", "CG1r1","CG2r1","CG3r1","CG4r1","CG5r1","CG6r1","CG7r1","CG8r1","CG9r1","CG10r1","CG11r1","CG12r1")
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG1r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG1r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG2r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG2r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG3r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG3r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG4r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG4r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG5r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG5r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG6r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG6r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG7r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG7r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG8r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG8r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG9r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG9r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG10r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG10r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG11r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG11r1)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG12r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$CG12r1)

dim(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py)
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate = aggregate(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py[,c(2:13)],by=list(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py$Sample), mean)
head(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate)
rownames(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate) <- set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate$Group.1
pheatmap::pheatmap(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate[,2:13])
colnames(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate) <- c("Sample","CG1r1","CG2r1","CG3r1","CG4r1","CG5r1","CG6r1","CG7r1","CG8r1","CG9r1","CG10r1","CG11r1","CG12r1")

#Take only heterozygosity samples
study_participants_info_het_as1B <- data.frame(paste0(study_participants_info_het$Sample1,"_as1B"))
colnames(study_participants_info_het_as1B) <- "Sample"
as1Bnewsample_R1_color <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/sample_color_as1B.txt")
colnames(as1Bnewsample_R1_color) <- c("Sample","Condition","Combine")

as1Bnewsample_R1_color_het <- merge(as1Bnewsample_R1_color, study_participants_info_het_as1B, by="Sample")
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col <- merge(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate, as1Bnewsample_R1_color_het, by="Sample")
rownames(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col) <- set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col$Combine
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col <- set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col[,c(2:13)]
head(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col)
#Remove undetermined
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col <- set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col[c(-48),]
sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col <- stack(as.matrix(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col))
head(sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col)
write.table(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col,"/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as1B_R1/set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col.txt", sep = "\t", quote = F, append = F)
sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col <- data.frame(sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col)
sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col <- splitstackshape::cSplit(sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col,"row","%")
sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col <- data.frame(sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col)
#python3 determine_CpG_position_reverse.py  enter-- /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_nospace.fa
as1Bref_Nnat_CGsitesR1 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as1B_R1/CGSites.txt")
colnames(as1Bref_Nnat_CGsitesR1) <- c("col", "sites")
sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_cg <- merge(sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_col, as1Bref_Nnat_CGsitesR1, by="col")
dim(sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_cg)



#Set1 and set2 R1 Cases 1-13
set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het <- merge(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py, as1Bnewsample_R1_color_het, by="Sample", all.y=T)
dim(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het)
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het$Sample)[c(1:24)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het_cases1_13.pdf

#Set1 and set2 R1 Control 1-11

plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het$Sample)[c(14:24)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het_control1_11.pdf

#Read 2 is also in  reverse as checked in IGV and fasta and BiQ alignmnet both for read and ref, 
#minus_as2B ref will be ok to use, assignment to CG position in reverse coordinataes will be Ok
#minus_as2B_R1
#Use python to compare results
python tabstsv_split_R1.py > set1nset2_meth_pattern_as2B_bypython_R1.txt
#Remove first column
grep "Sample" set1nset2_meth_pattern_as2B_bypython_R1.txt -v > set1nset2_meth_pattern_as2B_bypython_R1_re.txt

set1nset2_meth_pattern_bypython_minus_as2B_R1 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as2B_R1/set1nset2_meth_pattern_as2B_bypython_R1_re.txt", header = F, stringsAsFactors = F)
head(set1nset2_meth_pattern_bypython_minus_as2B_R1)
dim(set1nset2_meth_pattern_bypython_minus_as2B_R1)
rownames(set1nset2_meth_pattern_bypython_minus_as2B_R1) <- paste0(set1nset2_meth_pattern_bypython_minus_as2B_R1$V1, "%",
                                                                  set1nset2_meth_pattern_bypython_minus_as2B_R1$V2)
set1nset2_meth_pattern_bypython_minus_as2B_R1 <- set1nset2_meth_pattern_bypython_minus_as2B_R1[,-1]
head(set1nset2_meth_pattern_bypython_minus_as2B_R1)
set1nset2_meth_pattern_bypython_minus_as2B_R1[set1nset2_meth_pattern_bypython_minus_as2B_R1 == "x"] <- NA

#Remove NAs
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py <- na.omit(set1nset2_meth_pattern_bypython_minus_as2B_R1)
head(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py)
dim(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py) #Same number as R step
summary(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py <- set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py[,c(1,3:14)]
colnames(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py) <- c("Sample", "CG1r1","CG2r1","CG3r1","CG4r1","CG5r1","CG6r1","CG7r1","CG8r1","CG9r1","CG10r1","CG11r1","CG12r1")
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG1r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG1r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG2r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG2r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG3r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG3r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG4r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG4r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG5r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG5r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG6r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG6r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG7r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG7r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG8r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG8r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG9r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG9r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG10r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG10r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG11r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG11r1)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG12r1 <- as.numeric(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$CG12r1)

dim(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py)
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate = aggregate(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py[,c(2:13)],by=list(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py$Sample), mean)
head(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate)
rownames(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate) <- set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate$Group.1
pheatmap::pheatmap(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate[,2:13])
colnames(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate) <- c("Sample","CG1r1","CG2r1","CG3r1","CG4r1","CG5r1","CG6r1","CG7r1","CG8r1","CG9r1","CG10r1","CG11r1","CG12r1")

#Take only heterozygosity samples
study_participants_info_het_as2B <- data.frame(paste0(study_participants_info_het$Sample1,"_as2B"))
colnames(study_participants_info_het_as2B) <- "Sample"
as2Bnewsample_R1_color <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/sample_color_as2B.txt")
colnames(as2Bnewsample_R1_color) <- c("Sample","Condition","Combine")

as2Bnewsample_R1_color_het <- merge(as2Bnewsample_R1_color, study_participants_info_het_as2B, by="Sample")
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col <- merge(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate, as2Bnewsample_R1_color_het, by="Sample")
rownames(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col) <- set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col$Combine
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col <- set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col[,c(2:13)]
head(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col)
#Remove undetermined
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col <- set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col[c(-48),]
sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col <- stack(as.matrix(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col))
head(sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col)
write.table(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col,"/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as2B_R1/set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col.txt", sep = "\t", quote = F, append = F)
sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col <- data.frame(sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col)
sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col <- splitstackshape::cSplit(sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col,"row","%")
sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col <- data.frame(sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col)
#python3 determine_CpG_position_reverse.py  enter-- /media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1_nospace.fa
as2Bref_Nnat_CGsitesR1 <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as2B_R1/CGSites.txt")
colnames(as2Bref_Nnat_CGsitesR1) <- c("col", "sites")
sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_cg <- merge(sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_col, as2Bref_Nnat_CGsitesR1, by="col")
dim(sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_cg)



#Set1 and set2 R1 Cases 1-13
set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het <- merge(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py, as2Bnewsample_R1_color_het, by="Sample", all.y=T)
dim(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het)
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het$Sample)[c(1:13)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het_cases1_13.pdf

#Set1 and set2 R1 Control 1-11

plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het$Sample)[c(14:24)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}

#save as set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het_control1_11.pdf


#convert pdf to jpeg https://pdftoimage.com/

head(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het,2)
head(set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het,2)
#convert pdf to jpeg https://pdftoimage.com/

set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het <- rbind.data.frame(set1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_het, 
                                                                                  set1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_het)


dim(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het)
set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het <- data.frame(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het)
set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het <- set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het[order(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het$Sample),]
plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het$Sample)[c(1:26)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=8, nrow=4)
}

#save as set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het_cases1_13.pdf

#Set1 and set2 R1 Control 1-11

plot_list=list()
for (j in unique(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het$Sample)[c(14:24)]){
  samplecodes <- j
  print(samplecodes)
  altempframe <- set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het[which(set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het$Sample == samplecodes),]
  #print(tempframe)
  alpheatmapoutput <- pheatmap::pheatmap(altempframe[,2:14],color = colorRampPalette(c("darkblue", "white", "firebrick3"))(length(breaksListd)), breaks = breaksListd,fontsize = 7,clustering_distance_cols = "euclidean",cluster_rows = F,cluster_cols = F,clustering_method = "ward.D", border_color = "white",show_rownames = F, show_colnames = F,cellwidth = 2,legend = F,labels_col = T, main= samplecodes)
  plot_list[[samplecodes]] <- alpheatmapoutput[[4]]
  grid.arrange(grobs=plot_list, ncol=7, nrow=2)
}
nonsharedR1cpgs <- read.table("/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/minus_ref_subNant_NnS_as2B_R1/nonsharedR1cpgs.txt")
colnames(nonsharedR1cpgs) <- c("CpGs","sites")
nonsharedsset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_cg <- merge(sset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_cg, nonsharedR1cpgs, by.x="sites", all.y =T)
nonsharedsset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_cg <- merge(sset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_cg, nonsharedR1cpgs, by.x="sites", all.y =T)
nonsharedsset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_cg <- nonsharedsset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_cg[,c(2,3,4,5,1)]
nonsharedsset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_cg <- nonsharedsset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_cg[,c(2,3,4,5,1)]
#save as set1nset2_tab_ref_subNnat_minus_as1Bnas2B_R1_NnSseries_py_het_control1_11.pdf
dim(nonsharedsset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_cg)
dim(nonsharedsset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_cg)
dim(sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_cg)
dim(sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_cg)
sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries_py_aggregate_cg <- rbind.data.frame(nonsharedsset1nset2_tab_ref_subNnat_minus_as1B_R1_NnSseries_py_aggregate_cg,
                                                                                             nonsharedsset1nset2_tab_ref_subNnat_minus_as2B_R1_NnSseries_py_aggregate_cg,
                                                                                             sset1nset2_tab_ref_subNnat_minus_as1B_R2_NnSseries_py_aggregate_cg,
                                                                                           sset1nset2_tab_ref_subNnat_minus_as2B_R2_NnSseries_py_aggregate_cg)
library(data.table)
allele_sp_subNnat.case.ids <- c("N15_S7","N18_S8","N23_S13","N24_S14","N27_S15","N33_S19","N34_S20","N37_S23","N39_S25","N40_S26","N42_S28","N43_S29","N5_S4")
sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries.case.ids_list <- list()

for (myids.case in allele_sp_subNnat.case.ids){
  print(myids.case)
  temp3 <- sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries_py_aggregate_cg[sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries_py_aggregate_cg$row_1 %like% myids.case, ]
  temp3p <- ggplot(temp3, aes(x=as.numeric(sites), y=value, group=row_1))+
    geom_line(aes(color=row_1))+
    geom_vline(xintercept = unique(temp3$sites), colour = "gray", linetype="solid")+
    geom_point(aes(color=row_1), shape=20)+
    theme_classic()+ ylim(0,1)+
    scale_color_manual(values=c("darkorange","blue"))
  sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries.case.ids_list[[myids.case]] <- temp3p
}


pdf(file = "/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries.case.ids_list.pdf", height = 10, width = 8)
ggpubr::ggarrange(plotlist = sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries.case.ids_list, ncol=2, nrow=7, common.legend = F, labels=NULL, vjust = 1,hjust=-0.5,font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()


allele_sp_subNnat.control.ids <- c("S54","NC10_S34","NC12_S35","NC14_S37","NC19_S39","NC1_S30","NC3_S31","NC4_S32","NC7_S33","NC_S41","ND_S42")
sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries.control.ids_list <- list()

for (myids.control in allele_sp_subNnat.control.ids){
  print(myids.control)
  temp3 <- sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries_py_aggregate_cg[sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries_py_aggregate_cg$row_1 %like% myids.control, ]
  temp3p <- ggplot(temp3, aes(x=as.numeric(sites), y=value, group=row_1))+
    geom_line(aes(color=row_1))+
    geom_vline(xintercept = unique(temp3$sites), colour = "gray", linetype="solid")+
    geom_point(aes(color=row_1), shape=20)+
    theme_classic()+ ylim(0,1)+
    scale_color_manual(values=c("darkorange","blue"))
  sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries.control.ids_list[[myids.control]] <- temp3p
}


pdf(file = "/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries.control.ids_list.pdf", height = 10, width = 8)
ggpubr::ggarrange(plotlist = sset1nset2_tab_ref_subNnat_minus_as1Bas2B_R1R2_NnSseries.control.ids_list, ncol=2, nrow=7, common.legend = F, labels=NULL, vjust = 1,hjust=-0.5,font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

#convert pdf to jpeg https://pdftoimage.com/



