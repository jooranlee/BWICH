# Libraries required
library(clusterProfiler)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(DESeq2)
library(dplyr)
library(fgsea)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(ggrepel)
library(org.Hs.eg.db)
library(qusage)
library(readr)
library(readxl)
library(reshape)
library(reshape2)
library(SingleCellExperiment)
library(stringr)
library(tibble)
library(TCGAbiolinks)
library(tidyr)
library(UpSetR)
library(RColorBrewer)

##### Fig-S15 
# Define region of interest
chr <- "8"
start_pos <- 127000000
end_pos <- 128200000
samples <- c("A549", "CAKI2", "HCT116", "HepG2", "PANC1", "T47D", "G401", "RPMI_7951")
files <- c("A549_4DNFID68JQY9_score.bedgraph",
           "CAKI2_4DNFIJIPWD63_score.bedgraph",
           "HCT116_4DNFIXTAS6EE_score.bedgraph",
           "HepG2_4DNFICSTCJQZ_score.bedgraph",
           "PANC1_4DNFIFDGVWLU_score.bedgraph",
           "T47D_4DNFITSIPCSK_score.bedgraph", 
           "G401_4DNFI5VO3E1W_score.bedgraph", 
           "RPMI_7951_4DNFIOAA3ZEQ_score.bedgraph")

all_insulation <- data.frame()
for (i in seq_along(files)) {
  df <- read.table(files[i], header = FALSE)
  colnames(df) <- c("chr", "start", "end", "score")
  df <- df[df$chr == chr & df$start >= start_pos & df$end <= end_pos, ]
  df$cell_line <- samples[i]
  all_insulation <- rbind(all_insulation, df)
}

genes <- data.frame(
  gene = c("CCAT1", "MYC", "PVT1"),
  start = c(127207382, 127735434, 127794524),
  end = c(127219268, 127736623, 128101256),
  y = max(all_insulation$score) + 0.3,     
  ymin = -1,                             
  ymax = 2)

custom_colors <- setNames(brewer.pal(length(samples), "Dark2"), samples)
ggplot(all_insulation, aes(x = start, y = score, color = cell_line)) +
  geom_line(size = 1) +
  geom_rect(data = genes, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = gene),
            alpha = 0.2, color = "black")+
  geom_text(data = genes, inherit.aes = FALSE,
            aes(x = (start + end) / 2, y = y, label = gene),
            angle = 0, size = 4, vjust = -0.5) +
  theme_minimal() +
  scale_color_manual(values = custom_colors)+
  labs(title = "Insulation scores from Hi-C",
       x = "Genomic position on Chr8",
       y = "Insulation Index",
       color = "Cell Line",
       fill = "Gene") +
  theme(axis.text.x = element_text(), 
        axis.title = element_text(size = 13), 
        plot.title = element_text(face = 'bold', size = 15, hjust =0.5))


##### Fig-S16 
clinical <- read.csv('clinical2.tsv', sep = '\t')
matrix <- read_table('mRNA_COAD_pan_cancer.txt', col_names = T) %>% as.data.frame()
matrix$CCAT1 <- NULL
mut_matrix <- read.table('mutations.txt', header = T,  sep = "\t") 
mut_matrix$PALS2 <- NULL
mut_matrix_add <- read.table('mutations_2.txt', header = T, sep = "\t") 
mut_matrix <- merge(mut_matrix, mut_matrix_add, by = c('SAMPLE_ID'))
mut_matrix$STUDY_ID.y <- NULL
rownames(mut_matrix) <- mut_matrix$SAMPLE_ID
mut_matrix$SAMPLE_ID <- NULL
mut_matrix$STUDY_ID <- NULL
colnames(mut_matrix) <- c('STUDY_ID', colnames(mut_matrix)[2:58])
cna <- fread("COAD_cna_MYC_DIS3.txt", header = TRUE, data.table = FALSE)
colnames(cna) <- c("STUDY_ID", "SAMPLE_ID" ,"DIS3_cna","MYC_cna","CCAT1_cna","PVT1_cna")
mut_matrix2 <- merge(clinical, cna, by.x = c('Sample.ID') , by.y=c('SAMPLE_ID'))
mut_matrix2 <- merge(mut_matrix2, matrix, by.x = c('Sample.ID'), by.y = c('SAMPLE_ID'))
mut_matrix2$MYC_cna <- replace(mut_matrix2$MYC_cna, mut_matrix2$MYC_cna == '-1', 'Deletion')
mut_matrix2$MYC_cna <- replace(mut_matrix2$MYC_cna, mut_matrix2$MYC_cna == '0', 'Diploid')
mut_matrix2$MYC_cna <- replace(mut_matrix2$MYC_cna, mut_matrix2$MYC_cna == '1', 'Gain')
mut_matrix2$MYC_cna <- replace(mut_matrix2$MYC_cna, mut_matrix2$MYC_cna == '2', 'Amplification')
table(mut_matrix2$MYC_cna)# Define the desired order of the boxes
desired_order <- c("Deletion", "Diploid", "Gain", "Amplification")
mut_matrix2$MYC_cna <- factor(mut_matrix2$MYC_cna, levels = desired_order)
box_colors <- c('Deletion' = 'dodgerblue', 'Diploid' = "#eeb8b4", 'Gain' = "#d44b49", 'Amplification'= "#a21839") 
ggplot(mut_matrix2, aes(x=MYC_cna, y=as.numeric(MYC), MYC_cna)) +
  ylab("MYC: RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)")+
  xlab('MYC: Putative copy-number alterations from GISTIC') + 
  geom_boxplot()  + 
  theme(plot.title = element_text(size = 25, hjust =0.5, face = "bold"), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 10)) + 
  theme_classic()+
  theme(axis.text = element_text(size = 13)) + 
  geom_bracket(xmin=1,xmax=3 , label.size = 4, y.position = 25000,label="p=0.0284",inherit.aes = F)+
  geom_bracket(xmin=1,xmax=4 , label.size = 4, y.position = 32000,label="p=0.0087",inherit.aes = F)+ 
  geom_bracket(xmin=2,xmax=4 , label.size = 4, y.position = 33500,label="p=3.4e-13",inherit.aes = F)+
  geom_bracket(xmin=2,xmax=3, label.size = 4, y.position = 28000,label="p=1.5e-05",inherit.aes = F)


##### Fig-S17 
query <- GDCquery(
  project = c("TCGA-READ", "TCGA-COAD"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification")
GDCdownload(query)
Rnaseq.SE <- GDCprepare(query)
Matrix <- assay(Rnaseq.SE,"unstranded") %>% as.data.frame()

setwd("/Users/jooranlee/Library/CloudStorage/OneDrive-HKUSTConnect/Mphil/Basu_works/")
ccat1_tot <- read.csv('integrated_CCAT1_reads_table.csv')
ccat1_tot$SAMPLE_ID <- substr(ccat1_tot$sample_ID_4, 1, nchar(ccat1_tot$sample_ID_4) -1)
ccat1_tot <- merge(mut_matrix2[,c('Sample.ID', 'MYC_cna')],
                   ccat1_tot, by.x = c('Sample.ID'), by.y = c('SAMPLE_ID'), all.y = T )
ccat1_tot <- ccat1_tot %>% filter(project_ID == 'TCGA-READ' | project_ID == 'TCGA-COAD')
colnames(Matrix) <- substr(colnames(Matrix), 1, nchar(colnames(Matrix))-13)

MYC_exp <- t(Matrix['ENSG00000136997.21', ccat1_tot$Sample.ID ]) %>% as.data.frame()
MYC_exp$SAMPLE_ID <-rownames(MYC_exp)
colnames(MYC_exp) <- c('MYC_exp', 'SAMPLE_ID')
PVT1_exp <- t(Matrix['ENSG00000249859.12',ccat1_tot$Sample.ID ]) %>% as.data.frame()
PVT1_exp$SAMPLE_ID <-rownames(PVT1_exp)
colnames(PVT1_exp) <- c('PVT1_exp', 'SAMPLE_ID')
DIS3_exp <- t(Matrix['ENSG00000083520.15',ccat1_tot$Sample.ID ]) %>% as.data.frame()
DIS3_exp$SAMPLE_ID <-rownames(DIS3_exp)
colnames(DIS3_exp) <- c('DIS3_exp', 'SAMPLE_ID')
ddx21_exp <- t(Matrix['ENSG00000165732.13',ccat1_tot$Sample.ID ]) %>% as.data.frame()
ddx21_exp$SAMPLE_ID <-rownames(ddx21_exp)
colnames(ddx21_exp) <- c('DDX21_exp', 'SAMPLE_ID')
exosc3_exp <- t(Matrix['ENSG00000107371.14',ccat1_tot$Sample.ID ]) %>% as.data.frame()
exosc3_exp$SAMPLE_ID <-rownames(exosc3_exp)
colnames(exosc3_exp) <- c('EXOSC3_exp', 'SAMPLE_ID')

ccat1_tot <- merge(ccat1_tot, MYC_exp, by.x =c('Sample.ID'), by.y = c('SAMPLE_ID'))
ccat1_tot <- merge(ccat1_tot, PVT1_exp, by.x =c('Sample.ID'), by.y = c('SAMPLE_ID'))
ccat1_tot <- merge(ccat1_tot, DIS3_exp, by.x =c('Sample.ID'), by.y = c('SAMPLE_ID'))
ccat1_tot <- merge(ccat1_tot, ddx21_exp, by.x =c('Sample.ID'), by.y = c('SAMPLE_ID'))
ccat1_tot <- merge(ccat1_tot, exosc3_exp, by.x =c('Sample.ID'), by.y = c('SAMPLE_ID'))
ccat1_tot$PVT1_cpm <- ccat1_tot$PVT1_exp/ccat1_tot$reads_tot_mapped * 1000000
ccat1_tot$MYC_cpm <- ccat1_tot$MYC_exp/ccat1_tot$reads_tot_mapped * 1000000
ccat1_tot$DIS3_cpm <- ccat1_tot$DIS3_exp/ccat1_tot$reads_tot_mapped * 1000000
ccat1_tot$DDX21_cpm <- ccat1_tot$DDX21_exp/ccat1_tot$reads_tot_mapped * 1000000
ccat1_tot$EXOSC3_cpm <- ccat1_tot$EXOSC3_exp/ccat1_tot$reads_tot_mapped * 1000000
ccat1_tot <- ccat1_tot[!duplicated(ccat1_tot$sample_ID), ]

mut_matrix2_coad <- mut_matrix2 %>% filter(mut_matrix2$Tumor.Type != "Rectal Adenocarcinoma" |
                                             mut_matrix2$Tumor.Type != 'Rectum Adenocarcinoma' | 
                                             mut_matrix2$Tumor.Type != 'Rectal Adenocarcinoma, Mucinous Type')
ccat1_tot_coad <- ccat1_tot[ccat1_tot$project_ID == 'TCGA-COAD' ,]


scatter2 <- merge(mut_matrix2_coad, ccat1_tot_coad, 
                  by.x = 'Sample.ID', by.y = 'Sample.ID', all.y = TRUE)
scatter2[scatter2$type == 'normal', ]$Sample.ID 
scatter2 <- merge(scatter2, mut, by.x = c('Sample.ID'), by.y = c('SAMPLE_ID'), all.x = TRUE) 

scatter2$MYC_cna.x <- as.character(scatter2$MYC_cna.x)
scatter2$MYC_cna.x[is.na(scatter2$MYC_cna.x)] <- "NA"
scatter2$MYC_cna.x <- as.factor(scatter2$MYC_cna.x)
scatter2$DIS3 <- scatter2$DIS3.y
scatter2$DIS3 <- ifelse(scatter2$DIS3 == 'WT', NA, scatter2$DIS3)
scatter_dip2 <- merge(scatter2, diploid, by.x = c('Sample.ID'), by.y = c('SAMPLE_ID')) 

scatter_final <- rbind(scatter[,c('CCAT1_cpm', 'MYC_cpm', 'PVT1_cpm','type', 'MYC_cna.x', 'DIS3', 'Sample.ID', 'project_ID', 'sample_ID_4')],
                       scatter2[,c('CCAT1_cpm', 'MYC_cpm', 'type',  'PVT1_cpm','MYC_cna.x', 'DIS3', 'Sample.ID', 'project_ID', 'sample_ID_4')] )
scatter_final <- scatter_final[scatter_final$MYC_cna.x == 'Diploid' | 
                                 scatter_final$type == 'normal',]
scatter_final <- scatter_final[!duplicated(scatter_final$sample_ID_4), ]

scatter_long <- scatter_final %>%
  pivot_longer(cols = c(CCAT1_cpm, MYC_cpm, PVT1_cpm),
               names_to = "Gene", values_to = "CPM")
scatter_long$Gene <- recode(scatter_long$Gene,
                            "CCAT1_cpm" = "CCAT1",
                            "PVT1_cpm"  = "PVT1",
                            "MYC_cpm"   = "MYC")

ggplot(scatter_long, aes(x = project_ID, y = log2(CPM+1), fill = project_ID)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test",label.x = 1.5, 
                     size = 5, label = "p.signif") +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_classic(base_size = 14) +
  labs(title = "Expression of CCAT1, MYC, and PVT1",
       x = "", y = "log(CPM+1)") +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none", 
    axis.text =element_text(size = 12, face = "bold")
  )

##### Fig-S18 
complex_map <- tribble(
  ~complex, ~gene,
  "B-WICH complex", "DDX21",
  "B-WICH complex", "BAZ1B",
  "B-WICH complex", "DEK",
  "B-WICH complex", "ERCC6",
  "B-WICH complex", "MYBBP1A",
  "B-WICH complex", "MYO1C",
  "B-WICH complex", "SF3B1",
  "B-WICH complex", "SIRT7",
  "B-WICH complex", "SMARCA5",
  "CBC-ARS complex", "SRRT",
  "CBC-ARS complex", "ZC3H18",
  "Integrator", "INTS1",
  "Integrator", "INTS2",
  "Integrator", "INTS3",
  "Integrator", "INTS4",
  "Integrator", "INTS5",
  "Integrator", "INTS6",
  "Integrator", "INTS7",
  "Integrator", "INTS8",
  "Integrator", "INTS9",
  "Integrator", "INTS10",
  "Integrator", "INTS11",
  "Integrator", "INTS12",
  "HUSH complex", "MPHOSPH8",
  "HUSH complex", "PPHLN1",
  "HUSH complex", "TASOR",
  "PAXT/TRAMP/NEXT complex", "MTREX",
  "NEXT complex", "RBM7",
  "NEXT complex", "ZCCHC8",
  "PAXT complex", "PABPN1",
  "PAXT complex", "ZFC3H1",
  "TRAMP complex", "TENT4B",
  "TRAMP complex", "ZCCHC7",
  "RNA exosome", "DIS3",
  "RNA exosome", "EXOSC1",
  "RNA exosome", "EXOSC2",
  "RNA exosome", "EXOSC3",
  "RNA exosome", "EXOSC4",
  "RNA exosome", "EXOSC5",
  "RNA exosome", "EXOSC6",
  "RNA exosome", "EXOSC7",
  "RNA exosome", "EXOSC8",
  "RNA exosome", "EXOSC9",
  "RNA exosome", "EXOSC10",
  "RNA exosome cofactor", "C1D",
  "RNA exosome cofactor", "MPHOSPH6",
  "RNA exosome cofactor", "SETX",
  "RNA helicase", "AQR",
  "RNA helicase", "DDX1",
  "RNA helicase", "DDX19B",
  "RNA helicase", "DDX23",
  "RNA helicase", "DDX5",
  "RNA helicase", "DHX9",
  "RNase H1", "RNH1",
  "Others", "DGCR8",
  "Others", "PARN",
  "Others", "XRN2"
)
mut_matrix <- read.table('mutations.txt', header = T,  sep = "\t") 
mut_matrix$PALS2 <- NULL
mut_matrix_add <- read.table('mutations_2.txt', header = T, sep = "\t") 
mut_matrix <- merge(mut_matrix, mut_matrix_add, by = c('SAMPLE_ID'))
mut_matrix$STUDY_ID.y <- NULL
rownames(mut_matrix) <- mut_matrix$SAMPLE_ID
mut_matrix$STUDY_ID <- NULL
mut_matrix$SAMPLE_ID <- NULL
mut_matrix$STUDY_ID.x <- NULL

mut_matrix2 <- mut_matrix
mut_matrix2 <- mut_matrix2 %>%
  mutate(across(colnames(mut_matrix2)[1:57], ~ ifelse(. == "WT", 0, 1))) %>% as.data.frame()
mut_matrix2$SAMPLE_ID <- NULL
mut_matrix2$STUDY_ID.x<- NULL
mut_matrix2$STUDY_ID.x<- NULL
mut_matrix2 <- t(mut_matrix2)%>% as.data.frame()
mut_matrix2$sum <- rowSums(mut_matrix2)

mut_matrix3 <- mut_matrix2[,colnames(mut_matrix2) %in% diploid_coad_read$Sample.ID] #203
mut_matrix3$sum <- rowSums(mut_matrix3)

mut_matrix_final <- cbind(mut_matrix2[,c('sum')], mut_matrix3[,c('sum')])
rownames(mut_matrix_final) <- rownames(mut_matrix3)
mut_matrix$SAMPLE_ID <- rownames(mut_matrix)
diploid_coad_read <- merge(diploid_coad_read, mut_matrix, 
                           by.x = c('Sample.ID'), by.y = c('SAMPLE_ID'))
diploid_coad_read$MUT_Counts <- rowSums(diploid_coad_read[, which(names(diploid_coad_read) == "DIS3"):which(names(diploid_coad_read) == "INTS12")] != "WT")
diploid_coad_read <- ccat1_tot %>% filter(MYC_cna == 'Diploid') 

total <- colSums(mut_matrix3[1:203]) %>% as.data.frame()
colnames(total) <- c('sum')
count <- table(total$sum) %>% as.data.frame()
count <- count[1:17,]
ggplot(count, aes(x = Var1, y = Freq)) +
  theme_classic() + 
  geom_bar(stat = "identity", fill = "grey70") +
  coord_flip() + 
  geom_text(aes(label = Freq), hjust = -0.2, size = 5) +
  theme(plot.title = element_text(size = 15, hjust= 0.5, face = 'bold'), 
        axis.text = element_text(size =15), 
        axis.title = element_text(size = 15)) + 
  labs(x = "Number of mutations in RNAsuv genes", y = "Frequency")

##### Fig-4C (Top)
reshape_dis3 <- read.table('reshape_dis3_240709.txt')
reshape_myc <- read.table('reshape_myc_240709.txt')
reshape_pvt1 <- read.table('reshape_pvt1_240709.txt')
reshape_ccat1 <- read.table('reshape_ccat1_240709.txt')

nor_dis3 <- reshape_dis3[reshape_dis3$type == 'normal', ]
nor_myc <- reshape_myc[reshape_myc$type == 'normal', ]
nor_pvt1 <- reshape_pvt1[reshape_pvt1$type == 'normal', ]
nor_ccat1 <- reshape_ccat1[reshape_ccat1$type == 'normal', ]

diploid <- read.csv('all_tcga_cna.txt', sep = '\t')
diploid <- diploid[diploid$MYC == '0', ]

reshape_ccat1 <- rbind(reshape_ccat1[reshape_ccat1$sample_ID %in% diploid$SAMPLE_ID,], nor_ccat1)
reshape_myc <- rbind(reshape_myc[reshape_myc$sample_ID %in% diploid$SAMPLE_ID,], nor_myc)
reshape_pvt1 <- rbind(reshape_pvt1[reshape_pvt1$sample_ID %in% diploid$SAMPLE_ID,],nor_pvt1)
reshape_dis3 <- rbind(reshape_dis3[reshape_dis3$sample_ID %in% diploid$SAMPLE_ID,],nor_dis3)
reshape_ccat1 <- reshape_ccat1[!duplicated(reshape_ccat1), ]

reshape_ccat1$log2 <- log2(reshape_ccat1$value + 1)
reshape_myc$log2 <- log2(reshape_myc$value + 1)
reshape_pvt1$log2 <- log2(reshape_pvt1$value + 1)
reshape_dis3$log2 <- log2(reshape_dis3$value + 1)

reshape_ccat1$project <- substr(reshape_ccat1$project_ID, 6, nchar(reshape_ccat1$project_ID))
reshape_myc$project <- substr(reshape_myc$id_project, 6, nchar(reshape_myc$id_project))
reshape_pvt1$project <- substr(reshape_pvt1$id_project, 6, nchar(reshape_pvt1$id_project))
reshape_dis3$project <- substr(reshape_dis3$id_project, 6, nchar(reshape_dis3$id_project))

colnames(reshape_myc) <- colnames(reshape_ccat1)
colnames(reshape_dis3) <- colnames(reshape_ccat1)
colnames(reshape_pvt1) <- colnames(reshape_ccat1)

mean_myc <- reshape_myc %>% summarise(mean_pts=mean(log2))
reshape_myc_t <- reshape_myc %>% dplyr::filter(type == 'tumor')
ggplot(reshape_myc_t, aes(x = project, y = as.numeric(log2))) +
  ylab("MYC Expression, log2(CPM+1)") +
  xlab('') + 
  geom_boxplot(width = 0.4) + 
  stat_summary(
    fun.data = "mean_sdl", fun.args = list(mult = 1), 
    geom = "pointrange", color = "black"
  ) +
  ggtitle("Comparison of the MYC expressions in TCGA Pan Cancer") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 90),  # <- this line rotates x-axis labels
    axis.title = element_text(size = 10),
    text = element_text(size = 10)
  )

#### Fig-4C (bottom)
reshape_ccat1$log2 <- log2(reshape_ccat1$value + 1)
reshape_ccat1$project <- substr(reshape_ccat1$project_ID, 6, nchar(reshape_ccat1$project_ID))
ggplot(reshape_ccat1, aes(x = project, y = as.numeric(log2))) +
  ylab("CCAT1 Expression, log2(CPM+1)") +
  xlab('') + 
  geom_boxplot(width = 0.4) + 
  stat_summary(
    fun.data = "mean_sdl", fun.args = list(mult = 1), 
    geom = "pointrange", color = "black"
  ) +
  ggtitle("Comparison of the CCAT1 expressions in TCGA Pan Cancer") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 90),  # <- this line rotates x-axis labels
    axis.title = element_text(size = 10),
    text = element_text(size = 10)
  )

##### Fig-4D (Left)
survival_data <- ccat1_tot %>% filter(MYC_cna == 'Diploid')
survival_data <- survival_data[, c("Sample.ID", "case_ID","Osday", "vital_status", "MYC_cna")] 
survival_data$OS_STATUS <- ifelse(survival_data$vital_status == 'Alive', 1, 0)
myc <- read.table('reshape_myc_240709.txt')
survival_data <- merge(survival_data, myc[,c('sample_ID', 'value')], by.x = c('Sample.ID'), by.y = c('sample_ID'))
survival_data <- survival_data[!duplicated(survival_data$Sample.ID), ]
value <- quantile(survival_data$value,  0.5)
survival_data$myc_quant <- ifelse(survival_data$value < value,
                                  'low', 'high')
table(survival_data$myc_quant)
fit <- survfit(Surv(Osday, OS_STATUS) ~ myc_quant, data = survival_data)
p <- ggsurvplot(fit, data = survival_data,
                palette = c("red", "blue"),
                risk.table = FALSE, pval = TRUE)
##### Fig-4D (Right)
ggboxplot(
  survival_data, x = "myc_quant", y = "value", fill = "myc_quant",
  palette = c("low" = "#1f77b4", "high" = "#d62728"),
  add = "jitter"
) +
  stat_compare_means(method = "wilcox.test", label.y = max(survival_data$value) * 1.05) +
  labs(x = NULL, y = "MYC_cpm") +
  theme_classic(base_size = 14)

###### Fig4A
TCGA <- unique(reshape_ccat1[reshape_ccat1$type == 'normal',]$project)
df.fc <- data.frame(project = unique(reshape_ccat1[reshape_ccat1$type == 'normal',]$project))
df.fc$fold_change_CCAT1 <- 0
df.fc$fold_change_MYC <- 0
#CCAT1
for (i in TCGA){
  df <- reshape_ccat1[reshape_ccat1$project == i,]
  fc <- log2(df[df$type == 'tumor',]$value %>% mean()/
               df[df$type == 'normal',]$value %>% mean())
  df.fc[df.fc$project == i,]$fold_change_CCAT1 <- fc
}
#MYC
for (i in TCGA){
  df <- reshape_myc[reshape_myc$project == i,]
  fc <- log2(df[df$type == 'tumor',]$value %>% mean()/
               df[df$type == 'normal',]$value %>% mean())
  df.fc[df.fc$project == i,]$fold_change_MYC <- fc
}
colnames(df.fc) <- c('project', 'CCAT1', 'MYC')
df.melt <- melt(df.fc)
df.melt2 <- df.melt %>% filter(variable %in% c('CCAT1', 'MYC'))
ggplot(df.melt2, aes(x = project, y = value, fill = variable)) +
  geom_col_pattern(
    aes(pattern = variable),
    colour = "black",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    position = position_dodge2(preserve = 'single'),
  ) +theme_classic() +
  scale_fill_manual(
    values = c("MYC" = "#7CAE00", "CCAT1" = "#F8766D"),
    guide = guide_legend(override.aes = list(pattern = "none"))
  ) +
  labs(x = 'TCGA cohort', y = 'log2 Fold Change (tumor/normal)') +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size =10), 
    axis.text.y = element_text(size = 10)
  )

##### Fig4B
reshape_ccat1 <- reshape_ccat1 %>% filter(project_ID == 'TCGA-COAD' |project_ID == 'TCGA-READ')
reshape_myc <- reshape_myc %>% filter(project_ID == 'TCGA-COAD' |project_ID == 'TCGA-READ')
labels <- c("tumor" = "Diploid tumor (n=251)", "normal" = "Normal (n=49)")
g1 <- ggboxplot(reshape_ccat1, x = "type", y = "log2",
                color = "type", palette = c("tumor" = "black", "normal" = "red"),
                add = "jitter",  short.panel.labs = FALSE) + 
  theme(text = element_text(size = 10), 
        plot.title =element_text(hjust = 0.5, face ="bold", size = 12)) + 
  stat_compare_means(label.x  = 1.5) + 
  ggtitle('CCAT1')+ 
  ylab('log2 (CPM+1)')
g2 <- ggboxplot(reshape_myc, x = "type", y = "log2",
                color = "type", palette = c("tumor" = "black", "normal" = "red"),
                add = "jitter",  short.panel.labs = FALSE) + 
  theme(text = element_text(size = 10), 
        plot.title =element_text(hjust = 0.5, face ="bold", size = 12)) + 
  stat_compare_means(label.x  = 1.5)+ 
  ggtitle('MYC')+ 
  ylab('log2 (CPM+1)')
ggarrange(g1,g2, nrow = 1, common.legend = T, legend = "bottom")


##### Fig 4F
Matrix <- read.table('coad_read_table.txt')
ccat1_tot <- read.table('ccat1_tot.txt')
total <- read.table('total_df.csv')
ccat1_tot1 <- ccat1_tot
coad_matrix <- Matrix %>% as.data.frame()
colnames(coad_matrix) <- gsub("\\.", "-", colnames(coad_matrix) )
colnames(coad_matrix)[!colnames(coad_matrix) %in% ccat1_tot$sample_ID]
ensLookup <- gsub("\\.[0-9]*$", "", rownames(coad_matrix))
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db, key = ensLookup, columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
rownames(coad_matrix) <- ensLookup
coad_matrix <- coad_matrix[rownames(coad_matrix) %in% ens2symbol$ENSEMBL,]
ens2symbol <- ens2symbol[!duplicated(ens2symbol$ENSEMBL),]
rownames(ens2symbol) <- ens2symbol$ENSEMBL
ens2symbol <- ens2symbol[rownames(coad_matrix),]
all(ens2symbol$ENSEMBL == rownames(coad_matrix))
coad_matrix$Genes <- ens2symbol$SYMBOL
coad_matrix <- coad_matrix[!duplicated(coad_matrix$Genes),]
colnames(coad_matrix) <- gsub("\\.", "-", colnames(coad_matrix) )
coad_matrix <- na.omit(coad_matrix)
rownames(coad_matrix) <- coad_matrix$Genes
coad_matrix[] <- lapply(coad_matrix, as.numeric)
colData <- data.frame(samples = ccat1_tot$sample_ID, 
                      info = rep('samples', length(ccat1_tot$sample_ID)))
rownames(colData) <- colData$samples
colData$sample_id <-  substr(rownames(colData),1, 15)

normal <- ccat1_tot[ccat1_tot$type == 'normal',]$Sample.ID #49
colData$type <- ifelse(colData$sample_id %in% normal, 'normal', 'tumor')
tumor_diploid <- unique(ccat1_tot[ccat1_tot$MYC_cna=='Diploid', ]$sample_ID)
colData$type <- ifelse(colData$samples %in% tumor_diploid, 'tumor_diploid',colData$type)
colData[colData$type == 'normal',]$samples -> normal_sample_id

coad_matrix <- coad_matrix[,rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = coad_matrix, colData = colData, design = ~ type)
dds <- estimateSizeFactors(dds)
normalized_counts <- vst(dds)
RNAsurv <- read_excel('rna_surv_lists.xlsx')
normalized_counts <- assay(normalized_counts)
normalized_counts <- normalized_counts[rownames(normalized_counts) %in% c(RNAsurv$...4, 'DIS3', 'MYC', 'CCAT1') , ]

colnames(normalized_counts) <- substr(colnames(normalized_counts),1, 15)
no_mut <-rownames(total[total$sum ==0,])
mut <-  rownames(total[total$sum >0,])
normal_rna_surv <- normalized_counts[,colnames(normalized_counts) %in% normal] 
compare_rna_surv <- normalized_counts[,colnames(normalized_counts) %in% no_mut] 
mut_rna_surv <- normalized_counts[,colnames(normalized_counts) %in% mut] 
rna_surv_df <-  t(cbind(mut_rna_surv, normal_rna_surv, compare_rna_surv)) %>% as.data.frame()
rownames(rna_surv_df) <- gsub("\\.", "-", rownames(rna_surv_df) )
rna_surv_df$type <- ifelse(rownames(rna_surv_df) %in% normal, 'normal', 'tumor')
dds3_normal <- colData[colData$sample_id %in% normal,]$samples
dds3_tumor<- colData[colData$sample_id %in% c(mut,no_mut),]$samples
dds3_diploid_no_mut<- colData[colData$sample_id %in% c(no_mut),]$samples
dds3_diploid_mut<- colData[colData$sample_id %in% c(mut),]$samples
dds3 <- DESeqDataSetFromMatrix(countData = coad_matrix[,c(dds3_normal, dds3_tumor)], 
                               colData = colData[c(dds3_normal, dds3_tumor),], design = ~ type)
dds3 <- estimateSizeFactors(dds3)
dds3_counts <- vst(dds3)
dds3 <- DESeq(dds3)
res3 <- results(dds3, name = "type_tumor_diploid_vs_normal") %>% as.data.frame()
res3 <- res3[res3$padj < 0.05,]

significant_genes <- rownames(res3)
length(significant_genes)
ego <- enrichGO(gene = significant_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
ego[grep("ribosome|RNA|ribo", ego$Description),]-> lists
lists <- lists %>% top_n(Count, n=15)
dotplot(ego, showCategory=lists$Description, title="GO Biological Processes", font.size = 11) + 
  theme(plot.title = element_text(size = 13, hjust = 0.5, face = "bold"), 
        axis.title.x = element_text(size = 12)) + theme(
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11))

##### Fig 4E
library(ggrepel)
term_of_interest <- "ribonucleoprotein complex biogenesis"
go_genes <- ego@result %>%
  dplyr::filter(Description == term_of_interest) %>%
  dplyr::pull(geneID) %>%
  strsplit("/") %>%
  unlist()
go_genes_fc <- res3[rownames(res3) %in% go_genes, ]
go_genes_fc$gene <- rownames(go_genes_fc)
go_genes_fc$Significance <- ifelse(go_genes_fc$log2FoldChange > 0, "Upregulated", "Downregulated")

high_fc <- go_genes_fc %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)
low_fc <- go_genes_fc %>%
  arrange(log2FoldChange,desc = F) %>%
  head(3)

upregulated <- sum(go_genes_fc$log2FoldChange > 0)
downregulated <- sum(go_genes_fc$log2FoldChange < 0)
top_genes <- go_genes_fc[order(go_genes_fc$padj),][1:30,] 
top_genes <- rbind(top_genes, high_fc)
go_genes_fc[rownames(go_genes_fc) %in% RNAsurv$Gene,] -> enriched_surv
dim(enriched_surv)
top_genes <- rbind(top_genes, low_fc)
ggplot(go_genes_fc, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 4, color = "black") +  # Adjust label size
  scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "steelblue")) +
  xlab("log2FC") + 
  ylab("-log10(p.adj)") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 15),   
    axis.title.y = element_text(size = 15),   
    axis.text.x  = element_text(size = 15), 
    axis.text.y  = element_text(size = 15), 
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 15))

##### Fig-4G
Matrix<- read.table('TCGA_count_matrix_240301.txt', header = T) 
Matrix$Description <- NULL
Matrix <- Matrix[!duplicated(Matrix$Name),]
Matrix <- Matrix[!duplicated(rownames(Matrix)),]
Matrix <- Matrix[!is.na(Matrix$Name),]
rownames(Matrix) <- Matrix$Name
Matrix$Name <- NULL

coldata <- read.table('GSEA_column_names_TCGA_240301.txt')
rownames(coldata) <- coldata$c
coldata2 <- read.table('reshape_ccat1_240709.txt')
coldata <- merge(coldata, coldata2[,c('sample_ID', 'type')], by= c('sample_ID'))
dis_mut <- read.table('dis3_mut_for_GSEA_240405.txt')
table(dis_mut$dis_mut)
dis_mut$type1 <- ifelse(dis_mut$dis_mut == 'normal', 'normal', 'tumor')
dis_mut <- read_excel('dis3_mut.xlsx') %>% as.data.frame()
dis3 <- dis_mut$`Sample ID`
coldata$mut <- ifelse(coldata$sample_ID %in% dis3,'mut', 'wt')
coldata <- coldata[coldata$type == 'tumor', ]
coldata <- coldata[!duplicated(coldata$sample_ID), ]
rownames(coldata) <- coldata$c 
myc_quant <- read.csv('survival_myc_quant.csv')
coldata <- merge(coldata, myc_quant, by.x = c('sample_ID'), by.y = c('Sample.ID'))
coldata <- coldata[!duplicated(coldata$sample_ID), ]
rownames(coldata) <- coldata$c
Matrix <- Matrix[, rownames(coldata)]
Matrix <- Matrix %>% as.data.frame()

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = Matrix, colData = coldata, design = ~myc_quant)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
rm(keep)
vsd <- vst(dds, blind=FALSE)

#________________DE_analysis_____________#
dds <- DESeq(dds) 
resultsNames(dds) 
res <- results(dds, name="myc_quant_low_vs_high",alpha=0.05)
res <- as.data.frame(res) 
high_id <- coldata[coldata$myc_quant == 'high',] %>% rownames()
low_id <- coldata[coldata$myc_quant == 'low',] %>% rownames()

rna_surv <- read_excel('rna_surv_lists.xlsx', col_names = F) %>% as.data.frame()
colnames(rna_surv) <- c('gene', 'unit')
res$diff <- "NO"
res$diff <- ifelse(res$log2FoldChange>0, 'UP', 'DOWN')
res$delabel <- NA
res[rna_surv$gene,] %>% filter(abs(log2FoldChange) > 0.5) -> genes
res[rownames(genes), ]$delabel <- rownames(res[rownames(genes), ])
res$diff <- ifelse(abs(res$log2FoldChange) < 0.5 & res$padj < 0.05, 'no', res$diff)

pathways.hallmark <- read.gmt('h.all.v2023.2.Hs.symbols.gmt')
ranks <- res$stat
res$SYMBOL <- rownames(res)
names(ranks) <- res$SYMBOL
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) 
gene.in.pathway <- pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest(cols = c(SYMBOL)) %>% 
  inner_join(res, by="SYMBOL")

#______________________VISUALIZATION______________________________#

fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$pval <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
fgseaResTidy <- fgseaResTidy %>% 
  filter(abs(NES) > 1.2)
#find the hallmark gene sets where the MYC gene is involved 
pathways.hallmark <- split(pathways.hallmark$gene, pathways.hallmark$term)
listNames <- character(0)
element <- 'MYC'
for (i in seq_along(pathways.hallmark)) {
  if (element %in% pathways.hallmark[[i]]) {
    listNames <- c(listNames, names(pathways.hallmark)[i])
  }
}
fgseaResTidy$MYC <- ifelse(fgseaResTidy$pathway %in% listNames, 'MYC_related', '')
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  geom_text(aes(label = MYC)) +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")+ 
  guides(color = guide_colorbar(title = "Adjusted P-Value"))

##### Fig-S19B
library(reshape)

genes <- c()
p.values <- c()
fc.values <- c()
for (i in colnames(rna_surv_exp)[1:58]){ 
  k <- rna_surv_exp[rownames(rna_surv_exp) %in% no_mut,i]
  j <- rna_surv_exp[rownames(rna_surv_exp) %in% normal,i]
  fc <- log2(median(as.numeric(k))/median(as.numeric(j)))
  p <- wilcox.test(as.numeric(k), as.numeric(j))$p.value
  genes <- c(genes, i)
  p.values <- c(p.values, p)
  fc.values <- c(fc.values, fc)
}
rna_surv_df <- data.frame(gene = genes, 
                          log2FC = fc.values, 
                          p = p.values)
rna_surv_df$enriched_group <- ifelse(rna_surv_df$log2FC > 0 , 'tumor', 'normal')
rna_surv_df_sig <- rna_surv_df %>% filter(p < 0.05) %>%
  filter(enriched_group == 'normal')
melted_df <- melt(rna_surv_exp, id.vars = 'type', variable.name = colnames(rna_surv_exp)[1:51])  
melted_df  <- melted_df[melted_df$variable %in% rna_surv_df_sig$gene, ]
melted_df$type <- factor(melted_df$type, levels = c("tumor","normal"))
g1 <- ggboxplot(melted_df %>% filter(variable == 'RBM7'), x = "type", y = "value", 
                color = "type",
                add = "jitter",  short.panel.labs = FALSE,
                palette = c("tumor" = "red", "normal" = "black")) + 
  theme(text = element_text(size = 15), 
        plot.title =element_text(hjust = 0.5, 
                                 face ="bold", 
                                 size = 15)) + 
  stat_compare_means(size = 5, , label.x = 1.4, label.y = 8) + 
  ggtitle('RBM7') + 
  labs(y = 'Normalized RNA-seq expression')

g2 <- ggboxplot(melted_df %>% filter(variable == 'SETX'), x = "type", y = "value", 
                color = "type",
                add = "jitter",  short.panel.labs = FALSE,
                palette = c("tumor" = "red", "normal" = "black")) + 
  theme(text = element_text(size = 15), 
        plot.title =element_text(hjust = 0.5, 
                                 face ="bold", 
                                 size = 15)) + 
  stat_compare_means(size = 5, label.x = 1.4, label.y = 12.4) + 
  ggtitle('SETX') + 
  labs(y = 'Normalized RNA-seq expression')


g3 <- ggboxplot(melted_df %>% filter(variable == 'DDX19B'), x = "type", y = "value", 
                color = "type",
                add = "jitter",  short.panel.labs = FALSE,
                palette = c("tumor" = "red", "normal" = "black")) + 
  theme(text = element_text(size = 15), 
        plot.title =element_text(hjust = 0.5, 
                                 face ="bold", 
                                 size = 15)) + 
  stat_compare_means(size = 5, label.x = 1.4, label.y = 7.3) + 
  ggtitle('DDX19B') + 
  labs(y = 'Normalized RNA-seq expression')


g4 <- ggboxplot(melted_df %>% filter(variable == 'C1D'), x = "type", y = "value", 
                color = "type",
                add = "jitter",  short.panel.labs = FALSE,
                palette = c("tumor" = "red", "normal" = "black")) + 
  theme(text = element_text(size = 15), 
        plot.title =element_text(hjust = 0.5, 
                                 face ="bold", 
                                 size = 15)) + 
  stat_compare_means(size = 5, label.x = 1.4, label.y = 7.3) + 
  ggtitle('C1D') + 
  labs(y = 'Normalized RNA-seq expression')

g5 <- ggboxplot(melted_df %>% filter(variable == 'MYO1C'), x = "type", y = "value", 
                color = "type",
                add = "jitter",  short.panel.labs = FALSE,
                palette = c("tumor" = "red", "normal" = "black")) + 
  theme(text = element_text(size = 15), 
        plot.title =element_text(hjust = 0.5, 
                                 face ="bold", 
                                 size = 15)) + 
  stat_compare_means(size = 5, label.x = 1.2, label.y = 11) + 
  ggtitle('MYO1C') + 
  labs(y = 'Normalized RNA-seq expression')
ggarrange(g1, g2, g3 , g4, g5, common.legend =T,nrow = 2)

##### Fig-5D
melted_df2 <- melted_df %>%
  group_by(type, variable) %>%
  summarise(mean = mean(value)) %>% group_by(variable) %>%
  summarise(log2FC = log2(mean[type =='tumor']/mean[type == 'normal']))

ggplot(melted_df2, aes(variable, log2FC)) + 
  geom_bar(position="dodge",stat="identity") + 
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15)) + 
  labs(y = 'log2 FoldChange (Diploid Tumor / Normal)', 
       x = 'Genes') + theme_classic()

##### Fig-5A
no_mut2 <- no_mut[no_mut %in% rownames(rna_surv_exp) ]
mut2 <- mut[mut %in% rownames(rna_surv_exp)] 
samples <- c()
genes2 <- c()
p.values2 <- c()
fc.values2 <- c()
for (i in colnames(rna_surv_exp)[1:58]){ 
  n <- rna_surv_exp[,c(i, 'type')]
  n2 <- n[n$type == 'normal',]
  n2 <- as.numeric(n2[,i]) 
  for (j in c(no_mut2, mut2)) { 
    tumor_exp <- rna_surv_exp[j,i] %>% as.numeric()
    result <- ifelse(tumor_exp < mean(n2), 'downregulated', 'upregulated')
    fc <- log2(tumor_exp/mean(mean(n2)))
    samples <- c(samples, j)
    genes2 <- c(genes2, i)
    fc.values2 <- c(fc.values2, fc)
  }
}
rna_surv_df2 <- data.frame(samples = samples, 
                           gene = genes2, 
                           log2FC = fc.values2)
rna_surv_df2 <- rna_surv_df2[rna_surv_df2$gene != 'MYC' ,]
rna_surv_df2$enriched_type <- ifelse(rna_surv_df2$log2FC > 0, 'tumor', 'normal')
rna_surv_df2 %>%
  group_by(samples, enriched_type )%>%
  summarize(freq = n()) -> summary 

mut_matrix2_grouped <- rna_surv_df2 %>%
  filter(enriched_type == 'normal') %>%
  group_by(gene) %>%
  summarise(sample_count = n(), .groups = "drop") %>%
  left_join(complex_map, by = "gene") %>%
  mutate(
    complex = factor(complex),
    gene = factor(gene, levels = gene[order(complex, sample_count)])
  )
mut_matrix2_grouped %>%
  group_by(complex) %>%
  summarise(mean_sample_count = mean(sample_count), .groups = "drop")
mut_matrix2_grouped$direction <- "down"
mut_matrix3_grouped$direction <- "up"
colnames(mut_matrix2_grouped)[colnames(mut_matrix2_grouped) == "sample_count"] <- "value"
colnames(mut_matrix3_grouped)[colnames(mut_matrix3_grouped) == "sum"] <- "value"
mut_matrix2_grouped$value <- -mut_matrix2_grouped$value
combined_df <- rbind(mut_matrix2_grouped[, c("gene", "value", "complex", "direction")],
                     mut_matrix3_grouped[, c("gene", "value", "complex", "direction")])
gene_order <- combined_df %>%
  group_by(gene, complex) %>%
  summarise(total_freq = sum(abs(value)), .groups = "drop") %>%
  arrange(complex, desc(total_freq)) %>%
  pull(gene) %>%
  unique()
combined_df$gene <- factor(combined_df$gene, levels = gene_order)

ggplot(combined_df, aes(x = gene, y = value, fill = complex)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_text(aes(label = abs(value)), vjust = ifelse(combined_df$direction == "up", -0.3, 1.2), size = 2.5) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                               size = 8),
    plot.title = element_text(size = 13, hjust = 0.5, face = 'bold'),
    legend.title = element_blank()
  ) +
  labs(
    x = "",
    y = "Frequency",
  ) + geom_hline(yintercept = 0, color = "black", linewidth = 0.6) + 
  scale_fill_manual(values = c(
    "B-WICH complex" = "#1f77b4",
    "CBC-ARS complex" = "#ff7f0e",
    "Integrator" = "#2ca02c",
    "HUSH complex" = "#d62728",
    "PAXT/TRAMP/NEXT complex" = "#9467bd",
    "NEXT complex" = "#8c564b",
    "PAXT complex" = "#e377c2",
    "TRAMP complex" = "#7f7f7f",
    "RNA exosome" = "#bcbd22",
    "RNA exosome cofactor" = "#17becf",
    "RNA helicase" = "#aec7e8",
    "RNase H1" = "#ffbb78",
    "RNase H2" = "#98df8a",
    "Others" = "#c5b0d5"
  ))


##### Fig-5B (Top)
df1 <- total[total$sum == 0,]$sample_id
summary1 <- summary %>% filter(samples %in% df1) %>% 
  filter(enriched_type == 'normal') %>%
  as.data.frame()
summary1_freq <- summary1 %>%
  dplyr::mutate(freq = as.integer(as.character(freq))) %>%
  dplyr::count(freq)
ggplot(summary1_freq, aes(x = factor(freq), y = n)) + 
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  theme_classic() +
  labs(
    x = "Number of downregulated RNAsurv genes",
    y = "Number of TCGA samples",
    title = "Samples without any mutations in RNAsurv genes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
##### Fig-5B (Bottom)
df2 <- total[total$sum != 0,]$sample_id
summary2 <- summary %>% filter(samples %in% df2) %>% 
  filter(enriched_type == 'normal') %>%
  as.data.frame()
summary2_freq <- summary2 %>%
  dplyr::mutate(freq = as.integer(as.character(freq))) %>%
  dplyr::count(freq)

ggplot(summary2_freq, aes(x = factor(freq), y = n)) + 
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  theme_classic() +
  labs(
    x = "Number of downregulated RNAsurv genes",
    y = "Number of TCGA samples",
    title = "Samples with at least mutations in RNAsurv genes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

##### Fig-S19A
rna_surv_exp <- cbind(normal_rna_surv, compare_rna_surv, 
                      mut_rna_surv) %>% as.data.frame()
rownames(rna_surv_exp) <- gsub("\\.", "-", rownames(rna_surv_exp))
dim(rna_surv_exp) 
rna_surv_exp <- rna_surv_exp %>% t() %>% as.data.frame()
rownames(rna_surv_exp) <- gsub("\\.", "-", rownames(rna_surv_exp))
rna_surv_exp$type <- 'tumor' 
rna_surv_exp[rownames(rna_surv_exp) %in% normal,]$type <- 'normal'
no_mut2 <- no_mut[no_mut %in% rownames(rna_surv_exp) ] 
mut2 <- mut[mut %in% rownames(rna_surv_exp)] 
norm2 <- normal[normal %in% rownames(rna_surv_exp)] 
samples <- c()
genes2 <- c()
p.values2 <- c()
fc.values2 <- c()

for (i in c('RBM7', 'SETX', 'C1D', 'MYO1C', 'DDX19B')){ 
  n <- rna_surv_exp[c(no_mut2, mut2, norm2),i]
  n2 <- as.numeric(n) 
  for (j in c(no_mut2, mut2)){ 
    exp <- rna_surv_exp[j,i] %>% as.numeric()
    p <- wilcox.test(rna_surv_exp[c(no_mut2, mut2), i], rna_surv_exp[norm2, i])$p.value
    fc <- log2(exp/mean(mean(n2)))
    samples <- c(samples, j)
    genes2 <- c(genes2, i)
    fc.values2 <- c(fc.values2, fc)
    p.values2 <- c(p.values2, p)
  }
}
rna_surv_df2 <- data.frame(samples = samples, 
                           gene = genes2, 
                           p = p.values2,
                           log2FC = fc.values2)
rna_surv_df2 <- rna_surv_df2 %>% filter(p < 0.05)
rna_surv_df2 <- rna_surv_df2 %>% filter(log2FC < 0)
rna_surv_df2$enriched_type <- ifelse(rna_surv_df2$log2FC < 0 , 'downregulated', 'upregulated')
rna_surv_df2 %>%
  group_by(samples, enriched_type)%>%
  summarize(freq = n()) -> summary2
rna_surv_df2 <- rna_surv_df2 %>% filter(enriched_type == 'downregulated')
rna_surv_df3 <- as.data.frame(table(rna_surv_df2$samples, rna_surv_df2$gene))
colnames(rna_surv_df3) <- c("samples", "gene", "count")
rna_surv_df3 <- pivot_wider(rna_surv_df3,
                            names_from = gene,
                            values_from = count,
                            values_fill = 0) %>% as.data.frame()

main_bar_col <- c("violetred4")
sets_bar_col <- c("turquoise4")
matrix_col <- c("slateblue4")
shade_col <- c("wheat4")
text_scale_options1 <- c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5)
upset(rna_surv_df3, nsets = 5,
      sets = unique(rna_surv_df3$gene),
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col,
      order.by = "freq", 
      point.size = 4,  empty.intersections = "on",
      line.size = 1, 
      text.scale = text_scale_options1
)

##### Fig-5C
library(pheatmap)
mut_matrix <- read.table('mutations.txt', header = T,  sep = "\t") 
mut_matrix$PALS2 <- NULL
mut_matrix_add <- read.table('mutations_2.txt', header = T, sep = "\t") 
mut_matrix <- merge(mut_matrix, mut_matrix_add, by = c('SAMPLE_ID'))
mut_matrix$STUDY_ID.y <- NULL
mut_matrix$STUDY_ID.x <- NULL
dim(mut_matrix)
mut_matrix$STUDY_ID <- NULL
rownames(mut_matrix) <- mut_matrix$SAMPLE_ID

mut_matrix <- mut_matrix %>%
  mutate(across(colnames(mut_matrix)[2:58], ~ ifelse(. == "WT", 0, 1))) %>% as.data.frame()
bwich <- c('DDX21', 'MYBBP1A', 'BAZ1B', 'DEK', 'ERCC6', 'MYO1C', 'SF3B1', 'SIRT7', 'SMARCA5')
mut_matrix$sum <- rowSums(mut_matrix[2:58])
rna_surv_df_final<- data.frame(samples = samples, 
                               gene = genes2, 
                               p = p.values2,
                               log2FC = fc.values2)
rna_surv_df_final$significance <- ifelse(rna_surv_df_final$p < 0.05 , 'significant', 'non-significant')
rna_surv_df_final$enriched_type <- ifelse(rna_surv_df_final$log2FC < 0 , 'downregulated', 'upregulated')
mut_down_df <- merge(mut_matrix, rna_surv_df_final, by.x = c('SAMPLE_ID'), by.y = c('samples'))

# STEP 1: Define key genes of interest
key_genes <- c("RBM7", "C1D", "SETX", "MYO1C", 'DDX19B')
# STEP 2: Filter downregulated events only (for heatmap body)
down_df <- mut_down_df %>%
  filter(gene %in% key_genes,
         significance == "significant",
         enriched_type == "downregulated") %>%
  mutate(value = 1) %>%
  dplyr::select(SAMPLE_ID, gene, value)
# STEP 3: Build matrix: rows = genes, columns = samples, values = 0/1 (downregulated or not)
down_matrix <- down_df %>%
  pivot_wider(names_from = SAMPLE_ID, values_from = value, values_fill = list(value = 0)) %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(down_matrix) <- "numeric"
down_matrix[is.na(down_matrix)] <- 0
# STEP 4: Create mutation annotation (for top bar)
mutation_df <- mut_down_df %>%
  group_by(SAMPLE_ID) %>%
  summarise(mut_sum = max(sum)) %>%
  filter(mut_sum > 0) %>%
  mutate(mutated = "Mutated") %>%
  dplyr::select(SAMPLE_ID, mutated) %>%
  unique()
colnames(mutation_df) <-c('SAMPLE_ID', 'Mutation')
# Build annotation data frame across all samples in the matrix
mut_annot <- data.frame(SAMPLE_ID = colnames(down_matrix)) %>%
  left_join(mutation_df, by = "SAMPLE_ID") %>%
  mutate(Mutation = ifelse(is.na(Mutation), "Not Mutated", Mutation)) %>%
  column_to_rownames("SAMPLE_ID")
# STEP 5: Define custom annotation colors
ann_colors <- list(
  Mutation = c("Mutated" = "black", "Not Mutated" = "grey80")
)
# STEP 6: Plot with annotation
pheatmap(down_matrix,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         color = c("grey80", "firebrick"),
         annotation_col = mut_annot,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         fontsize_row = 13,
         main = "Sample clustering by Downregulated RNAsurv Genes")

##### Fig-4E
library(tidyr)

complex_counts <- data.frame()
complex_nondysregulated <- list() 
for (cx in unique(complex_map$complex)) {
  genes <- complex_map %>% filter(complex == cx) %>% pull(gene)
  genes <- intersect(genes, colnames(rna_surv_exp))
  
  if (length(genes) < 2) next
  
  expr_subset <- rna_surv_exp[, genes, drop = FALSE]
  expr_z <- sweep(expr_subset, 2, normal_means_all[genes], FUN = "-") / normal_sd_all[genes]
  
  sample_type <- rna_surv_exp$type
  expr_down <- expr_z < -1
  expr_down[sample_type != "tumor", ] <- FALSE
  has_down <- rowSums(expr_down) >= 1
  non_dysregulated_samples <- rownames(rna_surv_exp)[!has_down & sample_type == "tumor"]
  complex_nondysregulated[[cx]] <- non_dysregulated_samples
  
  n_total <- sum(has_down)
  n_top50 <- sum(has_down[names(has_down) %in% top50_samples])
  
  complex_counts <- rbind(complex_counts, data.frame(
    complex = cx,
    total_dysregulated = n_total,
    top50_dysregulated = n_top50
  ))
}
complex_nondysregulated_df <- stack(complex_nondysregulated)
colnames(complex_nondysregulated_df) <- c("sample", "complex")
complex_nondysregulated_df %>% filter(complex == 'B-WICH complex') -> non_bwich


plot_df <- complex_counts %>%
  pivot_longer(cols = c(total_dysregulated, top50_dysregulated),
               names_to = "category", values_to = "count") %>%
  mutate(category = factor(category, levels = c("total_dysregulated", "top50_dysregulated"),
                           labels = c("All Dysregulated", "Top 50%")))
plot_df <- plot_df %>%
  mutate(percent = count / 215 * 100)
plot_df2 <- plot_df %>% filter(category == 'All Dysregulated')
ggplot(plot_df2, aes(x = reorder(complex, -percent), y = percent, fill = complex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = sprintf("%.0f%%", percent)),
            position = position_dodge(width = 0.8),
            vjust = 0.3, size = 4) +
  ylim(0, 105) + coord_flip() + 
  labs(title = "", y = "Percent of Max Dysregulated", x = "") +
  #scale_fill_manual(values = c("All Dysregulated" = "grey70")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 14), legend.position = "none", 
        plot.title = element_text(face = 'bold', hjust = 0.5))+
  scale_fill_manual(values = c(
    "B-WICH complex" = "#1f77b4",
    "CBC-ARS complex" = "#ff7f0e",
    "Integrator" = "#2ca02c",
    "HUSH complex" = "#d62728",
    "PAXT/TRAMP/NEXT complex" = "#9467bd",
    "NEXT complex" = "#8c564b",
    "PAXT complex" = "#e377c2",
    "TRAMP complex" = "#7f7f7f",
    "RNA exosome" = "#bcbd22",
    "RNA exosome cofactor" = "#17becf",
    "RNA helicase" = "#aec7e8",
    "RNase H1" = "#ffbb78",
    "RNase H2" = "#98df8a",
    "Others" = "#c5b0d5"
  ))

##### Fig-S20
ggplot(plot_df2, aes(x = reorder(complex, -count), y = count, fill = complex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.8), 
            vjust = 0.3, size = 4) + ylim(0,210) + 
  labs(title = "", y = "Counts", x = "") +
  #scale_fill_manual(values = c("All Dysregulated" = "grey70")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3), 
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(), 
        legend.text = element_text(size = 13), legend.position = "none", 
        axis.text = element_text(size = 14), 
        plot.title = element_text(face = 'bold', hjust = 0.5)) +
  scale_fill_manual(values = c(
    "B-WICH complex" = "#1f77b4",
    "CBC-ARS complex" = "#ff7f0e",
    "Integrator" = "#2ca02c",
    "HUSH complex" = "#d62728",
    "PAXT/TRAMP/NEXT complex" = "#9467bd",
    "NEXT complex" = "#8c564b",
    "PAXT complex" = "#e377c2",
    "TRAMP complex" = "#7f7f7f",
    "RNA exosome" = "#bcbd22",
    "RNA exosome cofactor" = "#17becf",
    "RNA helicase" = "#aec7e8",
    "RNase H1" = "#ffbb78",
    "RNase H2" = "#98df8a",
    "Others" = "#c5b0d5"
  ))


##### Fig-S21
library(boot)
complex_summary_list <- list()

for (cx in unique(complex_map$complex)) {
  # Get genes for the complex
  genes <- complex_map %>% filter(complex == cx) %>% pull(gene)
  genes <- intersect(genes, colnames(rna_surv_exp))
  if (length(genes) < 2) next
  
  expr_subset <- rna_surv_exp[, genes, drop = FALSE]
  expr_z <- sweep(expr_subset, 2, normal_means_all[genes], FUN = "-") / normal_sd_all[genes]
  expr_z[rna_surv_exp$type != "tumor", ] <- NA
  down_count <- rowSums(expr_z < -1, na.rm = TRUE)
  common_samples <- intersect(rownames(expr_z), rownames(mut_matrix))
  mut_flag <- rep(NA, nrow(expr_z))
  names(mut_flag) <- rownames(expr_z)
  
  if (length(intersect(genes, colnames(mut_matrix))) > 0) {
    mut_flag[common_samples] <- apply(
      mut_matrix[common_samples, genes, drop = FALSE], 
      1, 
      function(x) any(x != "WT" & x != "" & !is.na(x))
    )
  }
  complex_df <- data.frame(
    Sample = rownames(expr_z),
    Complex = cx,
    Downregulated_Genes = down_count,
    Has_Mutation = mut_flag
  )
  complex_summary_list[[cx]] <- complex_df
}
complex_summary_all <- do.call(rbind, complex_summary_list)
rownames(complex_summary_all) <- NULL
scoring_df <- complex_summary_all %>%
  filter(!is.na(Has_Mutation))
complex_sizes <- complex_map %>%
  group_by(complex) %>%
  summarise(n_genes = n(), .groups = 'drop')
scoring_df <- scoring_df %>%
  left_join(complex_sizes, by = c("Complex" = "complex")) %>%
  mutate(
    norm_down = Downregulated_Genes / n_genes,
    mutation_binary = as.numeric(Has_Mutation)
  )
scoring_df <- scoring_df %>%
  left_join(scatter_final[, c("Sample.ID", "score")], by = c("Sample" = "Sample.ID")) %>%
  mutate(MYC_CCAT1_high = ifelse(score > threshold, 1, 0))
input_df <- scoring_df %>%
  select(MYC_CCAT1_high, norm_down, mutation_binary)
boot_fn <- function(data, indices) {
  d <- data[indices, ]
  fit <- glm(MYC_CCAT1_high ~ norm_down + mutation_binary, data = d, family = "binomial")
  return(coef(fit))
}

set.seed(123)
boot_res <- boot(input_df, boot_fn, R = 1000)
#### scatter plots #### ---
# Step 1: Merge CCAT1 and MYC expression into scoring_df
plot_df <- scoring_df %>%
  left_join(scatter_final[, c("Sample.ID", "CCAT1_cpm", "MYC_cpm")], 
            by = c("Sample" = "Sample.ID")) %>%
  mutate(
    CombinedScore = 0.63 * norm_down + (-0.72) * mutation_binary,  # weights from logistic regression
    log2_CCAT1 = log2(CCAT1_cpm + 1),
    log2_MYC = log2(MYC_cpm + 1)
  ) 
# 2. CCAT1 Faceted Plot
p1 <- ggplot(plot_df, aes(x = log2_CCAT1, y = CombinedScore)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_cor(method = "spearman", aes(label = sprintf("italic(R)==%.2f~italic(p)==%.2g", ..r.., ..p..)), 
           parse = TRUE, size = 4, label.x = 2, label.y = -0.22, color = 'red') +
  facet_wrap(~ Complex, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 11) +
  labs(x = "log2(CCAT1 CPM + 1)", y = "Combined Dysregulation Score") + 
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
# 3. MYC Faceted Plot
p2 <- ggplot(plot_df, aes(x = log2_MYC, y = CombinedScore)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_cor(method = "spearman", aes(label = sprintf("italic(R)==%.2f~italic(p)==%.2g", ..r.., ..p..)), 
           parse = TRUE, size = 4, color = 'red', label.x = 4, label.y = -0.24) +
  facet_wrap(~ Complex, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 11) +
  labs(x = "log2(MYC CPM + 1)", y = "Combined Dysregulation Score")+
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
ggarrange(p1, p2, ncol = 1, labels = c("A", "B"))


##### Fig-5F
# Filter B-WICH only
bwich_df <- plot_df %>% filter(Complex == "B-WICH complex")
# A. Scatter for CCAT1
p1 <- ggplot(bwich_df, aes(x = log2_CCAT1, y = CombinedScore)) +
  geom_point(alpha = 0.7, color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80") +
  stat_cor(method = "spearman", aes(label = sprintf("italic(R)==%.2f~italic(p)==%.2g", ..r.., ..p..)), 
           parse = TRUE, size = 5, color = 'red', label.x = 4, label.y = -0.24) +
  facet_wrap(~ Complex, scales = "free_y", ncol = 4) +
  labs(title = "B-WICH ~ CCAT1", x = expression(log[2]*"(CCAT1 CPM + 1)"), y = "Dysregulation Score") +
  theme_classic(base_size = 14) + theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
# B. Scatter for MYC
p2 <- ggplot(bwich_df, aes(x = log2_MYC, y = CombinedScore)) +
  geom_point(alpha = 0.7, color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80") +
  stat_cor(method = "spearman", aes(label = sprintf("italic(R)==%.2f~italic(p)==%.2g", ..r.., ..p..)), 
           parse = TRUE, size = 5, color = 'red', label.x = 5, label.y = -0.24) +
  facet_wrap(~ Complex, scales = "free_y", ncol = 4) +
  labs(title = "B-WICH ~ MYC", x = expression(log[2]*"(MYC CPM + 1)"), y = "Dysregulation Score") +
  theme_classic(base_size = 14) + 
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))

#RNA-exosome cofactor only 
rna_df <- plot_df %>% filter(Complex == "RNA exosome cofactor")
# A. Scatter for CCAT1
p1 <- ggplot(rna_df, aes(x = log2_CCAT1, y = CombinedScore)) +
  geom_point(alpha = 0.7, color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80") +
  stat_cor(method = "spearman", aes(label = sprintf("italic(R)==%.2f~italic(p)==%.2g", ..r.., ..p..)), 
           parse = TRUE, size = 5, color = 'red', label.x = 4, label.y = -0.24) +
  facet_wrap(~ Complex, scales = "free_y", ncol = 4) +
  labs(title = "RNA exosome cofactor ~ CCAT1", x = expression(log[2]*"(CCAT1 CPM + 1)"), y = "Dysregulation Score") +
  theme_classic(base_size = 14) + theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
# B. Scatter for MYC
p2 <- ggplot(rna_df, aes(x = log2_MYC, y = CombinedScore)) +
  geom_point(alpha = 0.7, color = "black", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80") +
  stat_cor(method = "spearman", aes(label = sprintf("italic(R)==%.2f~italic(p)==%.2g", ..r.., ..p..)), 
           parse = TRUE, size = 5, color = 'red', label.x = 5, label.y = -0.24) +
  facet_wrap(~ Complex, scales = "free_y", ncol = 4) +
  labs(title = "RNA exosome cofactor ~ MYC", x = expression(log[2]*"(MYC CPM + 1)"), y = "Dysregulation Score") +
  theme_classic(base_size = 14) + 
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )
ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))

