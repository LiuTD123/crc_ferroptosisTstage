rm(list = ls())

foldpath <- paste0("D:/workdir/23_Z002-L005/23patientsstatistic")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(GEOquery)
# library(idmap3)
options(timeout = 99999999)
library(readxl)
library(tableone)

# ----------------------------------GSE29621------------------
rm(list = ls())
GSEID <- "GSE29621"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")

# 临床信息
cli <- pData(gset)  #获取临床信息
save(gse29621cli,file = paste0(GSEID,"_cliall.RData"))

gse29621cli <- cli[,c("geo_accession",
                      "characteristics_ch1.1",
                      "characteristics_ch1.2",
                      "characteristics_ch1.6")]

# ----------------------------------GSE72970------------------
rm(list = ls())
GSEID <- "GSE72970"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")

# 临床信息
cli <- pData(gset)  #获取临床信息
save(gse72970cli,file = paste0(GSEID,"_cliall.RData"))

gse72970cli <- cli[,c("geo_accession",
                      "characteristics_ch1.1",
                      "characteristics_ch1.2",
                      "pt:ch1")]

# ----------------------------------GSE17536------------------
rm(list = ls())
GSEID <- "GSE17536"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")

# 临床信息
cli <- pData(gset)  #获取临床信息
save(gse17536cli,file = paste0(GSEID,"_cliall.RData"))

gse17536cli <- cli[,c("geo_accession",
                      "characteristics_ch1",
                      "characteristics_ch1.1",
                      "characteristics_ch1.3")]

# =====================
colnames(gse29621cli) <- c("ID","Gender","T Stage","AJCC Stage")
colnames(gse17536cli) <- c("ID","Age","Gender","AJCC Stage")
colnames(gse72970cli) <- c("ID","Gender","Age","T Stage")

gse17536cli$Age <- as.numeric(substr(gse17536cli$Age,6,999))
gse17536cli$Gender <- substr(gse17536cli$Gender,9,999)
gse17536cli$`AJCC Stage` <- as.numeric(substr(gse17536cli$`AJCC Stage`,13,999))

gse29621cli$Gender <- substr(gse29621cli$Gender,9,999)
gse29621cli$`T Stage` <- as.numeric(substr(gse29621cli$`T Stage`, 10,999))
gse29621cli$`AJCC Stage` <- as.numeric(substr(gse29621cli$`AJCC Stage`,21,999))

gse72970cli$Gender <- substr(gse72970cli$Gender,6,999)
gse72970cli$Age <- as.numeric(substr(gse72970cli$Age,6,999))
gse72970cli$`T Stage` <- substr(gse72970cli$`T Stage`,3,999)

save(gse17536cli,gse29621cli,gse72970cli,file = "geocli.RData")

# ==============
tcgacoad <- read.table("TCGA-COAD.GDC_phenotype.tsv.gz",sep = "\t",header = T,fill = TRUE)
tcgaread <- read.table("TCGA-READ.GDC_phenotype.tsv.gz",sep = "\t",header = T,fill = TRUE)

tcgacoad <- tcgacoad[,c("gender.demographic",
                        "age_at_initial_pathologic_diagnosis",
                        "pathologic_T",
                        "tumor_stage.diagnoses")]
colnames(tcgacoad) <- c("Gender","Age","T Stage","AJCC Stage")

tcgaread <- tcgaread[,c("gender.demographic",
                        "age_at_initial_pathologic_diagnosis",
                        "pathologic_T",
                        "tumor_stage.diagnoses")]
colnames(tcgaread) <- c("Gender","Age","T Stage","AJCC Stage")

tcga <- rbind(tcgacoad,tcgaread)
tcga$Gender <- ifelse(!tcga$Gender%in%c("male","female"),NA,
                      ifelse(tcga$Gender=="male","Male","Female"))
table(tcga$`AJCC Stage`)

tcga$Age <- as.numeric(tcga$Age)

tcga$`T Stage` <- ifelse(!startsWith(tcga$`T Stage`,"T"),NA,
                         substr(tcga$`T Stage`,2,2))

tcga$`AJCC Stage` <- ifelse(!startsWith(tcga$`AJCC Stage`,"stage"),NA,
                         substr(tcga$`AJCC Stage`,7,999))
tcga$`AJCC Stage` <- factor(tcga$`AJCC Stage`, levels = c("i",
                                                          "ii","iia","iib","iic",
                                                          "iii","iiia","iiib","iiic",
                                                          "iv","iva","ivb"),
                            labels = c("i",
                                       "ii","ii","ii","ii",
                                       "iii","iii","iii","iii",
                                       "iv","iv","iv"))
tcga$`AJCC Stage` <- as.numeric(tcga$`AJCC Stage`)

tcga$`T Stage` <- as.numeric(tcga$`T Stage`)
# ================
gse17536cli$ID <- "GSE17536"
gse17536cli$`T Stage` <- NA
gse17536cli$Gender <- ifelse(gse17536cli$Gender=="male","Male","Female")

gse29621cli$ID <- "GSE29621"
gse29621cli$Age <- NA
gse29621cli$Gender <- ifelse(gse29621cli$Gender=="MALE","Male","Female")

gse72970cli$ID <- "GSE72970"
gse72970cli$`AJCC Stage` <- NA
gse72970cli$`T Stage` <- as.numeric(gse72970cli$`T Stage`)

tcga$ID <- "TCGA"

all <- rbind(gse17536cli,
             gse29621cli,
             gse72970cli,
             tcga)
# ==============
outcome <- "ID"

factor_variable <- c("Gender",
                     "AJCC Stage",
                     "T Stage")

all_variable <- colnames(all)[2:5]

table_statistics <- CreateTableOne(vars = all_variable, 
                                   strata = outcome, 
                                   data = all, 
                                   factorVars = factor_variable)

table_statistics <- print(table_statistics, 
                          showAllLevels = TRUE ## 表示展示分类变量所有分类因子的结果
)

write.csv(table_statistics, '03.table_statistics_group_based_all.csv')
