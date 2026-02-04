options(timeout = Inf)

rm(list = ls())
foldpath <- "D:/workdir/23/15hospgroup"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

library(readxl)
library(stringr)

hospdata <- read.csv("hospdata.csv")
table(hospdata$t)

hospdata <- hospdata[!hospdata$t %in% c("δ˵\xc3\xf7",
                                        "δ\xbc\xfb\xb8\xb4\xb7\xa2",
                                        "\xbb\xb9\xc4\xc9\xca\xf5",
                                        "\xceޱ걾",
                                        ""),]
table(hospdata$t)

hospdata$t <- ifelse(hospdata$t %in% c("0","1","2",
                                       "T1","T2"),
                     "Low","High")
colnames(hospdata) <- c("ID","group")
# hospdata <- hospdata[-which(duplicated(hospdata$ID)),]

imagedata <- read_excel("D:/radiomicworkdir/CRC/outall.xlsx")
imagedata <- as.data.frame(imagedata)
colnames(imagedata)[1] <- "ID"

# data <- data[grepl("^(43[0-2]|43[3-4]\\d{1-3}1$|I6[0-4])", data$icd_code), ]

imagedata <- imagedata[!grepl("^TCGA", imagedata$ID), ]

patieninfo <- as.data.frame(imagedata[,"ID"])
colnames(patieninfo) <- "ID"
patieninfo$name <- gsub("[^A-Z|a-z]", "", 
                                patieninfo$ID
                                # substr(patieninfo$ID,5,99999)
                                )
patieninfo$name <- tolower(patieninfo$name)

# # 提取最后一个数字序列
# patieninfo$last_number <- sapply(patieninfo$ID, function(x) {
#   matches <- regmatches(x, gregexpr("-?[0-9]+", x))[[1]]
#   if (length(matches) > 0) {
#     matches[length(matches)]
#   } else {
#     NA
#   }
# })

pyinfo <- read.csv("pinfo.csv")
# pyinfo <- pyinfo[-which(duplicated(pyinfo$ID)),]
pyinfo <- merge(pyinfo,patieninfo,by="name")
colnames(pyinfo)[2:3] <- c("numbid","ID")

data <- merge(pyinfo,imagedata,by = "ID")
data <- data[,-c(1,2)]
colnames(data)[1] <- "ID"

data2 <- merge(hospdata,data,by = "ID")
data <- data2[-which(duplicated(data2$ID)),]

save(data,file = "data.RData")
