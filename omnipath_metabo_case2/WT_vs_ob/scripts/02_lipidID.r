#if(!require(devtools)) { install.packages("devtools") }
#install_github("lifs-tools/rgoslin")
# https://github.com/lifs-tools/goslin/blob/master/docs/README.adoc
# isValidLipidName("CoA(8:0)")
library("rgoslin")
library(readxl)

setwd("metabo-usecases/omnipath_metabo_case2/data")
# all lpids
data <- read_excel("ads2547_data_file_s1.xlsx", sheet = 14)
data = as.data.frame(data)

# Update lipid names
data[1:24,2] = gsub("Cer d(.+)_(.+)", "Cer(d\\1/\\2)", data[1:24,2])
data[591:600,2] = gsub("(.+) d(.+)_(.+)", "\\1(d\\2/\\3)", data[591:600,2])
data[668:1305,2] = gsub("TG\\((.+?)\\)\\((.+?)\\)\\((.+?)\\)", "TG(\\1_\\2_\\3)", data[668:1305,2])

#isValidLipidName("TG(12:0/12:0/22:0)[iso3]")
df <- parseLipidNames(data$CompoundName)

# acylcarnitine / acylcoA
data <- read_excel("ads2547_data_file_s1.xlsx", sheet = 15)
data = as.data.frame(data)
data[1:10,1] = gsub("Acylcoenzyme A (.+)", "CoA(\\1)",data[1:10,1])
data[11:19,1] = gsub("Acylcarnitine (.+)", "CAR(\\1)", data[11:19,1])

df <- parseLipidNames(data$CompoundName)
