##### R file created on July 2021 by Martina Reichmuth
#### National sequencing surveillance program on SARS-CoV-2 variants

### install packages
list.of.packages <- c("lubridate","downloader","xlsx","httr","nnet","effects","effects","splines","emmeans","reshape2","ggplot2","RColorBrewer","gridExtra","grid","tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
### load libraries:
library(lubridate)
library(downloader)
library(xlsx)
library(httr)
library(nnet)
library(effects)#nnet if effects()
library(splines)# if ns()
library(emmeans)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(tidyverse)



### load data:

## Swiss population size:
#url <- "https://www.bag.admin.ch/dam/bag/de/dokumente/mt/k-und-i/aktuelle-ausbrueche-pandemien/2019-nCoV/covid-19-basisdaten-bevoelkerungszahlen.xlsx.download.xlsx/Population_Size_BFS.xlsx"
#GET(url, write_disk(tf <- tempfile(fileext = ".xls")))
#CH_canton_age_population <- read.xlsx(tf,sheetIndex = 2, startRow = 0)
#colnames(CH_canton_age_population) <- c("canton", "gender", "age", "no_pop")
#remove(url)
cantons_ch <- c("CH", "AG","AI","AR","BE","BL","BS","FR","GE",
                "GL","GR", "JU","LU","NE","NW","OW","SG","SH",
                "SO","SZ","TG",  "TI","UR","VD","VS","ZG","ZH")

## Swiss SARS-CoV-2  metadata:
url <- readLines("https://www.covid19.admin.ch/api/data/context", warn=FALSE)
url <- gsub("^(.*)https:", "https:", url[20])
url <- gsub("\"", "",url)
url <- gsub(",", "",url)
download(url, dest=paste0("covid19_bag_", Sys.Date(),".zip"), mode="wb") 
unzip(paste0("covid19_bag_", Sys.Date(),".zip"), exdir = "./temp_data")


BAG_cases_canton <- read.csv("./temp_data/data/COVID19Cases_geoRegion.csv")
colnames(BAG_cases_canton)[2] <- "date"
colnames(BAG_cases_canton)[3] <- "cases_num"
BAG_cases_canton <- BAG_cases_canton[BAG_cases_canton$geoRegio %in% cantons_ch,]

BAG_test_canton <- read.csv("./temp_data/data/COVID19Test_geoRegion_PCR_Antigen.csv")
colnames(BAG_test_canton)[1] <- "date"
BAG_test_canton <- BAG_test_canton[BAG_test_canton$geoRegio %in% cantons_ch,]

BAG_test_canton_pcr <- BAG_test_canton[BAG_test_canton$nachweismethode=="PCR",]
colnames(BAG_test_canton_pcr)[2] <- "pcrtests_num"
colnames(BAG_test_canton_pcr)[3] <- "pcrtests_pos_num"

BAG_test_canton_antig <- BAG_test_canton[BAG_test_canton$nachweismethode=="Antigen_Schnelltest",]
colnames(BAG_test_canton_antig)[2] <- "antigtests_num"
colnames(BAG_test_canton_antig)[3] <- "antigtests_pos_num"
geoR <- grep("geoRegion",colnames(BAG_test_canton))
BAG_test_canton <- merge(BAG_test_canton_antig[,c(1,2,3,geoR)], BAG_test_canton_pcr[,c(1,2,3,geoR)],by=c("date","geoRegion"), all=TRUE )
BAG_test_canton<- BAG_test_canton %>% 
  rowwise() %>% #rowwise will make sure the sum operation will occur on each row
  mutate(tests_num = sum(antigtests_num,pcrtests_num, na.rm=TRUE))%>% 
  mutate(tests_pos_num = sum(antigtests_pos_num,pcrtests_pos_num, na.rm=TRUE))
BAG_test_canton[is.na(BAG_test_canton)] <- 0

BAG_Re_canton <- read.csv("./temp_data/data/COVID19Re_geoRegion.csv")
BAG_Re_canton <- BAG_Re_canton[BAG_Re_canton$geoRegio %in% cantons_ch,]

BAG_data <- c()
BAG_data <- merge(BAG_cases_canton[,c("date","geoRegion","cases_num", "pop")], BAG_test_canton[,c("date","geoRegion","tests_num", "tests_pos_num")],by=c("date","geoRegion"), all=TRUE )
BAG_data <- merge(BAG_data, BAG_Re_canton[,c("date","geoRegion", "median_R_mean", "median_R_highHPD", "median_R_lowHPD")],by=c("date","geoRegion"), all=TRUE )
BAG_data <- BAG_data[!is.na(BAG_data$cases_num),]

unlink("temp_data", recursive = TRUE)
unlink(paste0("covid19_bag_", Sys.Date(),".zip"), recursive = TRUE)
remove(BAG_cases_canton)
remove(BAG_test_canton)
remove(BAG_Re_canton)
remove(BAG_test_canton_antig)
remove(BAG_test_canton_pcr)

## Swiss SARS-CoV-2 sequencing metadata:
#variants_ch <- read.csv("https://raw.githubusercontent.com/CovSeqCH/SequencingStats_and_Variants/cr-dev/data/cases_seq_by_day_region.csv")

# new to get sequence data incl. canton:
#url <- GET("https://cov-spectrum.ethz.ch/api/resource/sample2?country=Switzerland&fields=date,division,pangolinLineage")
#url <- GET("https://cov-spectrum.ethz.ch/gisaid/api/v1/sample/aggregated?country=Switzerland&fields=date,division,pangoLineage&accessKey=9Cb3CqmrFnVjO3XCxQLO6gUnKPd")#Used since 17 Nov 2021

url <- GET("https://lapis.cov-spectrum.org/open/v1/sample/aggregated?country=Switzerland&fields=date,division,pangoLineage")#Used since Oct 2022
jsonRespParsed<- content(url,as="parsed", encoding="UTF-8") 
seq_ch <- jsonRespParsed%>%bind_rows#%>%select(date,division,pangolinLineage)# %>%subset(.,country %in% "Switzerland") #%>%
seq_ch <- seq_ch[,c("date","division","pangoLineage","count")]
#seq_ch1 <- jsonRespParsed%>%bind_rows%>%select(date,division,pangolinLineage)# %>%subset(.,country %in% "Switzerland") #%>%

seq_ch <- seq_ch[seq_ch$count>0&!is.na(seq_ch$count),]
seq_ch$country <- "CH"
seq_ch <- seq_ch[rep(row.names(seq_ch), seq_ch$count), c(1,2,3,5)]
#seq_ch <- seq_ch[!is.na(seq_ch$pangolinLineage),]

remove(url)
remove(jsonRespParsed)

### prepare data / data cleaning:
month_end <- as.numeric(format(Sys.Date(),"%m"))
year_end <- as.numeric(format(Sys.Date(),"%Y"))
if(month_end+1>9 &month_end+1<13){
  time_window <- c(as_date("2022-01-01"), as_date(paste0(year_end, month_end+1,"-01"))-1)
  
}
if(month_end==1 &year_end==2022){
  time_window <- c(as_date("2022-01-01"), as_date(paste0("2022-01-01"))-1)
  
}
if(month_end+1<10){
  time_window <- c(as_date("2022-01-01"), as_date(paste0(year_end,"-", month_end+1,"-01"))-1)
  
}
#time_window <- c(as_date("2021-01-01"), as_date(paste0("2021-", month_end+1,"-01"))-1)
month_start <- format(as.numeric(format(Sys.Date(),"%m"))-1, format="%m")
if(month_end==1){
  month_start <- 12
}
month_end <- format(month_end, format="%m")
period <- function(x) {
  c(as_date(format(x, paste0("%Y-",month_start,"-01"))),
    lubridate::floor_date(as_date(format(x, paste0("%Y-",month_end,"-01"))), unit = "month") - 1)
}
period_date <- period(Sys.Date())
if(month_start==12){
  period_date <- c(as_date("2021-12-01"), as_date("2021-12-31"))
}
period_days <- seq(period_date[1], period_date[2],1)

if(period_date[2]=="2021-09-30"){ #for September
period_date <- c(as_date("2021-08-30"), as_date("2021-10-03"))
period_days <- seq(period_date[1], period_date[2],1)
}

# after September starting looking at weeks not 1st to last date of month:
if(weekdays( period_date[1])== "Monday"){
  period_date <- c(period_date[1], floor_date(as.Date(period_date[2], "%m/%d/%Y"), unit="week"))
  
}
if(weekdays( period_date[1])!= "Monday"){
  period_date <- c(ceiling_date(as.Date(period_date[1], "%m/%d/%Y"), unit="week"), floor_date(as.Date(period_date[2], "%m/%d/%Y"), unit="week"))
  
}
period_days <- seq(period_date[1], period_date[2],1)


time_window[2] <- time_window[2]-1
#variants_ch <- subset(variants_ch, as_date(date) %in% seq(time_window[1],time_window[2],1))#Sys.Date()-14
BAG_data <- subset(BAG_data, as_date(date) %in% seq(time_window[1],time_window[2],1))
seq_ch <- subset(seq_ch, as_date(date) %in% seq(time_window[1],time_window[2],1))


### renaming functions / variables generation
## Renaming variants according to WHO
#variants_ch <- melt(variants_ch, id.vars=c("date", "region","cases","sequences"))

who_variant_names <- function(x){
  if(is.na(x)){return("undetermined")}
  else if(grepl("Alpha|alpha|B.1.1.7|Q.1|Q.2|Q.3|Q.4|Q.6",x,useBytes = TRUE)){return("Alpha")}
  else if(grepl("Beta|beta|B.1.351|B.1.351.1|B.1.351.2",x,useBytes = TRUE)){return("Beta")}
  else if(grepl("Gamma|gamma|P.1|P.1.1|P.1.2|P.1.3|P.1.4|P.1.5|P.1.6",x,useBytes = TRUE)){return("Gamma")}
  else if(grepl("Delta|delta|B.1.617.2|AY.1.1|AY.2|AY.3|AY.3.1|AY.4|AY.5|AY.5.1|AY.5.2|AY.6|AY.7|AY.7.1|AY.7.2|AY.8|AY.9|AY.10|AY.11|AY.12|AY.13|AY.14|AY.15|AY.16|AY.17|AY.18|AY.19|AY.20|AY.21|AY.22|AY.23|AY.24|AY.25|AY.26|AY.27|AY.28|AY.29|AY.30|AY.31|AY.32|AY.33|AY.34|AY.35|AY.36|AY.37",x,useBytes = TRUE)){return("Delta")}
  else if(grepl("C.37|Lambda|lambda",x)){return("Lambda")}
  #else if(grepl("C.36",x)){return("C.36*")}
  else if(grepl("Mu|mu|B.1.621|B.1.621.1|B.1.621.2|B.1.621.3",x,useBytes = TRUE)){return("Mu")}#useBytes = FALSE
  else if(grepl("B.1.1.318|AZ.2|AZ.",x)){return("B.1.1.318")}#,useBytes = FALSE
  #else if(grepl("B.1.1.529|BA.1|BA.2",x)){return("Omicron")}#,useBytes = FALSE
  else if(grepl("BA.2.75",x)){return("Omicron (BA.2.75)")}#,useBytes = FALSE
  else if(grepl("BA.2",x)){return("Omicron (BA.2)")}#,useBytes = FALSE
  else if(grepl("BA.1",x)){return("Omicron (BA.1)")}#,useBytes = FALSE
  else if(grepl("BA.3",x)){return("Omicron (BA.3)")}#,useBytes = FALSE
  else if(grepl("BA.4",x)){return("Omicron (BA.4)")}#,useBytes = FALSE
  else if(grepl("BA.5",x)){return("Omicron (BA.5)")}#,useBytes = FALSE
  else if(grepl("XE",x)){return("Omicron (BA.1 & BA.2)")}#,useBytes = FALSE
  else if(grepl("BA.2.12.1",x)){return("Omicron (BA.2.12.1)")}#,useBytes = FALSE
  else if(grepl("BQ.1|BQ.1.1|BQ.1.2|BQ.1.2|BQ.1.3|BQ.1.4",x)){return("Omicron (BQ.1)")}#,useBytes = FALSE
  else if(grepl("XBB|XBB.1",x)){return("Omicron (XBB)")}#,useBytes = FALSE
  else if(grepl("Unassigned",x)){return("undetermined")}#,useBytes = FALSE
  else{return("others")}
  #else if(x =="others"){return("others")}
  #else{return(x)} 
}
#https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/
#variants_ch$who_variants <- sapply(variants_ch$variable, who_variant_names)
seq_ch$who_variants <- sapply(seq_ch$pangoLineage, who_variant_names)#seq_ch$pangolinLineage
#lev <- c("Alpha",  "Beta",  "Gamma", "Delta","Lambda","Mu", "B.1.1.318", "Omicron", "others", "undetermined")

lev <- c("Alpha",  "Beta",  "Gamma", "Delta","Lambda","Mu", "B.1.1.318", "Omicron (BA.1)","Omicron (BA.2)", "Omicron (BA.3)","Omicron (BA.4)","Omicron (BA.5)","Omicron (BA.2.12.1)","Omicron (BA.1 & BA.2)", "Omicron (BA.2.75)","Omicron (BQ.1)","Omicron (XBB)","others", "undetermined")
seq_ch$who_variants <- factor(seq_ch$who_variants, levels = lev)
#variants_ch$variable <- NULL

## Divide Switzerland in 6 regions as following:
## region 1 «GE, NE, VD, VS»
## region 2 «BE, FR, JU»
## region 3 «AG, BL, BS, SO»
## region 4 «LU, NW, OW, SZ, UR, ZG»
## region 5 «AI, AR, GL, SG, SH, TG, ZH»
## region 6 «GR, TI».
from_cantons_to_abrev <- function(x){
  if(grepl("Aargau|Argovie|AG|Ag",x)){return("AG")}
  else if(grepl("Basel-Land|BL|Bl",x)){return("BL")}
  else if(grepl("Basel-Stadt|Basel|BS|Bs",x)){return("BS")}
  else if(grepl("Fribourg|Freiburg|FR|Ft",x)){return("FR")}
  else if(grepl("Genf|Geneva|Geneva|GE|Ge",x)){return("GE")}
  else if(grepl("Glarus|GL|Gl",x)){return("GL")}
  else if(grepl("Jura|JU|Ju",x)){return("JU")}
  else if(grepl("Lucerne|Luzern|LU|Lu",x)){return("LU")}
  else if(grepl("Nidwalden|NW|Nw",x)){return("NW")}
  else if(grepl("Obwalden|Obwald|OW|Ow",x)){return("OW")}
  else if(grepl("Schaffhausen|Schaffhouse|SH|Sh",x)){return("SH")}
  else if(grepl("Solothurn|Olten|SO|So",x)){return("SO")}
  else if(grepl("Schwyz|Schwytz|SZ|Sz",x)){return("SZ")}
  else if(grepl("Thurgau|Turgovia|Thurgovie|TG|Tg",x)){return("TG")}
  else if(grepl("Ticino|Tessin|TI|Ti",x)){return("TI")}
  else if(grepl("Uri|UR|Ur",x)){return("UR")}
  else if(grepl("Vaud|Waadt|VD|Vd|Lausanne|Yverdon-les-Bains",x)){return("VD")}
  else if(grepl("Valais|Wallis|VALAIS|VS|Vs|Sion",x)){return("VS")}
  else if(grepl("Zug|ZG|Zg|Zoug",x)){return("ZG")}
  else if(grepl("Zürich|Zurich|Zaerich|Zoerich|Zuerich|ZH|Zh",x)){return("ZH")}
  else if(grepl("Neuchâtel|Neuenburg|Neuchã¢Tel|NE|Ne",x)){return("NE")}#Ne brings problem
  else if(grepl("Appenzell Innerrrhoden|AI|Ai",x)){return("AI")}#AI brings problem due to VALAIS, Appenzell is in AI
  else if(grepl("Appenzell Ausserrhoden|Appenzell-Ausserrhoden|AR|Ar",x)){return("AR")}#Ar brings problem, attention not SG
  else if(grepl("Appenzell",x)){return("AI")}# Appenzell is in AI
  else if(grepl("Sankt Gallen|St. Gallen|St.Gallen|Saint-Gallen|Saint-Gall|SG|Sg|St Gall",x)){return("SG")}
  else if(grepl("Bern|BE|Be",x)){return("BE")}# Be brings problem , attention not SG
  else if(grepl("Graubünden|Graubã¼Nden|Grisons|Graubunden|GR|Gr",x)){return("GR")}#Gr brings problem
  else {return("Unknown")} 
}
seq_ch$canton  <- sapply(seq_ch$division, from_cantons_to_abrev)
#table(seq_ch$division[seq_ch$canton=="Unknown"])
#table(seq_ch$division,seq_ch$canton)
from_cantons_to_regions <- function(x){
  if(x== "AG"){return("region_3")}
  else if(x== "AI"){return("region_5")}
  else if(x== "AR"){return("region_5")}
  else if(x== "BE"){return("region_2")}
  else if(x== "BL"){return("region_3")}
  else if(x== "BS"){return("region_3")}
  else if(x== "FR"){return("region_2")}
  else if(x== "GE"){return("region_1")}
  else if(x== "GL"){return("region_5")}
  else if(x== "GR"){return("region_6")}
  else if(x== "JU"){return("region_2")}
  else if(x== "LU"){return("region_4")}
  else if(x== "NE"){return("region_1")}
  else if(x== "NW"){return("region_4")}
  else if(x== "OW"){return("region_4")}
  else if(x== "SG"){return("region_5")}
  else if(x== "SH"){return("region_5")}
  else if(x== "SO"){return("region_3")}
  else if(x== "SZ"){return("region_4")}
  else if(x== "TG"){return("region_5")}
  else if(x== "TI"){return("region_6")}
  else if(x== "UR"){return("region_4")}
  else if(x== "VD"){return("region_1")}
  else if(x== "VS"){return("region_1")}
  else if(x== "ZG"){return("region_4")}
  else if(x== "ZH"){return("region_5")}
  else {return("Unknown")} 
}
seq_ch$region  <- sapply(seq_ch$canton, from_cantons_to_regions)
seq_ch$region <- factor(seq_ch$region, levels = c("region_1","region_2","region_3","region_4","region_5","region_6", "Unknown"))
BAG_data$region  <- sapply(BAG_data$geoRegion, from_cantons_to_regions)
BAG_data$region <- factor(BAG_data$region, levels = c("region_1","region_2","region_3","region_4","region_5","region_6"))
#BAG_data$tests_num[BAG_data$geoRegion=="CH"]

region_names <- function(x){
  if(x== "0"){return("CH")}
  else if(x== "1"){return("region_1")}
  else if(x== "2"){return("region_2")}
  else if(x =="3"){return("region_3")}
  else if(x =="4"){return("region_4")}
  else if(x =="5"){return("region_5")}
  else if(x =="6"){return("region_6")}
}
#variants_ch$region <- sapply(variants_ch$region, region_names)
#variants_ch$region <- factor(variants_ch$region, levels = c("region_1","region_2","region_3","region_4","region_5","region_6","CH"))

setwd("./plots/")
dir_create <- paste0("",format(as_date(mean(period_date)), format = "%Y-%m"),"")
if (!dir.exists(dir_create)){
  dir.create(dir_create)
} 
setwd(paste0("./",format(as_date(mean(period_date)), format = "%Y-%m")))
dir_create <- paste0("pdf")
if (!dir.exists(dir_create)){
  dir.create(dir_create)
} 
dir_create <- paste0("png")
if (!dir.exists(dir_create)){
  dir.create(dir_create)
} 

setwd("../")
setwd("../")
setwd("./tables/")
dir_create <- paste0("./",format(as_date(mean(period_date)), format = "%Y-%m"))
if (!dir.exists(dir_create)){
  dir.create(dir_create)
} 
setwd("../")




