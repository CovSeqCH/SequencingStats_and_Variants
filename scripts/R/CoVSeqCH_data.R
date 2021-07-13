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
url <- "https://www.bag.admin.ch/dam/bag/de/dokumente/mt/k-und-i/aktuelle-ausbrueche-pandemien/2019-nCoV/covid-19-basisdaten-bevoelkerungszahlen.xlsx.download.xlsx/Population_Size_BFS.xlsx"
GET(url, write_disk(tf <- tempfile(fileext = ".xls")))
CH_canton_age_population <- read.xlsx(tf,sheetIndex = 2, startRow = 0)
colnames(CH_canton_age_population) <- c("canton", "gender", "age", "no_pop")

## Swiss SARS-CoV-2  metadata:
#html <- readLines("https://www.covid19.admin.ch/api/context")
#html <- gsub("^(.*)https:", "https:", html[18])
#url <- gsub("\"", "",html)
#download(url, dest=paste0("covid19_bag_", Sys.Date(),".zip"), mode="wb") 
#unzip(paste0("covid19_bag_", Sys.Date(),".zip"), exdir = "./")
#BAG_Re_canton <- read.csv("./COVID19Re_geoRegion.csv")
#BAG_test_canton <- read.csv("./COVID19Test_geoRegion_PCR_Antigen.csv")
#BAG_cases_canton <- read.csv("./COVID19Cases_geoRegion.csv")

## Swiss SARS-CoV-2 sequencing metadata:
variants_ch <- read.csv("https://raw.githubusercontent.com/CovSeqCH/SequencingStats_and_Variants/cr-dev/data/cases_seq_by_day_region.csv")

### prepare data / data cleaning:
month_end <- as.numeric(format(Sys.Date(),"%m"))
time_window <- c(as_date("2021-01-01"), as_date(paste0("2021-0", month_end+1,"-01"))-1)
month_start <- format(as.numeric(format(Sys.Date(),"%m"))-1, format="%m")
month_end <- format(month_end, format="%m")
period <- function(x) {
  c(as_date(format(x, paste0("%Y-",month_start,"-01"))),
    lubridate::floor_date(as_date(format(x, paste0("%Y-",month_end,"-01"))), unit = "month") - 1)
}
period_date <- period(Sys.Date())
period_days <- seq(period_date[1], period_date[2],1)


variants_ch <- subset(variants_ch, as_date(date) %in% seq(time_window[1],time_window[2],1))#Sys.Date()-14

#CH_canton_age_population <- subset(CH_canton_age_population, !canton %in% c("FL"))
#BAG_Re_canton <- subset(BAG_Re_canton, !geoRegion %in% c("FL", "CHFL") & as_date(date) %in% period_days)
#BAG_test_canton <- subset(BAG_test_canton, !geoRegion %in% c("FL", "CHFL") & as_date(datum) %in% period_days)
#BAG_cases_canton <- subset(BAG_cases_canton,!geoRegion %in% c("FL", "CHFL") &  as_date(datum) %in% period_days)
#metadata_ch <- variants_ch[names(variants_ch)  %in% c("region", "date","sequences", "cases")]

### renaming functions / variables generation
## Renaming variants according to WHO
#variants_data <- variants_ch[!names(variants_ch)  %in% c("sequences", "cases")]
variants_ch <- melt(variants_ch, id.vars=c("date", "region","cases","sequences"))
#variants_data <- melt(variants_data, id.vars=c("date", "region"))
who_variant_names <- function(x){
  if(x== "alpha"){return("Alpha")}
  else if(x== "beta"){return("Beta")}#else if(grepl("B.1.351",x)){return("beta")}
  else if(x== "gamma"){return("Gamma")}
  else if(x =="delta"){return("Delta")}
  else if(x =="C.36."){return("C.36*")}
  else if(x =="others"){return("others")}
  else{return(x)} 
}
variants_ch$who_variants <- sapply(variants_ch$variable, who_variant_names)
lev <- c("Alpha",  "Beta",  "Gamma", "Delta",  "C.36*",  "others")
variants_ch$who_variants <- factor(variants_ch$who_variants, levels = lev)
variants_ch$variable <- NULL

## Divide Switzerland in 6 regions as following:
## region 1 «GE, NE, VD, VS»
## region 2 «BE, FR, JU»
## region 3 «AG, BL, BS, SO»
## region 4 «LU, NW, OW, SZ, UR, ZG»
## region 5 «AI, AR, GL, SG, SH, TG, ZH»
## region 6 «GR, TI».
#CH_canton_age_population$regions <- sapply(CH_canton_age_population$canton, from_cantons_to_regions)
#BAG_Re_canton$regions  <- sapply(BAG_Re_canton$geoRegion, from_cantons_to_regions)
#BAG_test_canton$regions  <- sapply(BAG_test_canton$geoRegion, from_cantons_to_regions)
#BAG_cases_canton$regions  <- sapply(BAG_cases_canton$geoRegion, from_cantons_to_regions)

region_names <- function(x){
  if(x== "0"){return("CH")}
  else if(x== "1"){return("region_1")}
  else if(x== "2"){return("region_2")}
  else if(x =="3"){return("region_3")}
  else if(x =="4"){return("region_4")}
  else if(x =="5"){return("region_5")}
  else if(x =="6"){return("region_6")}
}
variants_ch$region <- sapply(variants_ch$region, region_names)
variants_ch$region <- factor(variants_ch$region, levels = c("region_1","region_2","region_3","region_4","region_5","region_6","CH"))
#metadata_ch$region <- sapply(metadata_ch$region, region_names)
#metadata_ch$region <- factor(metadata_ch$region, levels = c("region_1","region_2","region_3","region_4","region_5","region_6","CH"))


