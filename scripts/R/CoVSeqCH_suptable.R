
##### R file created on 21 May 2021 by Martina Reichmuth
#### National sequencing surveillance program on SARS-CoV-2 variants

### set working directory
setwd("/Users/mr19m223/Documents/COVID_projects/SequencingStats_and_Variants")#may need to be adapted
### load data:
source("./scripts/R/CoVSeqCH_data.R")
#period_date <- c(as_date("2021-11-01"), as_date("2021-12-05"))
BAGdata <- subset(BAG_data, as_date(date) %in% seq(period_date[1],period_date[2],1))
seqch <- subset(seq_ch, as_date(date) %in% seq(period_date[1],period_date[2],1))
#seqch1 <- subset(seq_ch1, as_date(date) %in% seq(period_date[1],period_date[2],1))

#BAGdata <- subset(BAG_data, as_date(date) %in% seq(as_date("2021-11-29"), as_date("2021-12-14"),1))
#seqch <- subset(seq_ch, as_date(date) %in% seq(as_date("2021-11-29"), as_date("2021-12-14"),1))

seqch_undetermined <- seqch[seqch$who_variants %in% "undetermined",]
seqch <- seqch[!seqch$who_variants %in% "undetermined",]

#test_seq  <- subset(seq_ch, as_date(date) %in% seq(as_date("2021-08-30"),as_date("2021-10-03"),1))


### Supplementary table creation:
table <- as.data.frame(matrix(ncol=46, nrow= 33))
rownames(table) <- c("CH",
                     "region_1",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_1"])),
                     "region_2",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_2"])),
                     "region_3",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_3"])),
                     "region_4",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_4"])),
                     "region_5",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_5"])),
                     "region_6",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_6"]))
)

for (c in unname(unlist(rownames(table)))){
  #population size:
  if(grepl("region", c)){
    table[,1][rownames(table) == c] <- sum(unique(BAGdata$pop[BAGdata$region %in% c]))
  }
  else{
    table[,1][rownames(table) == c] <- unique(BAGdata$pop[BAGdata$geoRegion %in% c])
  }
  #confirmed cases:
  if(grepl("region", c)){
    table[,2][rownames(table) == c] <- sum(na.omit(BAGdata$cases_num[BAGdata$region %in% c]))
  }
  else{
    table[,2][rownames(table) == c] <- sum(na.omit(BAGdata$cases_num[BAGdata$geoRegion %in% c]))
  }
  #14d-incidence
  period_14 <- length(period_days)/14
  if(grepl("region", c)){
    table[,3][rownames(table) == c] <- round(sum(na.omit(BAGdata$cases_num[BAGdata$region %in% c]))/sum(unique(BAGdata$pop[BAGdata$region %in% c]))*period_14*10^5)
  }
  else{
    table[,3][rownames(table) == c] <- round(sum(na.omit(BAGdata$cases_num[BAGdata$geoRegion %in% c]))/unique(BAGdata$pop[BAGdata$geoRegion %in% c])*period_14*10^5)
  }
  #effective reproduction number taken from BAG 
  if(grepl("region", c)){
    table[,4][rownames(table) == c] <- paste0(format(round(median(na.omit(BAGdata$median_R_mean[BAGdata$region %in% c])),2), nsmall = 2)," (",format(round(median(na.omit(BAGdata$median_R_lowHPD[BAGdata$region %in% c])),2), nsmall = 2),"-",format(round(median(na.omit(BAGdata$median_R_highHPD[BAGdata$region %in% c])),2), nsmall = 2),")")
  }
  else{
    table[,4][rownames(table) == c] <- paste0(format(round(median(na.omit(BAGdata$median_R_mean[BAGdata$geoRegion %in% c])),2), nsmall = 2)," (",format(round(median(na.omit(BAGdata$median_R_lowHPD[BAGdata$geoRegion %in% c])),2), nsmall = 2),"-",format(round(median(na.omit(BAGdata$median_R_highHPD[BAGdata$geoRegion %in% c])),2), nsmall = 2),")")
  }
  # number of tests
  if(grepl("region", c)){
    table[,5][rownames(table) == c] <- sum(na.omit(BAGdata$tests_num[BAGdata$region %in% c]))
  }
  else{
    table[,5][rownames(table) == c] <- sum(na.omit(BAGdata$tests_num[BAGdata$geoRegion %in% c]))
  }
  # incidence of tests for 14d over period
  if(grepl("region", c)){
    table[,6][rownames(table) == c] <- round(sum(na.omit(BAGdata$tests_num[BAGdata$region %in% c]))/sum(unique(BAGdata$pop[BAGdata$region %in% c]))*10^5)
  }
  else{
    table[,6][rownames(table) == c] <- round(sum(na.omit(BAGdata$tests_num[BAGdata$geoRegion %in% c]))/unique(BAGdata$pop[BAGdata$geoRegion %in% c])*10^5)
  }
  
  # test positivity
  if(grepl("region", c)){
    table[,7][rownames(table) == c] <- format(round(sum(na.omit(BAGdata$tests_pos_num[BAGdata$region %in% c]))/sum((na.omit(BAGdata$tests_num[BAGdata$region %in% c])))*10^2,2), nsmall = 2)
  }
  else{
    table[,7][rownames(table) == c] <- format(round(sum(na.omit(BAGdata$tests_pos_num[BAGdata$geoRegion %in% c]))/sum(na.omit(BAGdata$tests_num[BAGdata$geoRegion %in% c]))*10^2,2), nsmall = 2)
  }
  
  # number of sequenced samples 
  if(grepl("CH", c)){
    table[,8][rownames(table) == c] <- length(seqch$who_variants[seqch$country %in% c])
  }
  else if(grepl("region", c)){
    table[,8][rownames(table) == c] <- length(seqch$who_variants[seqch$region %in% c])
  }
  else{
    table[,8][rownames(table) == c] <- length(seqch$who_variants[seqch$canton %in% c])
  }
  
  # Proportion of cases that have been sequenced
  if(grepl("CH", c)){
    table[,9][rownames(table) == c] <- format(round(length(na.omit(seqch$who_variants[seqch$country %in% c]))/sum(BAGdata$cases_num[BAGdata$geoRegion %in% c])*100,1), nsmall = 1)
  }
  else if(grepl("region", c)){
    table[,9][rownames(table) == c] <- format(round(length(na.omit(seqch$who_variants[seqch$region %in% c]))/sum(BAGdata$cases_num[BAGdata$region %in% c])*100,1), nsmall = 1)
  }
  else{
    table[,9][rownames(table) == c] <- format(round(length(na.omit(seqch$who_variants[seqch$canton %in% c]))/sum(BAGdata$cases_num[BAGdata$geoRegion %in% c])*100,1), nsmall = 1)
  }
  
  
  # number of  "alpha" 
  if(grepl("CH", c)){
    table[,10][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Alpha"]))
  }
  else if(grepl("region", c)){
    table[,10][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Alpha"]))
  }
  else{
    table[,10][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Alpha"]))
  }
  
  # Percentage of  "alpha" to all sequences
  if(table[,10][rownames(table) == c]!=0 & !is.na(table[,10][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,10][rownames(table) == c]), table[,8][rownames(table) == c])  # binomial 95% confidence intervals
    table[,11][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,11][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,11][rownames(table) == c] <- "-"
  }
  
  # Number of  beta 
  if(grepl("CH", c)){
    table[,12][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Beta"]))
  }
  else if(grepl("region", c)){
    table[,12][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Beta"]))
  }
  else{
    table[,12][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Beta"]))
  }
  
  # Proportion of  beta to all sequences
  if(table[,12][rownames(table) == c]!=0 & !is.na(table[,12][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,12][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,13][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,13][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,13][rownames(table) == c] <- "-"
  }
  
  # Number of gamma 
  if(grepl("CH", c)){
    table[,14][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Gamma"]))
  }
  else if(grepl("region", c)){
    table[,14][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Gamma"]))
  }
  else{
    table[,14][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Gamma"]))
  }
  
  # Proportion of  Gamma to all sequences
  if(table[,14][rownames(table) == c]!=0 & !is.na(table[,14][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,14][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,15][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,15][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,15][rownames(table) == c] <- "-"
  }
  
  # Number of delta 
  if(grepl("CH", c)){
    table[,16][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Delta"]))
  }
  else if(grepl("region", c)){
    table[,16][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Delta"]))
  }
  else{
    table[,16][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Delta"]))
  }

  # Proportion of  Delta to all sequences
  if(table[,16][rownames(table) == c]!=0 & !is.na(table[,16][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,16][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,17][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,17][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,17][rownames(table) == c] <- "-"
  }
  
  # Number of Lambda 
  if(grepl("CH", c)){
    table[,18][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Lambda"]))
  }
  else if(grepl("region", c)){
    table[,18][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Lambda"]))
  }
  else{
    table[,18][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Lambda"]))
  }
  
  # Proportion of  Lambda to all sequences
  if(table[,18][rownames(table) == c]!=0 & !is.na(table[,18][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,18][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,19][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,19][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,19][rownames(table) == c] <- "-"
  }
  
  # Number of B.1.1.318 
  if(grepl("CH", c)){
    table[,20][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "B.1.1.318"]))
  }
  else if(grepl("region", c)){
    table[,20][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "B.1.1.318"]))
  }
  else{
    table[,20][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "B.1.1.318"]))
  }
  
  # Proportion of  B.1.1.318 to all sequences
  if(table[,20][rownames(table) == c]!=0 & !is.na(table[,20][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,20][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,21][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,21][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,21][rownames(table) == c] <- "-"
  }
  
  # Number of Mu 
  if(grepl("CH", c)){
    table[,22][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Mu"]))
  }
  else if(grepl("region", c)){
    table[,22][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Mu"]))
  }
  else{
    table[,22][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Mu"]))
  }
  
  # Proportion of  Mu to all sequences
  if(table[,22][rownames(table) == c]!=0 & !is.na(table[,22][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,22][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,23][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,23][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,23][rownames(table) == c] <- "-"
  }
  # Number of Omicron BA.1
  if(grepl("CH", c)){
    table[,24][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (BA.1)"]))
  }
  else if(grepl("region", c)){
    table[,24][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (BA.1)"]))
  }
  else{
    table[,24][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (BA.1)"]))
  }
  
  # Proportion of  Omicron to all sequences
  if(table[,24][rownames(table) == c]!=0 & !is.na(table[,24][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,24][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,25][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,25][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,25][rownames(table) == c] <- "-"
  }
  # Number of Omicron BA.2
  if(grepl("CH", c)){
    table[,26][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (BA.2)"]))
  }
  else if(grepl("region", c)){
    table[,26][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (BA.2)"]))
  }
  else{
    table[,26][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (BA.2)"]))
  }
  
  # Proportion of  Omicron to all sequences
  if(table[,26][rownames(table) == c]!=0 & !is.na(table[,26][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,26][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,27][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,27][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,27][rownames(table) == c] <- "-"
  }
  # Number of Omicron BA.3
  if(grepl("CH", c)){
    table[,28][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (BA.3)"]))
  }
  else if(grepl("region", c)){
    table[,28][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (BA.3)"]))
  }
  else{
    table[,28][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (BA.3)"]))
  }
  # Proportion of  Omicron to all sequences
  if(table[,28][rownames(table) == c]!=0 & !is.na(table[,28][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,28][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,29][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,29][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,29][rownames(table) == c] <- "-"
  }
  # Number of Omicron BA.4
  if(grepl("CH", c)){
    table[,30][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (BA.4)"]))
  }
  else if(grepl("region", c)){
    table[,30][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (BA.4)"]))
  }
  else{
    table[,30][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (BA.4)"]))
  }
  # Proportion of  Omicron to all sequences
  if(table[,30][rownames(table) == c]!=0 & !is.na(table[,30][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,30][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,31][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,31][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,31][rownames(table) == c] <- "-"
  }
  # Number of Omicron BA.5
  if(grepl("CH", c)){
    table[,32][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (BA.5)"]))
  }
  else if(grepl("region", c)){
    table[,32][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (BA.5)"]))
  }
  else{
    table[,32][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (BA.5)"]))
  }
  # Proportion of  Omicron to all sequences
  if(table[,32][rownames(table) == c]!=0 & !is.na(table[,32][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,32][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,33][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,33][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,33][rownames(table) == c] <- "-"
  }
  
  # Number of Omicron BA.2.12.1
  if(grepl("CH", c)){
    table[,34][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (BA.2.12.1)"]))
  }
  else if(grepl("region", c)){
    table[,34][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (BA.2.12.1)"]))
  }
  else{
    table[,34][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (BA.2.12.1)"]))
  }
  # Proportion of  Omicron to all sequences
  if(table[,34][rownames(table) == c]!=0 & !is.na(table[,34][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,34][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,35][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,35][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,35][rownames(table) == c] <- "-"
  }
  
  
  
  # Number of Omicron BA.1 & BA.2
  if(grepl("CH", c)){
    table[,36][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (BA.1 & BA.2)"]))
  }
  else if(grepl("region", c)){
    table[,36][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (BA.1 & BA.2)"]))
  }
  else{
    table[,36][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (BA.1 & BA.2)"]))
  }
  # Proportion of  Omicron to all sequences
  if(table[,36][rownames(table) == c]!=0 & !is.na(table[,36][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,36][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,37][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,37][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,37][rownames(table) == c] <- "-"
  }
  
  # Number of Omicron BA.2.75
  if(grepl("CH", c)){
    table[,38][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (BA.2.75)"]))
  }
  else if(grepl("region", c)){
    table[,38][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (BA.2.75)"]))
  }
  else{
    table[,38][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (BA.2.75)"]))
  }
  # Proportion of  Omicron to all sequences
  if(table[,38][rownames(table) == c]!=0 & !is.na(table[,38][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,38][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,39][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,39][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,39][rownames(table) == c] <- "-"
  }
  # Number of Omicron BQ.1
  if(grepl("CH", c)){
    table[,40][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (BQ.1)"]))
  }
  else if(grepl("region", c)){
    table[,40][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (BQ.1)"]))
  }
  else{
    table[,40][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (BQ.1)"]))
  }
  # Proportion of  Omicron to all sequences
  if(table[,40][rownames(table) == c]!=0 & !is.na(table[,40][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,40][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,41][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,41][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,41][rownames(table) == c] <- "-"
  }
  # Number of Omicron XBB
  if(grepl("CH", c)){
    table[,42][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "Omicron (XBB)"]))
  }
  else if(grepl("region", c)){
    table[,42][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "Omicron (XBB)"]))
  }
  else{
    table[,42][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "Omicron (XBB)"]))
  }
  # Proportion of  Omicron to all sequences
  if(table[,42][rownames(table) == c]!=0 & !is.na(table[,42][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,42][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,43][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,43][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,43][rownames(table) == c] <- "-"
  }
  
  
  # Number of "other"
  if(grepl("CH", c)){
    table[,44][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$country %in% c & seqch$who_variants  %in% "others"]))
  }
  else if(grepl("region", c)){
    table[,44][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$region %in% c & seqch$who_variants  %in% "others"]))
  }
  else{
    table[,44][rownames(table) == c] <- length(na.omit(seqch$who_variants[seqch$canton %in% c & seqch$who_variants  %in% "others"]))
  }
  
  # Proportion of  others to all sequences
  if(table[,44][rownames(table) == c]!=0 & !is.na(table[,44][rownames(table) == c])){
    interval <- binom.test(as.numeric(table[,44][rownames(table) == c]), table[,8][rownames(table) == c])  #binomial 95% confidence intervals
    table[,45][rownames(table) == c] <- paste0(format(round(interval$estimate*100,1), nsmall = 1), " (",format(round(interval$conf.int[1]*100,1), nsmall = 1),"-",format(round(interval$conf.int[2]*100,1), nsmall = 1),")")
    if(interval$estimate<"0.001" & interval$estimate!="0"){
      table[,45][rownames(table) == c] <- "<0.1"
    }
  }
  else{
    table[,45][rownames(table) == c] <- "-"
  }
  # Number of "undetermined"
  if(grepl("CH", c)){
    table[,46][rownames(table) == c] <- length(na.omit(seqch_undetermined$who_variants[seqch_undetermined$country %in% c & seqch_undetermined$who_variants  %in% "undetermined"]))
  }
  else if(grepl("region", c)){
    table[,46][rownames(table) == c] <- length(na.omit(seqch_undetermined$who_variants[seqch_undetermined$region %in% c & seqch_undetermined$who_variants  %in% "undetermined"]))
  }
  else{
    table[,46][rownames(table) == c] <- length(na.omit(seqch_undetermined$who_variants[seqch_undetermined$canton %in% c & seqch_undetermined$who_variants  %in% "undetermined"]))
  }
  
}

colnames(table) <- c("Population size",
                     "Confirmed cases", "14-day incidence of confirmed cases (per 100,000)",
                     "Effective reproduction number (median, range)", "Number of tests (PCR and antigen)",
                     "Incidence of tests (per 100,000)","Test positivity (%)",
                     "Sequenced samples","Proportion sequenced (%)",
                     "Alpha","Percentage Alpha (95% CI)",
                     "Beta","Percentage Beta (95% CI)",
                     "Gamma","Percentage Gamma (95% CI)",
                     "Delta","Percentage Delta (95% CI)",
                     "Lambda","Percentage Lambda (95% CI)",
                     "B.1.1.318","Percentage B.1.1.318 (95% CI)",
                     "Mu","Percentage Mu (95% CI)",
                     "Omicron (BA.1)","Percentage Omicron (BA.1) (95% CI)",
                     "Omicron (BA.2)","Percentage Omicron (BA.2) (95% CI)",
                     "Omicron (BA.3)","Percentage Omicron (BA.3) (95% CI)",
                     "Omicron (BA.4)","Percentage Omicron (BA.4) (95% CI)",
                     "Omicron (BA.5)","Percentage Omicron (BA.5) (95% CI)",
                     "Omicron (BA.2.12.1)","Percentage Omicron (BA.2.12.1) (95% CI)",
                     "Omicron (BA.1 & BA.2)","Percentage Omicron (BA.1 & BA.2) (95% CI)",
                     "Omicron (BA.2.75)","Percentage Omicron (BA.2.75) (95% CI)",
                     "Omicron (BA.2.75)","Percentage Omicron (BQ.1) (95% CI)",
                     "Omicron (BA.2.75)","Percentage Omicron (XBB) (95% CI)",
                     "Other variants","Percentage other variants (95% CI)",
                     "Undetermined sequences, excluded from further analysis")
#sum(table$`Confirmed cases`[grepl("region",rownames(table))])==sum(table$`Confirmed cases`[grepl("CH",rownames(table))])
#sum(seqch$canton=="Unknown")+sum(table$`Sequenced samples`[grepl("region",rownames(table))])==sum(table$`Sequenced samples`[grepl("CH",rownames(table))])# 


rownames(table) <- c("Switzerland overall",
                     "Region 1",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_1"])),
                     "Region 2",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_2"])),
                     "Region 3",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_3"])),
                     "Region 4",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_4"])),
                     "Region 5",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_5"])),
                     "Region 6",
                     unique(na.omit(BAGdata$geoRegion[BAGdata$region=="region_6"]))
)

table[is.na(table)] <- "-"
# output table:
#period_date1 <- period_date
#period_date <- period(Sys.Date())
#period_days <- seq(period_date[1], period_date[2],1)
table <- table[-c(4),]

perioddate <- paste0(unique(format(period_date, format="%b")),collapse = "_")
write.xlsx(table, paste0("./tables/",format(as_date(mean(period_date)), format = "%Y-%m"),"/sup_table_overview_",perioddate,".xlsx"),sheetName=paste0("",period_date[1], " to ", period_date[2]), 
           col.names=TRUE, row.names=TRUE, append=FALSE)


table <- table[grepl("Region|Switzerland", rownames(table)),]
table<- table[,-c(1:7)]

write.xlsx(table, paste0("./tables/",format(as_date(mean(period_date)), format = "%Y-%m"),"/regional_table_",perioddate,".xlsx"),sheetName=paste0("",period_date[1], " to ", period_date[2]), 
           col.names=TRUE, row.names=TRUE, append=FALSE)


remove(table)

