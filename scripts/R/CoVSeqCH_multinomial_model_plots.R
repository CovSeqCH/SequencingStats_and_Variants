
#### Combating SARS-CoV-2 variants in Switzerland
## Created by M Reichmuth on June 2021

## help from here: https://github.com/tomwenseleers/newcovid_belgium/blob/main/functions/simulate_overdispersion.R

### set working directory
setwd("/Users/mr19m223/Documents/COVID_projects/SequencingStats_and_Variants")#may need to be adapted
### load data:
seqch <- subset(seq_ch, as_date(date) %in% seq(time_window[1],period_date[2],1))
seqch$who_variants[seqch$who_variants %in% "Mu"] <- "others"
seqch$who_variants[seqch$who_variants  %in% "B.1.1.318"] <- "others"
seqch <- seqch[!seqch$who_variants  %in% "undetermined",]

lev <- names(table(seqch$who_variants)[table(seqch$who_variants)>0])
seqch <- seqch[seqch$who_variants %in% lev,]

lev_variants <- c("Alpha",  "Beta",  "Gamma", "Delta",
                  "Omicron (BA.1)","Omicron (BA.2)",
                  "Omicron (BA.1 & BA.2)",
                  "Omicron (BA.3)","Omicron (BA.4)","Omicron (BA.5)",
                  "Omicron (BA.2.75)","Omicron (BQ.1)",
                  "Omicron (XBB)","Omicron (XBB.1.5)",
                  "Omicron (rest)",
                  "others", "undetermined")
color_variants <- c("#E29391", "#E5A993", "#E8C195","#EBD898",
                             "#D0A4CC", "#BFC9A3", "#BFA9F1",
                             "#ACDDC5", "#CE9CAF", "#B0D9A1",
                             "#ABDBDC","#ADCBA5", 
                             "#A2D0C2", "#A42416",
                             "#35358c",
                             "#5e5e5d", "#FFFFFF")
color_variants <- setNames(color_variants, lev_variants)
#variants included in data 
color_variants <- color_variants[names(color_variants) %in% lev]
lev_variants <- names(color_variants)
seqch$who_variants <- factor(seqch$who_variants, levels = lev_variants)

seqch$date <- as_date(seqch$date)
seqch$date_num <- as.numeric(seqch$date)

seqch$week <- paste0(year(seqch$date),"-",week(seqch$date))
seqch$week <- gsub("53","52", seqch$week)
seqch$week_day <- parse_date_time(paste(seqch$week, '-Sun'), "%Y-%W-%a")
seqch$week_day <- as_date(seqch$week_day)

variant_week <- seqch %>% group_by(week_day,region,who_variants)  %>% summarise(variants= length(who_variants))
seq_week <- variant_week %>% group_by(week_day,region)  %>% summarise(total= sum(variants))
variant_week <- merge(variant_week, seq_week,by=c("region","week_day"))
variant_week_ch <- seqch %>% group_by(week_day,who_variants)  %>% summarise(variants= length(who_variants))
seq_week_ch <- variant_week_ch %>% group_by(week_day)  %>% summarise(total= sum(variants))
variant_week_ch <- merge(variant_week_ch, seq_week_ch,by=c("week_day"))
variant_week_ch$region <- "CH"
variant_week <- rbind(variant_week_ch, variant_week)
remove(seq_week_ch)
remove(variant_week_ch)
remove(seq_week)

variant_week$total <- as.numeric(variant_week$total)
variant_week$variants <- as.numeric(variant_week$variants)
variant_week <- variant_week[variant_week$total!=0,]
variant_week$who_variants <- factor(variant_week$who_variants, levels = lev_variants)
#variant_week$week_day <-as.POSIXct( paste(variant_week$week,1,  sep = "-" ), format = "%Y-%U-%u" )
variant_week$week_day <- as_date(variant_week$week_day)

lower <- numeric()
upper <- numeric()
conf <- numeric()
sampled_ci <- function(data){
  for(i in 1:length(data$who_variants)) {
    int <- binom.test(data$variants[i], data$total[i])#int <- binom.test(data$value[i], data$sequences[i])
    lower[i] <- int$conf.int[1]
    upper[i] <- int$conf.int[2]
    conf[i] <- int$estimate
  }
  data <- cbind(data, conf=conf, lower = lower, upper = upper)
 return(data)
}
variant_week <- sampled_ci(variant_week)
variant_week <- as.data.frame(variant_week)
variant_week$region <- gsub("_"," ", variant_week$region)
variant_week$region <- gsub("r","R", variant_week$region)
variant_week$region <- factor(variant_week$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6", "CH"))
remove(conf)
remove(upper)
remove(lower)


### prepare data for multinomial function and plot
## reorganize dataset
variants_reg <- seqch[!seqch$region %in% c("Unknown"),]
variants_ch <- seqch


## overall Swiss data
variants_ch$who_variants <- factor(variants_ch$who_variants, levels = lev_variants)

## regional data
variants_reg <- variants_reg[grepl("egion",variants_reg$region),]
variants_reg$region <- gsub("_"," ", variants_reg$region)
variants_reg$region <- gsub("r","R", variants_reg$region)
variants_reg$region <- factor(variants_reg$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6"))

variants_reg$who_variants <- factor(variants_reg$who_variants, levels = lev_variants)
variants_reg$region <- factor(variants_reg$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6"))


min_time = min(as.numeric(variants_ch$date))
variants_ch$date_num <- as.numeric(as.numeric(variants_ch$date)-min_time)
variants_reg$date_num <- as.numeric(as.numeric(variants_reg$date)-min_time)
time_window_num <- as.numeric(seq(time_window[1],time_window[2], by=1))-min_time

variants_ch$who_variants <-factor(variants_ch$who_variants, levels = lev_variants)
variants_reg$who_variants <-factor(variants_reg$who_variants, levels = lev_variants)

### Different multinomial models:
mnom_date_spline <- multinom(who_variants ~ ns(date_num, df=2), data = variants_ch)
mnom_date_reg_spline <- multinom(who_variants ~ ns(date_num, df=2) + region, data = variants_reg)

## predict for mnom_week model
predict_mnom_date =  as.data.frame(emmeans(mnom_date_spline,~ who_variants,
                                           by="date_num",  mode="prob",
                                           at=list(date_num=time_window_num)))
predict_mnom_date <- predict_mnom_date[,c("date_num","who_variants", "prob","lower.CL", "upper.CL")]
predict_mnom_date$who_variants <- factor(predict_mnom_date$who_variants, levels = lev_variants)
predict_mnom_date$date <- as_date(predict_mnom_date$date_num+min_time)
colnames(predict_mnom_date) <- c("date_num","who_variants", "prob","lower", "upper", "date")

## predict for mnom_week_reg model
predict.eff_date_reg <- predict.eff_date_reg_backup <- effects::Effect(c("date_num","region"), mnom_date_reg_spline,level=0.95,
                                                                       xlevels=list(date_num=time_window_num))
predict_date_reg <- as.data.frame(predict.eff_date_reg$x)
predict_date_reg <- predict_date_reg[rep(seq_len(nrow(predict_date_reg)), times=length(lev_variants)),]
predict.eff_date_reg <- cbind(reshape2::melt(predict.eff_date_reg$prob), 
                              reshape2::melt(predict.eff_date_reg$lower.prob),
                              reshape2::melt(predict.eff_date_reg$upper.prob))
predict.eff_date_reg <- cbind(predict_date_reg,predict.eff_date_reg)

predict.eff_date_reg <- predict.eff_date_reg[,-c(3,6,7,9,10)]
colnames(predict.eff_date_reg) <- c("date_num","region", "who_variants", "prob","lower", "upper")
predict.eff_date_reg$who_variants <- gsub("prob.", "", predict.eff_date_reg$who_variants)
predict.eff_date_reg$who_variants <- gsub("Omicron..","Omicron (",predict.eff_date_reg$who_variants)
predict.eff_date_reg$who_variants <- gsub("\\.$",")",predict.eff_date_reg$who_variants)
predict.eff_date_reg$who_variants <- factor(predict.eff_date_reg$who_variants, levels = lev_variants)
predict.eff_date_reg$date <- as_date(predict.eff_date_reg$date_num+min_time)

### Figures
## prepare ploting:
yscaling <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "e", l)
  l <- gsub("\\+", "", l)
  l <- gsub("e", "10^", l)
  parse(text=l)
}
g_legend <- function(a.gplot,num){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]] 
  if(num==1){
    legend$grobs[[1]]$grobs[[1]] <-  editGrob(legend$grobs[[1]]$grobs[[1]], gp=gpar(fill ="transparent",col="transparent"))
  }
  if(num==2){
    legend$grobs[[1]]$grobs[[1]] <-  editGrob(legend$grobs[[1]]$grobs[[1]], gp=gpar(fill ="transparent",shape="transparent"))
    legend$grobs[[2]]$grobs[[1]] <-  editGrob(legend$grobs[[2]]$grobs[[1]], gp=gpar(fill ="transparent",shape="transparent"))
  }
  legend
} 

## overall Swiss multi-nominal figure
variant_week_ch <- variant_week[na.omit(variant_week$region %in%"CH"),]
variant_week_ch$who_variants <- factor(variant_week_ch$who_variants, levels = lev_variants)
variants_plot_model <- ggplot() + 
  geom_line(data=predict_mnom_date, aes(x = date, y = prob, color = who_variants))+#predict.eff_date_reg[predict.eff_date_reg$region=="Region 1",]
  geom_ribbon(data=predict_mnom_date, aes(x = date, y = prob, ymin = lower,ymax = upper,fill=who_variants),alpha=0.4)+
  geom_errorbar(data= variant_week[na.omit(variant_week$region%in%"CH"),], aes(x = week_day, ymin=lower, ymax=upper,color=who_variants), width=.1) +
  geom_point(data= variant_week[na.omit(variant_week$region%in%"CH"),],aes(x = week_day, y=conf, color = who_variants))+
  geom_rect(aes(xmin = as_date(max(seqch$date)), ymin = 0, xmax = time_window[2], ymax = 1), fill= "#e8e8e8", colour= "transparent", alpha=0.4)+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window[1],time_window[2]+1)))+
  scale_color_manual(values= color_variants,name="SARS-CoV-2 variants", lev_variants) +
  scale_fill_manual(values= color_variants,name="SARS-CoV-2 variants", lev_variants) +
  theme_minimal()+
  theme(plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))+
  scale_y_continuous(limits = c(0,1))+
  coord_cartesian(expand = FALSE)+
  labs(tag=bquote(.("")),subtitle = "Switzerland", x = "", y =bquote("Proportion of SARS-CoV-2 variants"))
#ggsave(variants_plot_model, filename = paste0("./plots/",format(period_date[2],"%Y-%m"),"/png/multinomial_variants_ch_",format(period_date[2],"%Y-%m"), ".png"), height = 3, width = 6,  bg = "transparent")

variants_plot_model_log <- variants_plot_model+
  coord_cartesian(ylim = c(10^-3, 10^0))+
  theme(legend.position = "none")+
  scale_y_continuous(trans='log10',labels=yscaling,limits = c(min(predict_mnom_date$lower),1))+
  labs(y =bquote("Proportion of SARS-CoV-2 variants (log10-scale)"))
variants_legend <- g_legend(variants_plot_model,1)
variants_plot_model <- variants_plot_model + theme(legend.position = "none")

## regional Swiss multi-nominal figure
regional_variants_plot_model <- NULL
regional_variants_plot_model <-  ggplot(predict.eff_date_reg) + 
  geom_line(data= predict.eff_date_reg, aes(x = date, y = prob, color = who_variants))+
  geom_ribbon(data= predict.eff_date_reg, aes(x= date, ymin = lower,ymax = upper,fill=who_variants),alpha=0.02)+
  geom_errorbar(data= variant_week[na.omit(variant_week$region!="CH"),], aes(x = week_day, ymin=lower, ymax=upper, color = who_variants), width=.1) +
  geom_point(data= variant_week[na.omit(variant_week$region!="CH"),],aes(x = week_day, y=conf, color = who_variants))+#, size = conf
  geom_rect(aes(xmin = as_date(max(seqch$date)), ymin = 0, xmax = time_window[2], ymax = 1), fill= "#e8e8e8", colour= "transparent", alpha=0.4)+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window[1],time_window[2]+15)))+
  scale_color_manual(values= color_variants,name="SARS-CoV-2 variants") +
  scale_fill_manual(values= color_variants,name="SARS-CoV-2 variants") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))+
  facet_wrap(~ region, ncol=2)+
  scale_y_continuous(limits = c(0,1))+
  coord_cartesian(expand = FALSE)+
  labs(tag=bquote(.("")), x = "", y =bquote("Proportion of SARS-CoV-2 variants"))
regional_variants_plot_model_log <- regional_variants_plot_model+
  coord_cartesian(ylim = c(10^-3, 1))+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window[1],time_window[2]+1)))+
  scale_y_continuous(trans='log10',labels=yscaling, limits = c(min(predict.eff_date_reg$lower),1))+
  labs(y =bquote("Proportion of SARS-CoV-2 variants (log10-scale)"))
com_regional_ch_variants<- grid.arrange(variants_plot_model, variants_plot_model_log, 
                                        variants_legend,regional_variants_plot_model, 
                                        regional_variants_plot_model_log, nrow = 2)
log_regional_ch_variants<- grid.arrange(variants_plot_model_log,
                                        regional_variants_plot_model_log,
                                        variants_legend, nrow = 1)

ggsave(com_regional_ch_variants, filename = paste0("./plots/",format(period_date[2],"%Y-%m"),"/png/multinomial_variants_ch_reg_",format(period_date[2],"%Y-%m"), ".png"), height = 9, width = 27,  bg = "transparent")
ggsave(com_regional_ch_variants, filename = paste0("./plots/",format(period_date[2],"%Y-%m"),"/pdf/multinomial_variants_ch_reg_",format(period_date[2],"%Y-%m"), ".pdf"), height = 9, width = 27,  bg = "transparent")
ggsave(log_regional_ch_variants, filename = paste0("./plots/",format(period_date[2],"%Y-%m"),"/png/multinomial_variants_ch_reg_log_",format(period_date[2],"%Y-%m"), ".png"), height = 6, width = 20,  bg = "transparent")

