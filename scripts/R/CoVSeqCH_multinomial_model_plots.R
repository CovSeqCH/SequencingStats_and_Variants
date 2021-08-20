
#### Combating SARS-CoV-2 variants in Switzerland
## Created by M Reichmuth on June 2021

## help from here: https://github.com/tomwenseleers/newcovid_belgium/blob/main/functions/simulate_overdispersion.R

### set working directory
setwd("/Users/mr19m223/Documents/COVID_projects/SequencingStats_and_Variants")#may need to be adapted
### load data:
#source("./scripts/R/CoVSeqCH_data.R")

seqch <- subset(seq_ch, as_date(date) %in% seq(time_window[1],period_date[2],1))
### binomial confidence intervals, will be showed by week
#variants_ch$date <- as_date(variants_ch$date)
#variants_ch$date_num <- as.numeric(variants_ch$date)
#variants_ch$week <- week(variants_ch$date)
seqch$date <- as_date(seqch$date)
seqch$date_num <- as.numeric(seqch$date)
seqch$week <- week(seqch$date)

variant_week <- seqch %>% group_by(week,region,who_variants)  %>% summarise(variants= length(who_variants))
seq_week <- variant_week %>% group_by(week,region)  %>% summarise(total= sum(variants))
#variant_week <- aggregate(value~region+week+who_variants, variants_ch, sum)
#seq_week <- aggregate(sequences~region+week+who_variants, variants_ch, sum)
#variant_week<- merge(variant_week, seq_week,by=c("region","week","who_variants"))
variant_week <- merge(variant_week, seq_week,by=c("region","week"))
variant_week_ch <- seqch %>% group_by(week,who_variants)  %>% summarise(variants= length(who_variants))
seq_week_ch <- variant_week_ch %>% group_by(week)  %>% summarise(total= sum(variants))
variant_week_ch <- merge(variant_week_ch, seq_week_ch,by=c("week"))
variant_week_ch$region <- "CH"
variant_week <- rbind(variant_week_ch, variant_week)
remove(seq_week_ch)
remove(variant_week_ch)
remove(seq_week)

variant_week$total <- as.numeric(variant_week$total)
variant_week$variants <- as.numeric(variant_week$variants)
variant_week <- variant_week[variant_week$total!=0,]
#variant_week$sequences <- as.numeric(variant_week$sequences)
#variant_week$value <- as.numeric(variant_week$value)
#variant_week <- variant_week[variant_week$sequences!=0,]
variant_week$who_variants <- factor(variant_week$who_variants, levels = lev)

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
variant_week$week_day <- as_date(as.Date(paste(2021, variant_week$week, 1, sep="-"),format= "%Y-%U-%u"))
variant_week$region <- gsub("_"," ", variant_week$region)
variant_week$region <- gsub("r","R", variant_week$region)
variant_week$region <- factor(variant_week$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6", "CH"))
variant_week$week_day <- as_date(as.Date(paste(2021, variant_week$week, 1, sep="-"), "%Y-%U-%u"))
remove(conf)
remove(upper)
remove(lower)


### prepare data for multinomial function and plot
## reorganize dataset
#variants_reg <- variants_ch[variants_ch$region!="CH",]
#variants_ch <- variants_ch[variants_ch$region=="CH",]
variants_reg <- seqch[seqch$region!="Unknown",]
variants_ch <- seqch


## overall Swiss data
#variants_ch$region <- factor(variants_ch$region, levels = c("CH"))
#variants_ch <- variants_ch[variants_ch$value>0,]
#variants_ch <- variants_ch[rep(row.names(variants_ch), variants_ch$value), c(2,6,7)]
#variants_ch$who_variants <- factor(variants_ch$who_variants, levels = lev)

## regional data
variants_reg <- variants_reg[grepl("egion",variants_reg$region),]
variants_reg$region <- gsub("_"," ", variants_reg$region)
variants_reg$region <- gsub("r","R", variants_reg$region)
variants_reg$region <- factor(variants_reg$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6"))
#variants_reg <- variants_reg[variants_reg$value>0,]
#variants_reg <- variants_reg[rep(row.names(variants_reg), variants_reg$value), c(2,6,7)]
#variants_reg$who_variants <- factor(variants_reg$who_variants, levels = lev)
#variants_reg$region <- factor(variants_reg$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6"))


### Different multinomial models:
mnom_date_spline <- multinom(who_variants ~ ns(date_num, df=2), data = variants_ch)
mnom_date_reg_spline <- multinom(who_variants ~ ns(date_num, df=2) + region, data = variants_reg)

## predict for mnom_week model
predict_date <- data.frame(date_num = seq((min(variants_ch$date_num)),as.numeric(time_window[2]),1))

predict.eff_date <- Effect("date_num", mnom_date_spline,level=0.95,
                           xlevels=list(date_num=seq(min(variants_ch$date_num),as.numeric(time_window[2]),1)))
predict.eff_date <- cbind(predict_date, melt(predict.eff_date$prob), melt(predict.eff_date$lower.prob),melt(predict.eff_date$upper.prob))
predict.eff_date <- predict.eff_date[,-c(2,5,6,8,9)]
colnames(predict.eff_date) <- c("date_num","who_variants", "prob","lower", "upper")
predict.eff_date$who_variants <- gsub("prob.", "", predict.eff_date$who_variants)
predict.eff_date$who_variants <- sapply(predict.eff_date$who_variants, who_variant_names)
#table(predict.eff_date$who_variants)
predict.eff_date$who_variants <- factor(predict.eff_date$who_variants, levels = lev)
predict.eff_date$date <- as_date(predict.eff_date$date_num)
remove(predict_date)


## predict for mnom_week_reg model
predict_date_reg <- data.frame(date_num=rep(unique(seq(min(variants_reg$date_num),as.numeric(time_window[2]),1)), length(unique(variants_reg$region))), 
                               region=rep(levels(variants_reg$region),each=length(unique(seq((min(variants_reg$date_num)),as.numeric(time_window[2]),1)))))

predict.eff_date_reg <- Effect(c("date_num","region"), mnom_date_reg_spline,level=0.95,
                               xlevels=list(date_num=seq(min(variants_reg$date_num),as.numeric(time_window[2]),1)))
predict.eff_date_reg <- cbind(predict_date_reg[,c(1,2)],melt(predict.eff_date_reg$prob), melt(predict.eff_date_reg$lower.prob),melt(predict.eff_date_reg$upper.prob))
predict.eff_date_reg <- predict.eff_date_reg[,-c(3,6,7,9,10)]
colnames(predict.eff_date_reg) <- c("date_num","region", "who_variants", "prob","lower", "upper")
predict.eff_date_reg$who_variants <- gsub("prob.", "", predict.eff_date_reg$who_variants)
predict.eff_date_reg$who_variants <- sapply(predict.eff_date_reg$who_variants, who_variant_names)
#table(predict.eff_date_reg$who_variants)
predict.eff_date_reg$who_variants <- factor(predict.eff_date_reg$who_variants, levels = lev)
predict.eff_date_reg$region <- factor(predict.eff_date_reg$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6"))
predict.eff_date_reg$date <- as_date(predict.eff_date_reg$date_num)
remove(predict_date_reg)

remove(mnom_date_reg_spline)
remove(mnom_date_spline)

### Figures
## prepare ploting:
col_9 <- (brewer.pal(9,"Set1"))

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
variants_plot_model <- ggplot() + 
  geom_line(data=predict.eff_date, aes(x = date, y = prob, color = who_variants))+
  geom_ribbon(data=predict.eff_date, aes(x = date, y = prob, ymin = lower,ymax = upper,fill=who_variants),alpha=0.4)+
  geom_errorbar(data= variant_week[na.omit(variant_week$region=="CH"),], aes(x = week_day, ymin=lower, ymax=upper, color = who_variants), width=.1) +
  geom_point(data= variant_week[na.omit(variant_week$region=="CH"),],aes(x = week_day, y=conf, color = who_variants))+
  geom_rect(aes(xmin = as_date(max(seqch$date)), ymin = 0, xmax = time_window[2], ymax = 1), fill= col_9[9], colour= "transparent", alpha=0.4)+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window[1],time_window[2]+1)))+
  scale_color_manual(values= col_9[2:9],name="SARS-CoV-2 variants") +
  scale_fill_manual(values= col_9[2:9],name="SARS-CoV-2 variants") +
  theme_minimal()+
  theme(plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))+
  scale_y_continuous(limits = c(0,1))+
  coord_cartesian(expand = FALSE)+
  labs(tag=bquote(.("")),subtitle = "Switzerland", x = "", y =bquote("Proportion of SARS-CoV-2 variants"))
variants_plot_model_log <- variants_plot_model+
  coord_cartesian(ylim = c(10^-3, 10^0))+
  theme(legend.position = "none")+
  scale_y_continuous(trans='log10',labels=yscaling,limits = c(min(predict.eff_date$lower),1))+
  labs(y =bquote("Proportion of SARS-CoV-2 variants (log10-scale)"))
variants_legend <- g_legend(variants_plot_model,1)
variants_plot_model <- variants_plot_model + theme(legend.position = "none")

## reginal Swiss multi-nominal figure
regional_variants_plot_model <- ggplot(predict.eff_date_reg) + 
  geom_line(data= predict.eff_date_reg, aes(x = date, y = prob, color = who_variants))+
  geom_ribbon(data= predict.eff_date_reg, aes(x= date, ymin = lower,ymax = upper,fill=who_variants),alpha=0.4)+
  geom_errorbar(data= variant_week[na.omit(variant_week$region!="CH"),], aes(x = week_day, ymin=lower, ymax=upper, color = who_variants), width=.1) +
  geom_point(data= variant_week[na.omit(variant_week$region!="CH"),],aes(x = week_day, y=conf, color = who_variants))+#, size = conf
  geom_rect(aes(xmin = as_date(max(seqch$date)), ymin = 0, xmax = time_window[2], ymax = 1), fill= col_9[9], colour= "transparent", alpha=0.008)+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window[1],time_window[2]+15)))+
  scale_color_manual(values= col_9[2:9],name="SARS-CoV-2 variants") +
  scale_fill_manual(values= col_9[2:9],name="SARS-CoV-2 variants") +
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
  scale_y_continuous(trans='log10',labels=yscaling, limits = c(min(predict.eff_date$lower),1))+
  labs(y =bquote("Proportion of SARS-CoV-2 variants (log10-scale)"))
com_regional_ch_variants<- grid.arrange(variants_plot_model, variants_plot_model_log, 
                                        variants_legend,regional_variants_plot_model, 
                                        regional_variants_plot_model_log, nrow = 2)
ggsave(com_regional_ch_variants, filename = paste0("./plots/",format(period_date[2],"%Y-%m"),"/png/multinomial_variants_ch_reg_",format(period_date[2],"%Y-%m"), ".png"), height = 9, width = 18,  bg = "transparent")
ggsave(com_regional_ch_variants, filename = paste0("./plots/",format(period_date[2],"%Y-%m"),"/pdf/multinomial_variants_ch_reg_",format(period_date[2],"%Y-%m"), ".pdf"), height = 9, width = 18,  bg = "transparent")

