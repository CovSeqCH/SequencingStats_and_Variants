
#### Combating SARS-CoV-2 variants in Switzerland
## Created by M Reichmuth on June 2021

## help from here: https://github.com/tomwenseleers/newcovid_belgium/blob/main/functions/simulate_overdispersion.R #https://zepel.io/blog/how-to-create-a-new-branch-in-github/#create-branch-command-line

### set working directory
setwd("/Users/mr19m223/Documents/COVID_projects/Swiss_sequencing_project")
### load data:
source("./R/CoVSeqCH_data.R")


### binomial confidence intervals, will be showed by week
variants_ch$date <- as_date(variants_ch$date)
variants_ch$date_num <- as.numeric(variants_ch$date)
variants_ch$week <- week(variants_ch$date)

variant_week <- aggregate(value~region+week+who_variants, variants_ch, sum)
seq_week <- aggregate(sequences~region+week+who_variants, variants_ch, sum)
variant_week<- merge(variant_week, seq_week,by=c("region","week","who_variants"))
remove(seq_week)
variant_week$sequences <- as.numeric(variant_week$sequences)
variant_week$value <- as.numeric(variant_week$value)
variant_week <- variant_week[variant_week$sequences!=0,]
variant_week$who_variants <- factor(variant_week$who_variants, levels = c(lev[-c(c36,others)],lev[c(c36,others)]))

lower <- numeric()
upper <- numeric()
conf <- numeric()
sampled_ci <- function(data){
  for(i in 1:length(data$who_variants)) {
    int <- binom.test(data$value[i], data$sequences[i])
    lower[i] <- int$conf.int[1]
    upper[i] <- int$conf.int[2]
    conf[i] <- int$estimate
  }
  data <- cbind(data, conf=conf, lower = lower, upper = upper)
 return(data)
}
variant_week <- sampled_ci(variant_week)

variant_week$week_day <- as_date(as.Date(paste(2021, variant_week$week, 1, sep="-"),format= "%Y-%U-%u"))
variant_week$region <- gsub("_"," ", variant_week$region)
variant_week$region <- gsub("r","R", variant_week$region)
variant_week$region <- factor(variant_week$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6", "CH"))
variant_week$week_day <- as_date(as.Date(paste(2021, variant_week$week, 1, sep="-"), "%Y-%U-%u"))
remove(conf)
remove(upper)
remove(lower)


### prepare data for multinomial function and plot
## overall Swiss data
variants_ch <- variants_ch[variants_ch$region=="CH",]
variants_ch$region <- factor(variants_ch$region, levels = c("CH"))
variants_ch <- variants_ch[variants_ch$value>0,]
variants_ch1 <- variants_ch[rep(row.names(variants_ch), variants_ch$value), c(2,6,7)]
variants_ch1$who_variants <- factor(variants_ch1$who_variants, levels = c(lev[-c(c36,others)],lev[c(c36,others)]))

## regional data
variants_reg <- variants_ch[variants_ch$region!="CH",]
variants_reg$region <- gsub("_"," ", variants_reg$region)
variants_reg$region <- gsub("r","R", variants_reg$region)
variants_reg$region <- factor(variants_reg$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6"))
variants_reg <- variants_reg[variants_reg$value>0,]
variants_reg1 <- variants_reg[rep(row.names(variants_reg), variants_reg$value), c(2,6,7)]
variants_reg1$who_variants <- factor(variants_reg1$who_variants, levels = c(lev[-c(c36,others)],lev[c(c36,others)]))
variants_reg1$region <- factor(variants_reg1$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6"))


### Different multi-nominal models:
mnom_date_spline <- multinom(who_variants ~ ns(date_num, df=2), data = variants_ch1)
mnom_date_reg_spline <- multinom(who_variants ~ ns(date_num, df=2) + region, data = variants_reg1)

## new data frame for prediction
predict_date <- data.frame(date_num = seq((min(variants_ch$date_num)),as.numeric(time_window[2]),1))
predict_date_reg <- data.frame(date_num=rep(unique(seq((min(variants_reg1$date_num)),as.numeric(time_window[2]),1)), length(unique(variants_reg1$region))), 
                               region=rep(levels(variants_reg1$region),each=length(unique(seq((min(variants_reg1$date_num)),as.numeric(time_window[2]),1)))))

## predict for mnom_week model
predict.eff_date <- Effect("date_num", mnom_date_spline,level=0.95,
                           xlevels=list(date_num=seq(min(variants_ch$date_num),as.numeric(time_window[2]),1)))
predict.eff_date_plot <- cbind(predict_date, melt(predict.eff_date$prob), melt(predict.eff_date$lower.prob),melt(predict.eff_date$upper.prob))
predict.eff_date_plot <- predict.eff_date_plot[,-c(2,5,6,8,9)]
colnames(predict.eff_date_plot) <- c("date_num","who_variants", "prob","lower", "upper")
predict.eff_date_plot$who_variants <- gsub("prob.", "", predict.eff_date_plot$who_variants)
predict.eff_date_plot$who_variants <- sapply(predict.eff_date_plot$who_variants, who_variant_names)
predict.eff_date_plot$who_variants <- factor(predict.eff_date_plot$who_variants, levels = c(lev[-c(c36,others)],lev[c(c36,others)]))
predict.eff_date_plot$date <- as_date(predict.eff_date_plot$date_num)

## predict for mnom_week_reg model
predict.eff_date_reg <- Effect(c("date_num","region"), mnom_date_reg_spline,level=0.95,
                               xlevels=list(date_num=seq(min(variants_reg$date_num),as.numeric(time_window[2]),1)))
predict.eff_date_reg_plot <- cbind(predict_date_reg[,c(1,2)],melt(predict.eff_date_reg$prob), melt(predict.eff_date_reg$lower.prob),melt(predict.eff_date_reg$upper.prob))
predict.eff_date_reg_plot <- predict.eff_date_reg_plot[,-c(3,6,7,9,10)]
colnames(predict.eff_date_reg_plot) <- c("date_num","region", "who_variants", "prob","lower", "upper")
predict.eff_date_reg_plot$who_variants <- gsub("prob.", "", predict.eff_date_reg_plot$who_variants)
predict.eff_date_reg_plot$who_variants <- sapply(predict.eff_date_reg_plot$who_variants, who_variant_names)
predict.eff_date_reg_plot$who_variants <- factor(predict.eff_date_reg_plot$who_variants, levels = c(lev[-c(c36,others)],lev[c(c36,others)]))
predict.eff_date_reg_plot$region <- gsub("_"," ", predict.eff_date_reg_plot$region)
predict.eff_date_reg_plot$region <- gsub("r","R", predict.eff_date_reg_plot$region)
predict.eff_date_reg_plot$region <- factor(predict.eff_date_reg_plot$region, levels = c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6"))
predict.eff_date_reg_plot$date <- as_date(predict.eff_date_reg_plot$date_num)



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
  geom_line(data=predict.eff_date_plot, aes(x = date, y = prob, color = who_variants))+
  geom_ribbon(data=predict.eff_date_plot, aes(x = date, y = prob, ymin = lower,ymax = upper,fill=who_variants),alpha=0.4)+
  geom_errorbar(data= variant_week[variant_week$region=="CH",], aes(x = week_day, ymin=lower, ymax=upper, color = who_variants), width=.1) +
  geom_point(data= variant_week[variant_week$region=="CH",],aes(x = week_day, y=conf, color = who_variants))+
  geom_rect(aes(xmin = as_date(max(variants_ch$date)), ymin = 0, xmax = time_window[2], ymax = 1), fill= col_9[9], colour= "transparent", alpha=0.2)+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window[1],time_window[2]+1)))+
  scale_color_manual(values= col_9[2:7],name="SARS-CoV-2 variants") +
  scale_fill_manual(values= col_9[2:7],name="SARS-CoV-2 variants") +
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
  scale_y_continuous(trans='log10',labels=yscaling,limits = c(min(predict.eff_date_plot$lower),10^0))+
  labs(y =bquote("Proportion of SARS-CoV-2 variants (log10-scale)"))
variants_legend <- g_legend(variants_plot_model,1)
variants_plot_model <- variants_plot_model + theme(legend.position = "none")

## reginal Swiss multi-nominal figure
regional_variants_plot_model <- ggplot(predict.eff_date_reg_plot) + 
  geom_line(data= predict.eff_date_reg_plot, aes(x = date, y = prob, color = who_variants))+
  geom_ribbon(data= predict.eff_date_reg_plot, aes(x= date, ymin = lower,ymax = upper,fill=who_variants),alpha=0.4)+
  geom_errorbar(data= variant_week[variant_week$region!="CH",], aes(x = week_day, ymin=lower, ymax=upper, color = who_variants), width=.1) +
  geom_point(data= variant_week[variant_week$region!="CH",],aes(x = week_day, y=conf, color = who_variants))+#, size = conf
  geom_rect(aes(xmin = as_date(max(variants_ch$date)), ymin = 0, xmax = time_window[2], ymax = 1), fill= col_9[9], colour= "transparent", alpha=0.008)+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window[1],time_window[2]+1)))+
  scale_color_manual(values= col_9[2:7],name="SARS-CoV-2 variants") +
  scale_fill_manual(values= col_9[2:7],name="SARS-CoV-2 variants") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))+
  facet_wrap(~ region, ncol=2)+
  scale_y_continuous(limits = c(0,1))+
  coord_cartesian(expand = FALSE)+
  labs(tag=bquote(.("")), x = "", y =bquote("Proportion of SARS-CoV-2 variants"))
regional_variants_plot_model_log <- regional_variants_plot_model+
  coord_cartesian(ylim = c(10^-3, 10^0))+
  scale_y_continuous(trans='log10',labels=yscaling, limits = c(10^-3,10^0))+
  labs(y =bquote("Proportion of SARS-CoV-2 variants (log10-scale)"))
com_regional_ch_variants<- grid.arrange(variants_plot_model, variants_plot_model_log, 
                                        variants_legend,regional_variants_plot_model, 
                                        regional_variants_plot_model_log, nrow = 2)
ggsave(com_regional_ch_variants, filename = paste0("../../plots/com_regional_ch_variants_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 9, width = 14,  bg = "transparent")
SequencingStats_and_Variants/scripts/R
