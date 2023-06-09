---
title: "Data prepara Template SS3"
subtitle: "Alternative Analysis to incorporate in Krill Stock Assessment Model 48.1 SubArea"
author: "Mardones, M; Watters, G.; Kinzey. D."
date:  "13 April, 2023"
linkcolor: blue
output:
  html_document:
    keep_md: true
    toc: true
    toc_deep: 3
    toc_float:
      collapsed: false
      smooth_scroll: false
    theme: cosmo
    fontsize: 0.9em
    linestretch: 1.7
    html-math-method: katex
    self-contained: true
    code-tools: true
editor_options: 
  markdown: 
    wrap: 72
---


```r
rm(list = ls())
knitr::opts_chunk$set(echo = TRUE,
                      message = FALSE,
                      warning = FALSE,
                      fig.align = 'center',
                      dev = 'jpeg',
                      dpi = 300)
#XQuartz is a mess, put this in your onload to default to cairo instead
options(bitmapType = "cairo") 
# (https://github.com/tidyverse/ggplot2/issues/2655)
# Lo mapas se hacen mas rapido
```



```r
library(tidyverse)
library(ggridges)
```

# Background

The following document intends to carry out a complementary
methodological analysis to correlate environmental variables with the
population dynamics of krill (*Euphausia superba*), in this case, with a
biological component like lengths from fishery monitoring.


# Methodology

Load data from Doug Kinzey code



Get data specific object



But, this structure is usefull for `SS3` template in `.dat`. 


```r
dat.all <- # net length frequencies including zero hauls
           read.csv(file=paste("AMLR_Krill_LFD_data.csv",sep=""),sep=",",
	   header=T,stringsAsFactors=F)
  dat.all$length[dat.all$length==0] <- NA # replace 840 zero lengths
  dat.all$amount[dat.all$amount==0] <- NA # replaced 683 zero amounts
dat.l <- subset(dat.all,dat.all$leg=="A" | dat.all$leg =="D") # A and D legs only
dat.l <- subset(dat.l,dat.l$amlr_area=="EI" | dat.l$amlr_area =="SA" |
                dat.l$amlr_area=="WA" | dat.l$amlr_area =="JI")
dat.l <-cbind(Year=substr(dat.l$AMLR_Cruise,5,8),dat.l)
dat.l <- dat.l[-which(is.na(dat.l$amount)),]
```

Expand frecuency data related length, in this case `amount` column have frecuency that we need expand to whole data frame. 


```r
df <- dat.l %>% 
  type.convert(as.is = TRUE) %>% 
  uncount(amount)
```



```r
jzstrata <- ggplot(df ,
                   aes(x=length, 
                       y = as.factor(Year), 
                       fill=amlr_area ))+
  geom_density_ridges(stat = "binline", bins = 30, 
                      scale = 1.9, 
                      draw_baseline = FALSE,
                      alpha=0.9)+
  facet_wrap(.~amlr_area , ncol=7) +   
  geom_vline(xintercept = 40, color = "red")+
  scale_x_continuous(breaks = seq(from = 10, to = 80, 
                                  by = 10))+
  scale_y_discrete(breaks = seq(from = 1991, 
                                to = 2011, by = 1))+
  scale_fill_viridis_d(name="Strata",
                       option="H")+
  
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 2))+
  xlab("Length (mm.)")+
  ylab("")
jzstrata
```

<img src="index_files/figure-html/unnamed-chunk-6-1.jpeg" style="display: block; margin: auto;" />
by leg


```r
legtrata <- ggplot(df ,
                   aes(x=length, 
                       y = as.factor(Year), 
                       fill=leg ))+
  geom_density_ridges(stat = "binline", bins = 30, 
                      scale = 1.9, 
                      draw_baseline = FALSE,
                      alpha=0.9)+
  facet_wrap(.~amlr_area , ncol=7) +   
  geom_vline(xintercept = 40, color = "red")+
  scale_x_continuous(breaks = seq(from = 10, to = 80, 
                                  by = 10))+
  scale_y_discrete(breaks = seq(from = 1991, 
                                to = 2011, by = 1))+
  scale_fill_viridis_d(name="Leg",
                       option="H")+
  
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 2))+
  xlab("Length (mm.)")+
  ylab("")
legtrata
```

<img src="index_files/figure-html/unnamed-chunk-7-1.jpeg" style="display: block; margin: auto;" />


Matrurity stage


```r
mattrata <- ggplot(df %>% 
                     filter(!maturity %in% c("FEM", "MALE", "UNKN")),
                   aes(x=length, 
                       y = as.factor(Year), 
                       fill=amlr_area ))+
  geom_density_ridges(stat = "binline", bins = 30, 
                      scale = 1.9, 
                      draw_baseline = FALSE,
                      alpha=0.9)+
  facet_wrap(.~maturity , ncol=7) +   
  geom_vline(xintercept = 40, color = "red")+
  scale_x_continuous(breaks = seq(from = 10, to = 80, 
                                  by = 10))+
  scale_y_discrete(breaks = seq(from = 1991, 
                                to = 2011, by = 1))+
  scale_fill_viridis_d(name="Strata",
                       option="G")+
  
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 2))+
  xlab("Length (mm.)")+
  ylab("")
mattrata
```

<img src="index_files/figure-html/unnamed-chunk-8-1.jpeg" style="display: block; margin: auto;" />
Give template structure by year.



```r
# cut in order
df$catlon <- cut(x = df$length, 
                 breaks = seq(0,70,2),
                 labels = seq(0,68,2),
                 right = FALSE)

dft <- table(df$Year, df$catlon)
```



```r
write.csv(dft, "lenghtAMLR19912011.csv", sep = ",", row.names = TRUE)
```

