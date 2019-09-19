# LD decay plot

system("~/programs/plink --bfile data/SNP_chip/sheep_geno_imputed_04092019 --sheep --r2 --ld-window-r2 0 --ld-window 5000 --ld-window-kb 5000 --out data/ld_decay/sheep_ld")

# run from command line
# cat data/ld_decay/sheep_ld.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > data/ld_decay/sheep_ld_summary 

library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
source("../sheep_ID/theme_clean.R")
sheep_ld <- fread("data/ld_decay/sheep_ld.ld")


# calculate distance between snps
sheep_ld[, dist_kb := (abs(BP_A - BP_B)/1000)] 

# dplyr approach
# sheep_ld_df <- sheep_ld %>% 
#                     mutate(dist_kb := (abs(BP_A - BP_B)/1000)) %>% 
#                     dplyr::filter(CHR_A == CHR_B) %>% 
#                     arrange(dist_kb)

# data.table_approach
# add dist
sheep_ld[, dist_kb := (abs(BP_A - BP_B)/1000)] 
# filter
sheep_ld <- sheep_ld[CHR_A == CHR_B] 
# sort
sheep_ld <- sheep_ld[order(dist_kb)]
# cut in groups
sheep_ld$distc <- cut(sheep_ld$dist_kb, seq(from=min(sheep_ld$dist_kb) - 0.001,to= max(sheep_ld$dist_kb), by=100))


# summarise per group
# sheep_ld_sum <- sheep_ld %>% group_by(distc) %>% summarise(mean=mean(R2),median=median(R2))
sheep_ld_sum <- sheep_ld %>% group_by(distc) %>% summarise(mean=mean(R2),median=median(R2),
                                                           ci_95_low = quantile(R2, probs = 0.025),
                                                           ci_95_high = quantile(R2, probs = 0.975),
                                                           ci_80_low = quantile(R2, probs = 0.1),
                                                           ci_80_high = quantile(R2, probs = 0.9),
                                                           ci_50_low = quantile(R2, probs = 0.250),
                                                           ci_50_high = quantile(R2, probs = 0.750))

dfr1 <-sheep_ld_sum %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

p1 <- ggplot(data=dfr1, aes(x=start,y=mean))+
    geom_errorbar(aes(ymin = ci_95_low, ymax = ci_95_high), width = 0,color = "#9ecae1", size = 0.5) +
    geom_errorbar(aes(ymin = ci_80_low, ymax = ci_80_high), width = 0,color = "#4292c6", size = 1) +
    geom_errorbar(aes(ymin = ci_50_low, ymax = ci_50_high), width = 0,color = "#08306b", size = 2) +
    geom_line(aes(x=start,y=mean),size=0.7,alpha=1,colour="black") +
    geom_point(colour="black",fill="#d9d9d9", shape = 21, size =3, stroke = 0.5) +
    theme_clean() +
    ylab("LD(r2)") +
    xlab("Distance (Kb)") +
   # theme(legend.position="right") +
    ggtitle("LD decay in Soay sheep", subtitle = "Mean and 50%, 80% and 95% CI")
p1
ggsave("figs/ld_decay_mean_uncertainty.jpg", p1, width = 7, height = 5)


p2 <- ggplot()+
    geom_point(data=dfr1,aes(x=start,y=median),size=0.4,colour="grey20")+
    geom_line(data=dfr1,aes(x=start,y=median),size=0.3,alpha=0.5,colour="grey40")+
    labs(x="Distance (Kilobases)",y=expression(LD~(r^{2})))+
    #scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8"))+
    theme_clean()
p2


