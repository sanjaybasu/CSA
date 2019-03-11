setwd("~/Documents/OneDrive - Leland Stanford Junior University/Documents/Epi/Research/SDH/CSA")
#install.packages('stargazer')
#install.packages('ggplot2')
load('merged.RData')

mergedt$eligible = ((mergedt$INDFMPIR<2) | (mergedt$HIQ031D==17))& (mergedt$BMXBMI>=25)
table(mergedt$eligible)
mergedsub = mergedt[which(mergedt$eligible==1),]
age = mergedsub$RIDAGEYR
female = (mergedsub$RIAGENDR=="Female")
black = (mergedsub$RIDRETH1=="Non-Hispanic Black")
black[is.na(black)]=0
hisp = (mergedsub$RIDRETH1=="Other Hispanic")
hisp[is.na(hisp)]=0
income = mergedsub$INDFMPIR
bmi = mergedsub$BMXBMI
bmi[is.na(bmi)]=mean(na.omit(bmi))
totchol = mergedsub$LBXTC
totchol[is.na(totchol)]=mean(na.omit(totchol))
hdlc = mergedsub$LBDHDD
hdlc[is.na(hdlc)]=mean(na.omit(hdlc))
dm = (mergedsub$DIQ010==1) | 
  (mergedsub$LBXGH>=6.5) | 
  (mergedsub$LBXGLU>=126)
dm[is.na(dm)]=0
sercreat = mergedsub$LBXSCR
sercreat[is.na(sercreat)]=mean(na.omit(sercreat))
uralbcreat =  mergedsub$URDACT
uralbcreat[is.na(uralbcreat)]=mean(na.omit(uralbcreat))
sysbp = mergedsub$BPXSY1
sysbp[is.na(sysbp)]=mean(na.omit(sysbp))
rxbp = (mergedsub$BPQ040A==1)
rxbp[is.na(rxbp)]=0
cvdhist = (mergedsub$MCQ160F==1)|(mergedsub$MCQ160F==1)
cvdhist[is.na(cvdhist)]=0
hgba1c = mergedsub$LBXGH
hgba1c[is.na(hgba1c)]=mean(na.omit(hgba1c))
cursmoke = (mergedsub$SMQ040==1 | mergedsub$SMQ040==2)
cursmoke[is.na(cursmoke)]=0
library(stringr)
statin = grepl("STATIN",mergedsub$RXDDRUG)
statin[is.na(statin)]=0
oralrx = grepl("METFORMIN",mergedsub$RXDDRUG)|grepl("GLIPIZIDE",mergedsub$RXDDRUG)
oralrx[is.na(oralrx)]=0
anticoag = grepl("WARFARIN",mergedsub$RXDDRUG)
anticoag[is.na(anticoag)]=0


devtools::install_github("timfolsom/hei")
library(hei)
fped1314 <- get_fped("2013/2014", "both")
diet1314 <- get_diet("2013/2014", "both")
demo1314 <- get_demo("2013/2014")
eligfped = fped1314[fped1314$SEQN %in% mergedsub$SEQN,]
eligdiet = diet1314[diet1314$SEQN %in% mergedsub$SEQN,]
eligdemo = demo1314[demo1314$SEQN %in% mergedsub$SEQN,]
heiscoreelig = hei(eligfped,eligdiet,eligdemo)
hei = heiscoreelig$HEI

# Table 2 ####
dat = data.frame(cbind(age,female,black,hisp,income,hei,bmi,sysbp,totchol,hdlc,dm,hgba1c,sercreat,uralbcreat,cursmoke,cvdhist,rxbp,statin,oralrx,anticoag))
stargazer::stargazer(dat,type="text", out="Table 2.html")

source("CSArun.R")
#risktabs = CSArun(1)

# iters = 10000
# eventtabs = matrix(0,6*iters,4)
# start_time <- Sys.time()
# for (i in 1:iters){
# eventtabs[((6*i-5):(6*i)),1:4] = CSArun(i)
# }
# end_time <- Sys.time()
# end_time-start_time
# 
# # eventtab
# # cols = #rx1_deltatab_10,rx1_deltatab_life,rx2_deltatab_10,rx2_deltatab_life
# # rows = ascvd_events_10,dminc_events_10,esrd_events_10,neuro_events_10,retin_events_10
# save(eventtabs, file="events.RData")
# 
# intervention1 = c(rep("CSA",iters*6),rep("cash",iters*6))
# outcomeit1 = c(rep("ASCVD",iters),rep("DM inc",iters),rep("ESRD",iters),rep("Neuro",iters),rep("Retin",iters),rep("Deaths",iters))
# outcome1 = rep(outcomeit1,2)
# 
# rx1_ascvd_evts_prev_10 = (eventtabs[seq(1,dim(eventtabs)[1],6),1])
# rx1_dminc_evts_prev_10 = (eventtabs[seq(2,dim(eventtabs)[1],6),1])
# rx1_esrd_evts_prev_10 = (eventtabs[seq(3,dim(eventtabs)[1],6),1])
# rx1_neuro_evts_prev_10 = (eventtabs[seq(4,dim(eventtabs)[1],6),1])
# rx1_retin_evts_prev_10 = (eventtabs[seq(5,dim(eventtabs)[1],6),1])
# rx1_mort_evts_prev_10 = (eventtabs[seq(6,dim(eventtabs)[1],6),1])
# 
# rx2_ascvd_evts_prev_10 = (eventtabs[seq(1,dim(eventtabs)[1],6),3])
# rx2_dminc_evts_prev_10 = (eventtabs[seq(2,dim(eventtabs)[1],6),3])
# rx2_esrd_evts_prev_10 = (eventtabs[seq(3,dim(eventtabs)[1],6),3])
# rx2_neuro_evts_prev_10 = (eventtabs[seq(4,dim(eventtabs)[1],6),3])
# rx2_retin_evts_prev_10 = (eventtabs[seq(5,dim(eventtabs)[1],6),3])
# rx2_mort_evts_prev_10 = (eventtabs[seq(6,dim(eventtabs)[1],6),3])
# 
# rx1_ascvd_evts_prev_life = (eventtabs[seq(1,dim(eventtabs)[1],6),2])
# rx1_dminc_evts_prev_life = (eventtabs[seq(2,dim(eventtabs)[1],6),2])
# rx1_esrd_evts_prev_life = (eventtabs[seq(3,dim(eventtabs)[1],6),2])
# rx1_neuro_evts_prev_life = (eventtabs[seq(4,dim(eventtabs)[1],6),2])
# rx1_retin_evts_prev_life = (eventtabs[seq(5,dim(eventtabs)[1],6),2])
# 
# rx2_ascvd_evts_prev_life = (eventtabs[seq(1,dim(eventtabs)[1],6),4])
# rx2_dminc_evts_prev_life = (eventtabs[seq(2,dim(eventtabs)[1],6),4])
# rx2_esrd_evts_prev_life = (eventtabs[seq(3,dim(eventtabs)[1],6),4])
# rx2_neuro_evts_prev_life = (eventtabs[seq(4,dim(eventtabs)[1],6),4])
# rx2_retin_evts_prev_life = (eventtabs[seq(5,dim(eventtabs)[1],6),4])
# 
# summary(rx1_ascvd_evts_prev_life)
# summary(rx1_dminc_evts_prev_life)
# summary(rx1_esrd_evts_prev_life)
# summary(rx1_neuro_evts_prev_life)
# summary(rx1_retin_evts_prev_life)
# 
# summary(rx2_ascvd_evts_prev_life)
# summary(rx2_dminc_evts_prev_life)
# summary(rx2_esrd_evts_prev_life)
# summary(rx2_neuro_evts_prev_life)
# summary(rx2_retin_evts_prev_life)
# 
# 
# 
# 
# result_10 = c(rx1_ascvd_evts_prev_10,rx1_dminc_evts_prev_10,rx1_esrd_evts_prev_10,rx1_neuro_evts_prev_10,rx1_retin_evts_prev_10,rx1_mort_evts_prev_10,
#                        rx2_ascvd_evts_prev_10,rx2_dminc_evts_prev_10,rx2_esrd_evts_prev_10,rx2_neuro_evts_prev_10,rx2_retin_evts_prev_10,rx2_mort_evts_prev_10)
# 
# 
# result_life = c(rx1_ascvd_evts_prev_life,rx1_dminc_evts_prev_life,rx1_esrd_evts_prev_life,rx1_neuro_evts_prev_life,rx1_retin_evts_prev_life,
#                 rx2_ascvd_evts_prev_life,rx2_dminc_evts_prev_life,rx2_esrd_evts_prev_life,rx2_neuro_evts_prev_life,rx2_retin_evts_prev_life)
# 
# 
# 
# plotdat_10 = data.frame(intervention = intervention1,
#                         outcome = outcome1,
#                         result = -result_10)
# plotdat_life = data.frame(intervention = intervention1,
#                         outcome = outcome1,
#                         result = -result_life)
# library(ggplot2)
# # Figure 2 ####
# plotdat_10$outcome2 = factor(plotdat_10$outcome, levels=c("ASCVD", "DM inc", "ESRD", "Neuro", "Retin", "Deaths"))
# p = ggplot(data = plotdat_10, aes(x=outcome2, y=result, fill=intervention)) +
#   geom_boxplot()
# p +scale_fill_brewer(palette="Accent")+ ylim(0, 250)+
#   ggtitle("Outcomes averted over 10-year horizon") +
#   xlab("Outcomes") + ylab("Number averted")
# 
# q <- ggplot(subset(plotdat_life, outcome %in% c("ASCVD", "DM inc", "ESRD", "Neuro", "Retin")), aes(x=outcome, y=result, fill=intervention)) +
#   geom_boxplot()
# q +scale_fill_brewer(palette="Accent")+ ylim(0, 250)+
#   ggtitle("Outcomes averted over life-course horizon") +
#   xlab("Outcomes") + ylab("Number averted")



# Table 3 #####

iters = 1
source("CSArun.R")

qalytabs = matrix(0,6*iters,6)
start_time <- Sys.time()
for (i in 1:iters){
  qalytabs[((6*i-5):(6*i)),1:6] = CSArun(i)
}
end_time <- Sys.time()
end_time-start_time

# qalytab
# cols = #rx1_deltatab_10,rx1_deltatab_life,rx2_deltatab_10,rx2_deltatab_life
# rows = ascvd_qalys_10,dminc_qalys_10,esrd_qalys_10,neuro_qalys_10,retin_qalys_10
# save(qalytabs, file="qalys.RData")
# 
# 
base_qalys_ascvd_10 = summary(qalytabs[seq(1,dim(qalytabs)[1],6),1])
base_qalys_dminc_10 = summary(qalytabs[seq(1,dim(qalytabs)[1],6),2])
base_qalys_esrd_10 = summary(qalytabs[seq(1,dim(qalytabs)[1],6),3])
base_qalys_retin_10 = summary(qalytabs[seq(1,dim(qalytabs)[1],6),4])
base_qalys_neuro_10 = summary(qalytabs[seq(1,dim(qalytabs)[1],6),5])
base_qalys_allmort_10 = summary(qalytabs[seq(1,dim(qalytabs)[1],6),6])
base_qalys_ascvd_10+base_qalys_dminc_10+base_qalys_esrd_10+
      base_qalys_retin_10+base_qalys_neuro_10+base_qalys_allmort_10

rx1_qalys_ascvd_10 = summary(qalytabs[seq(2,dim(qalytabs)[1],6),1])
rx1_qalys_dminc_10 = summary(qalytabs[seq(2,dim(qalytabs)[1],6),2])
rx1_qalys_esrd_10 = summary(qalytabs[seq(2,dim(qalytabs)[1],6),3])
rx1_qalys_retin_10 = summary(qalytabs[seq(2,dim(qalytabs)[1],6),4])
rx1_qalys_neuro_10 = summary(qalytabs[seq(2,dim(qalytabs)[1],6),5])
rx1_qalys_allmort_10 = summary(qalytabs[seq(2,dim(qalytabs)[1],6),6])
rx1_qalys_ascvd_10+rx1_qalys_dminc_10+rx1_qalys_esrd_10+
  rx1_qalys_retin_10+rx1_qalys_neuro_10+rx1_qalys_allmort_10

rx2_qalys_ascvd_10 = summary(qalytabs[seq(3,dim(qalytabs)[1],6),1])
rx2_qalys_dminc_10 = summary(qalytabs[seq(3,dim(qalytabs)[1],6),2])
rx2_qalys_esrd_10 = summary(qalytabs[seq(3,dim(qalytabs)[1],6),3])
rx2_qalys_retin_10 = summary(qalytabs[seq(3,dim(qalytabs)[1],6),4])
rx2_qalys_neuro_10 = summary(qalytabs[seq(3,dim(qalytabs)[1],6),5])
rx2_qalys_allmort_10 = summary(qalytabs[seq(3,dim(qalytabs)[1],6),6])
rx2_qalys_ascvd_10+rx2_qalys_dminc_10+rx2_qalys_esrd_10+
  rx2_qalys_retin_10+rx2_qalys_neuro_10+rx2_qalys_allmort_10

base_qalys_ascvd_life = summary(qalytabs[seq(4,dim(qalytabs)[1],6),1])
base_qalys_dminc_life = summary(qalytabs[seq(4,dim(qalytabs)[1],6),2])
base_qalys_esrd_life = summary(qalytabs[seq(4,dim(qalytabs)[1],6),3])
base_qalys_retin_life = summary(qalytabs[seq(4,dim(qalytabs)[1],6),4])
base_qalys_neuro_life = summary(qalytabs[seq(4,dim(qalytabs)[1],6),5])
base_qalys_allmort_life = summary(qalytabs[seq(4,dim(qalytabs)[1],6),6])
base_qalys_ascvd_life+base_qalys_dminc_life+base_qalys_esrd_life+
  base_qalys_retin_life+base_qalys_neuro_life+base_qalys_allmort_life

rx1_qalys_ascvd_life = summary(qalytabs[seq(5,dim(qalytabs)[1],6),1])
rx1_qalys_dminc_life = summary(qalytabs[seq(5,dim(qalytabs)[1],6),2])
rx1_qalys_esrd_life = summary(qalytabs[seq(5,dim(qalytabs)[1],6),3])
rx1_qalys_retin_life = summary(qalytabs[seq(5,dim(qalytabs)[1],6),4])
rx1_qalys_neuro_life = summary(qalytabs[seq(5,dim(qalytabs)[1],6),5])
rx1_qalys_allmort_life = summary(qalytabs[seq(5,dim(qalytabs)[1],6),6])
rx1_qalys_ascvd_life+rx1_qalys_dminc_life+rx1_qalys_esrd_life+
  rx1_qalys_retin_life+rx1_qalys_neuro_life+rx1_qalys_allmort_life

rx2_qalys_ascvd_life = summary(qalytabs[seq(6,dim(qalytabs)[1],6),1])
rx2_qalys_dminc_life = summary(qalytabs[seq(6,dim(qalytabs)[1],6),2])
rx2_qalys_esrd_life = summary(qalytabs[seq(6,dim(qalytabs)[1],6),3])
rx2_qalys_retin_life = summary(qalytabs[seq(6,dim(qalytabs)[1],6),4])
rx2_qalys_neuro_life = summary(qalytabs[seq(6,dim(qalytabs)[1],6),5])
rx2_qalys_allmort_life = summary(qalytabs[seq(6,dim(qalytabs)[1],6),6])
rx2_qalys_ascvd_life+rx2_qalys_dminc_life+rx2_qalys_esrd_life+
rx2_qalys_retin_life+rx2_qalys_neuro_life+rx2_qalys_allmort_life

iters = 10000
source("CSArun.R")

costtabs = matrix(0,6*iters,8)
start_time <- Sys.time()
for (i in 1:iters){
  costtabs[((6*i-5):(6*i)),1:8] = CSArun(i)
}
end_time <- Sys.time()
end_time-start_time

# costtab
# cols = #rx1_deltatab_10,rx1_deltatab_life,rx2_deltatab_10,rx2_deltatab_life
# rows = ascvd_costs_10,dminc_costs_10,esrd_costs_10,neuro_costs_10,retin_costs_10
save(costtabs, file="costs.RData")
#load("costs.RData")
base_costs_int_10 = summary(costtabs[seq(1,dim(costtabs)[1],6),1])
base_costs_ascvd_10 = summary(costtabs[seq(1,dim(costtabs)[1],6),2])
base_costs_dminc_10 = summary(costtabs[seq(1,dim(costtabs)[1],6),3])
base_costs_esrd_10 = summary(costtabs[seq(1,dim(costtabs)[1],6),4])
base_costs_retin_10 = summary(costtabs[seq(1,dim(costtabs)[1],6),5])
base_costs_neuro_10 = summary(costtabs[seq(1,dim(costtabs)[1],6),6])
base_costs_ag_10 = summary(costtabs[seq(1,dim(costtabs)[1],6),7])
base_costs_prod_10 = summary(costtabs[seq(1,dim(costtabs)[1],6),8])

base_costs_int_10+base_costs_ascvd_10+base_costs_dminc_10+base_costs_esrd_10+
  base_costs_retin_10+base_costs_neuro_10
base_costs_ag_10+base_costs_prod_10

rx1_costs_int_10 = summary(costtabs[seq(2,dim(costtabs)[1],6),1])
rx1_costs_ascvd_10 = summary(costtabs[seq(2,dim(costtabs)[1],6),2])
rx1_costs_dminc_10 = summary(costtabs[seq(2,dim(costtabs)[1],6),3])
rx1_costs_esrd_10 = summary(costtabs[seq(2,dim(costtabs)[1],6),4])
rx1_costs_retin_10 = summary(costtabs[seq(2,dim(costtabs)[1],6),5])
rx1_costs_neuro_10 = summary(costtabs[seq(2,dim(costtabs)[1],6),6])
rx1_costs_ag_10 = summary(costtabs[seq(2,dim(costtabs)[1],6),7])
rx1_costs_prod_10 = summary(costtabs[seq(2,dim(costtabs)[1],6),8])

rx1_costs_int_10+rx1_costs_ascvd_10+rx1_costs_dminc_10+rx1_costs_esrd_10+
  rx1_costs_retin_10+rx1_costs_neuro_10
rx1_costs_ag_10+rx1_costs_prod_10

rx2_costs_int_10 = summary(costtabs[seq(3,dim(costtabs)[1],6),1])
rx2_costs_ascvd_10 = summary(costtabs[seq(3,dim(costtabs)[1],6),2])
rx2_costs_dminc_10 = summary(costtabs[seq(3,dim(costtabs)[1],6),3])
rx2_costs_esrd_10 = summary(costtabs[seq(3,dim(costtabs)[1],6),4])
rx2_costs_retin_10 = summary(costtabs[seq(3,dim(costtabs)[1],6),5])
rx2_costs_neuro_10 = summary(costtabs[seq(3,dim(costtabs)[1],6),6])
rx2_costs_ag_10 = summary(costtabs[seq(3,dim(costtabs)[1],6),7])
rx2_costs_prod_10 = summary(costtabs[seq(3,dim(costtabs)[1],6),8])

rx2_costs_int_10+rx2_costs_ascvd_10+rx2_costs_dminc_10+rx2_costs_esrd_10+
  rx2_costs_retin_10+rx2_costs_neuro_10
rx2_costs_ag_10+rx2_costs_prod_10

base_costs_int_life = summary(costtabs[seq(4,dim(costtabs)[1],6),1])
base_costs_ascvd_life = summary(costtabs[seq(4,dim(costtabs)[1],6),2])
base_costs_dminc_life = summary(costtabs[seq(4,dim(costtabs)[1],6),3])
base_costs_esrd_life = summary(costtabs[seq(4,dim(costtabs)[1],6),4])
base_costs_retin_life = summary(costtabs[seq(4,dim(costtabs)[1],6),5])
base_costs_neuro_life = summary(costtabs[seq(4,dim(costtabs)[1],6),6])
base_costs_ag_life = summary(costtabs[seq(4,dim(costtabs)[1],6),7])
base_costs_prod_life = summary(costtabs[seq(4,dim(costtabs)[1],6),8])

base_costs_int_life+base_costs_ascvd_life+base_costs_dminc_life+base_costs_esrd_life+
  base_costs_retin_life+base_costs_neuro_life
base_costs_ag_life+base_costs_prod_life

rx1_costs_int_life = summary(costtabs[seq(5,dim(costtabs)[1],6),1])
rx1_costs_ascvd_life = summary(costtabs[seq(5,dim(costtabs)[1],6),2])
rx1_costs_dminc_life = summary(costtabs[seq(5,dim(costtabs)[1],6),3])
rx1_costs_esrd_life = summary(costtabs[seq(5,dim(costtabs)[1],6),4])
rx1_costs_retin_life = summary(costtabs[seq(5,dim(costtabs)[1],6),5])
rx1_costs_neuro_life = summary(costtabs[seq(5,dim(costtabs)[1],6),6])
rx1_costs_ag_life = summary(costtabs[seq(5,dim(costtabs)[1],6),7])
rx1_costs_prod_life = summary(costtabs[seq(5,dim(costtabs)[1],6),8])

rx1_costs_int_life+rx1_costs_ascvd_life+rx1_costs_dminc_life+rx1_costs_esrd_life+
  rx1_costs_retin_life+rx1_costs_neuro_life
rx1_costs_ag_life+rx1_costs_prod_life

rx2_costs_int_life = summary(costtabs[seq(6,dim(costtabs)[1],6),1])
rx2_costs_ascvd_life = summary(costtabs[seq(6,dim(costtabs)[1],6),2])
rx2_costs_dminc_life = summary(costtabs[seq(6,dim(costtabs)[1],6),3])
rx2_costs_esrd_life = summary(costtabs[seq(6,dim(costtabs)[1],6),4])
rx2_costs_retin_life = summary(costtabs[seq(6,dim(costtabs)[1],6),5])
rx2_costs_neuro_life = summary(costtabs[seq(6,dim(costtabs)[1],6),6])
rx2_costs_ag_life = summary(costtabs[seq(6,dim(costtabs)[1],6),7])
rx2_costs_prod_life = summary(costtabs[seq(6,dim(costtabs)[1],6),8])

rx2_costs_int_life+rx2_costs_ascvd_life+rx2_costs_dminc_life+rx2_costs_esrd_life+
  rx2_costs_retin_life+rx2_costs_neuro_life
rx2_costs_ag_life+rx2_costs_prod_life

icer_cash = ((rx1_costs_int_life+rx1_costs_ascvd_life+rx1_costs_dminc_life+rx1_costs_esrd_life+
    rx1_costs_retin_life+rx1_costs_neuro_life-(rx1_costs_ag_life+rx1_costs_prod_life))-(base_costs_int_life+base_costs_ascvd_life+base_costs_dminc_life+base_costs_esrd_life+
                                                  base_costs_retin_life+base_costs_neuro_life))/((rx2_qalys_ascvd_life+rx2_qalys_dminc_life+rx2_qalys_esrd_life+
                                                                                                    rx2_qalys_retin_life+rx2_qalys_neuro_life+rx2_qalys_allmort_life)-(base_qalys_ascvd_life+base_qalys_dminc_life+base_qalys_esrd_life+
                                                                                                                                                                         base_qalys_retin_life+base_qalys_neuro_life+base_qalys_allmort_life))
icer_CSA = ((rx2_costs_int_life+rx2_costs_ascvd_life+rx2_costs_dminc_life+rx2_costs_esrd_life+
    rx2_costs_retin_life+rx2_costs_neuro_life-(rx2_costs_ag_life+rx2_costs_prod_life))-(base_costs_int_life+base_costs_ascvd_life+base_costs_dminc_life+base_costs_esrd_life+
                                                  base_costs_retin_life+base_costs_neuro_life))/((rx1_qalys_ascvd_life+rx1_qalys_dminc_life+rx1_qalys_esrd_life+
                                                                                                    rx1_qalys_retin_life+rx1_qalys_neuro_life+rx1_qalys_allmort_life)-(base_qalys_ascvd_life+base_qalys_dminc_life+base_qalys_esrd_life+
                                                                                                                                                                         base_qalys_retin_life+base_qalys_neuro_life+base_qalys_allmort_life))

