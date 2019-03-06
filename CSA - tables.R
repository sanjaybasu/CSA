setwd("~/Documents/Epi/Research/SDH/CSA")
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

# iters = 100
# eventtabs = matrix(0,5*iters,4)
# start_time <- Sys.time()
# for (i in 1:iters){
# eventtabs[((5*i-4):(5*i)),1:4] = CSArun(i)  
# }
# end_time <- Sys.time()
# end_time-start_time
# 
# # eventtab
# # cols = #rx1_deltatab_10,rx1_deltatab_life,rx2_deltatab_10,rx2_deltatab_life
# # rows = ascvd_events_10,dminc_events_10,esrd_events_10,neuro_events_10,retin_events_10
# save(eventtabs, file="events.RData")
# 
# intervention = c(rep("CSA",iters*5),rep("cash",iters*5))
# outcomeit = c(rep("ASCVD",iters),rep("DM inc",iters),rep("ESRD",iters),rep("Neuro",iters),rep("Retin",iters))
# outcome = rep(outcomeit,2)
# 
# rx1_ascvd_evts_prev_10 = (eventtabs[seq(1,dim(eventtabs)[1],5),1])
# rx1_dminc_evts_prev_10 = (eventtabs[seq(2,dim(eventtabs)[1],5),1])
# rx1_esrd_evts_prev_10 = (eventtabs[seq(3,dim(eventtabs)[1],5),1])
# rx1_neuro_evts_prev_10 = (eventtabs[seq(4,dim(eventtabs)[1],5),1])
# rx1_retin_evts_prev_10 = (eventtabs[seq(5,dim(eventtabs)[1],5),1])
# 
# rx2_ascvd_evts_prev_10 = (eventtabs[seq(1,dim(eventtabs)[1],5),3])
# rx2_dminc_evts_prev_10 = (eventtabs[seq(2,dim(eventtabs)[1],5),3])
# rx2_esrd_evts_prev_10 = (eventtabs[seq(3,dim(eventtabs)[1],5),3])
# rx2_neuro_evts_prev_10 = (eventtabs[seq(4,dim(eventtabs)[1],5),3])
# rx2_retin_evts_prev_10 = (eventtabs[seq(5,dim(eventtabs)[1],5),3])
# 
# rx1_ascvd_evts_prev_life = (eventtabs[seq(1,dim(eventtabs)[1],5),2])
# rx1_dminc_evts_prev_life = (eventtabs[seq(2,dim(eventtabs)[1],5),2])
# rx1_esrd_evts_prev_life = (eventtabs[seq(3,dim(eventtabs)[1],5),2])
# rx1_neuro_evts_prev_life = (eventtabs[seq(4,dim(eventtabs)[1],5),2])
# rx1_retin_evts_prev_life = (eventtabs[seq(5,dim(eventtabs)[1],5),2])
# 
# rx2_ascvd_evts_prev_life = (eventtabs[seq(1,dim(eventtabs)[1],5),4])
# rx2_dminc_evts_prev_life = (eventtabs[seq(2,dim(eventtabs)[1],5),4])
# rx2_esrd_evts_prev_life = (eventtabs[seq(3,dim(eventtabs)[1],5),4])
# rx2_neuro_evts_prev_life = (eventtabs[seq(4,dim(eventtabs)[1],5),4])
# rx2_retin_evts_prev_life = (eventtabs[seq(5,dim(eventtabs)[1],5),4])
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
# result_10 = c(rx1_ascvd_evts_prev_10,rx1_dminc_evts_prev_10,rx1_esrd_evts_prev_10,rx1_neuro_evts_prev_10,rx1_retin_evts_prev_10,
#                        rx2_ascvd_evts_prev_10,rx2_dminc_evts_prev_10,rx2_esrd_evts_prev_10,rx2_neuro_evts_prev_10,rx2_retin_evts_prev_10)
# 
# 
# result_life = c(rx1_ascvd_evts_prev_life,rx1_dminc_evts_prev_life,rx1_esrd_evts_prev_life,rx1_neuro_evts_prev_life,rx1_retin_evts_prev_life,
#                 rx2_ascvd_evts_prev_life,rx2_dminc_evts_prev_life,rx2_esrd_evts_prev_life,rx2_neuro_evts_prev_life,rx2_retin_evts_prev_life)
# 
# 
# 
# plotdat_10 = data.frame(intervention = intervention,
#                         outcome = outcome,
#                         result = -result_10)
# plotdat_life = data.frame(intervention = intervention,
#                         outcome = outcome,
#                         result = -result_life)
# library(ggplot2)
# # Figure 2 ####
# p <- ggplot(data = plotdat_10, aes(x=outcome, y=result, fill=intervention)) + 
#   geom_boxplot()
# p +scale_fill_brewer(palette="Accent")+ ylim(0, 250)+ 
#   ggtitle("Outcomes averted over 10-year horizon") +
#   xlab("Outcomes") + ylab("Number averted")
# 
# q <- ggplot(data = plotdat_life, aes(x=outcome, y=result, fill=intervention)) + 
#   geom_boxplot()
# q +scale_fill_brewer(palette="Accent")+ ylim(0, 250)+ 
#   ggtitle("Outcomes averted over life-course horizon") +
#   xlab("Outcomes") + ylab("Number averted")

iters = 100
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
save(qalytabs, file="qalys.RData")

base_qalys_ascvd_10 = summary(qalytabs[seq(1,dim(eventtabs)[1],6),1])
base_qalys_dminc_10 = summary(qalytabs[seq(1,dim(eventtabs)[1],6),2])
base_qalys_esrd_10 = summary(qalytabs[seq(1,dim(eventtabs)[1],6),3])
base_qalys_retin_10 = summary(qalytabs[seq(1,dim(eventtabs)[1],6),4])
base_qalys_neuro_10 = summary(qalytabs[seq(1,dim(eventtabs)[1],6),5])
base_qalys_allmort_10 = summary(qalytabs[seq(1,dim(eventtabs)[1],6),6])


