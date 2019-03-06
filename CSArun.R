
CSArun <- function(iter) {
  setwd("~/Documents/Epi/Research/SDH/CSA")
  load('merged.RData')
  mergedt$eligible = ((mergedt$INDFMPIR<2) | (mergedt$HIQ031D==17)) & (mergedt$BMXBMI>=25)
  table(mergedt$eligible)
  mergedsub = mergedt[which(mergedt$eligible==1),]
  age = mergedsub$RIDAGEYR
  female = (mergedsub$RIAGENDR=="Female")
  black = (mergedsub$RIDRETH1=="Non-Hispanic Black")
  black[is.na(black)]=0
  hisp = (mergedsub$RIDRETH1=="Other Hispanic")
  hisp[is.na(hisp)]=0
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
  ascvd_10   = ( 1.0*female / (1.0 + exp( - (
    -12.823110 +
      0.106501 * as.numeric(age) +
      0.432440 * as.numeric(black) +
      0.000056 * (as.numeric(sysbp) ^ 2) +
      0.017666 * as.numeric(sysbp) +
      0.731678 * as.numeric(rxbp) +
      0.943970 * as.numeric(dm) +
      1.009790 * as.numeric(cursmoke) +
      0.151318 * (as.numeric(totchol) / as.numeric(hdlc)) +
      -0.008580 * as.numeric(age) * as.numeric(black) +
      -0.003647 * as.numeric(sysbp) * as.numeric(rxbp) +
      0.006208 * as.numeric(sysbp) * as.numeric(black) +
      0.152968 * as.numeric(black) * as.numeric(rxbp) +
      -0.000153 * as.numeric(age) * as.numeric(sysbp) +
      0.115232 * as.numeric(black) * as.numeric(dm) +
      -0.092231 * as.numeric(black) * as.numeric(cursmoke) +
      0.070498 * as.numeric(black) * (as.numeric(totchol) / as.numeric(hdlc)) +
      -0.000173 * as.numeric(black)  * as.numeric(sysbp) * as.numeric(rxbp) +
      -0.000094 * as.numeric(age) * as.numeric(sysbp) * as.numeric(black)))))+
    (1.0*(1-female) / (1.0 + exp( - (
      -11.679980 +
        0.064200 * as.numeric(age) +
        0.482835 * as.numeric(black) +
        -0.000061 * (as.numeric(sysbp) ^ 2) +
        0.038950 * as.numeric(sysbp) +
        2.055533 * as.numeric(rxbp) +
        0.842209 * as.numeric(dm) +
        0.895589 * as.numeric(cursmoke) +
        0.193307 * (as.numeric(totchol) / as.numeric(hdlc)) +
        -0.014207 * as.numeric(sysbp) * as.numeric(rxbp) +
        0.011609 * as.numeric(sysbp) * as.numeric(black) +
        -0.119460 * as.numeric(rxbp) * as.numeric(black) +
        0.000025 * as.numeric(age) * as.numeric(sysbp) +
        -0.077214 * as.numeric(black) * as.numeric(dm) +
        -0.226771 * as.numeric(black) * as.numeric(cursmoke) +
        -0.117749 * (as.numeric(totchol) / as.numeric(hdlc)) * as.numeric(black) +
        0.004190 * as.numeric(black) * as.numeric(rxbp) * as.numeric(sysbp) +
        -0.000199 * as.numeric(black) * as.numeric(age) * as.numeric(sysbp)))))
    dminc_10 = (mean(c(.031,.068))*(age<45)+
                mean(c(.109,.068))*(age>=45 & age<65)+
                mean(c(.094,.068))*(age>=65))*100/1000*10
  nbetax = (-1.938e-02*age+
              -1.129e-02*female+
              -8.812e-02*black+
              2.338e-01*hisp+
              1.483e-01*cursmoke+
              3.027e-03*sysbp+
              -7.952e-02*rxbp+
              -1.256e-01*oralrx+
              3.199e-02*anticoag+
              -2.164e-02*cvdhist+
              1.369e-01*hgba1c+
              -1.112e-03*totchol+
              6.289e-03*hdlc+
              8.609e-01*sercreat+
              3.615e-04*uralbcreat)
  esrdrisk_10 = 1-.973^exp(nbetax-mean(nbetax,na.rm=T))
  rbetax = (2.285e-02*age+
              2.264e-01*female+
              -1.677e-01*black+
              8.243e-03*sysbp+
              6.393e-02*rxbp+
              -2.349e-01*oralrx+
              1.127e-01*cvdhist+
              1.449e-01*hgba1c+
              -1.676e-04*totchol+
              5.447e-03*hdlc+
              6.947e-01*sercreat+
              1.992e-04*uralbcreat)
  retinrisk_10 = 1 - .92^exp(rbetax-mean(rbetax,na.rm=T))
  ebetax = (0.0302237*age+
              -0.1868000*female+
              -0.0944841*black+
              0.0045609*sysbp+
              0.1819157*rxbp+
              -0.2574724*oralrx+
              0.2667152*cvdhist+
              0.1886579*hgba1c+
              0.0021850*totchol+
              -0.0053887*hdlc+
              0.6044183*sercreat)
  neurorisk_10 = 1 - .87^exp(ebetax-mean(ebetax,na.rm=T))
  abetax = (6.703e-02*age+
              -1.529e-01*female+
              -2.393e-02*black+
              5.399e-01*cursmoke+
              -2.988e-03*sysbp+
              8.776e-02*rxbp+
              -2.681e-01*statin+
              4.036e-01*anticoag+
              -0.2574724*oralrx+
              5.888e-01*cvdhist+
              1.659e-01*hgba1c+
              -9.478e-04*totchol+
              -4.378e-03*hdlc+
              3.597e-01*sercreat+
              3.889e-04*uralbcreat)
  mortrisk_10 = 1 - .93^exp(abetax-mean(abetax,na.rm=T))
     lifeexp = (female==0)*(black==0)*(hisp==0)*( # NH WM (non-Hispanic White Male)
      62.1*(age<20)+
        57.3*(age>=20 & age<25)+
        52.6*(age>=25 & age<30)+
        47.9*(age>=30 & age<35)+
        43.3*(age>=35 & age<40)+
        38.7*(age>=40 & age<45)+
        34.2*(age>=45 & age<50)+
        29.8*(age>=50 & age<55)+
        25.7*(age>=55 & age<60)+
        21.7*(age>=60 & age<65)+
        18*(age>=65 & age<70)+
        14.4*(age>=70 & age<75)+
        11.2*(age>=75 & age<80)+
        8.3*(age>=80)
    )+
    (female==1)*(black==0)*(hisp==0)*( # NH WF
      66.7*(age<20)+
        61.8*(age>=20 & age<25)+
        56.9*(age>=25 & age<30)+
        52.1*(age>=30 & age<35)+
        47.3*(age>=35 & age<40)+
        42.5*(age>=40 & age<45)+
        37.9*(age>=45 & age<50)+
        33.3*(age>=50 & age<55)+
        28.9*(age>=55 & age<60)+
        24.7*(age>=60 & age<65)+
        20.5*(age>=65 & age<70)+
        16.6*(age>=70 & age<75)+
        13*(age>=75 & age<80)+
        9.7*(age>=80)
    )+
    (female==0)*(black==1)*(hisp==0)*( # NH BM
      58.3*(age<20)+
        53.6*(age>=20 & age<25)+
        49*(age>=25 & age<30)+
        44.5*(age>=30 & age<35)+
        40*(age>=35 & age<40)+
        35.6*(age>=40 & age<45)+
        31.2*(age>=45 & age<50)+
        27*(age>=50 & age<55)+
        23*(age>=55 & age<60)+
        19.4*(age>=60 & age<65)+
        16.3*(age>=65 & age<70)+
        13.2*(age>=70 & age<75)+
        10.5*(age>=75 & age<80)+
        8.1*(age>=80)  
    )+
    (female==1)*(black==1)*(hisp==0)*( # NH BF
      64.2*(age<20)+
        59.3*(age>=20 & age<25)+
        54.5*(age>=25 & age<30)+
        49.7*(age>=30 & age<35)+
        44.9*(age>=35 & age<40)+
        40.3*(age>=40 & age<45)+
        35.7*(age>=45 & age<50)+
        31.3*(age>=50 & age<55)+
        27.2*(age>=55 & age<60)+
        23.3*(age>=60 & age<65)+
        19.5*(age>=65 & age<70)+
        16*(age>=70 & age<75)+
        12.7*(age>=75 & age<80)+
        9.8*(age>=80)      
    )+
    (female==0)*(black==0)*(hisp==1)*( # HM
      65*(age<20)+
        60.1*(age>=20 & age<25)+
        55.4*(age>=25 & age<30)+
        50.7*(age>=30 & age<35)+
        45.9*(age>=35 & age<40)+
        41.2*(age>=40 & age<45)+
        36.6*(age>=45 & age<50)+
        32.1*(age>=50 & age<55)+
        27.7*(age>=55 & age<60)+
        23.6*(age>=60 & age<65)+
        19.7*(age>=65 & age<70)+
        16*(age>=70 & age<75)+
        12.6*(age>=75 & age<80)+
        9.5*(age>=80)    
    )+
    (female==1)*(black==0)*(hisp==1)*( # HF
      70*(age<20)+
        65*(age>=20 & age<25)+
        60.2*(age>=25 & age<30)+
        55.3*(age>=30 & age<35)+
        50.4*(age>=35 & age<40)+
        45.5*(age>=40 & age<45)+
        40.7*(age>=45 & age<50)+
        36.1*(age>=50 & age<55)+
        31.5*(age>=55 & age<60)+
        27*(age>=60 & age<65)+
        22.8*(age>=65 & age<70)+
        18.7*(age>=70 & age<75)+
        14.8*(age>=75 & age<80)+
        11.2*(age>=80)       
    )
  
  ascvd_rate = -log(1-ascvd_10)/10
  ascvd_life = 1-exp(-ascvd_rate*lifeexp)
  
  dminc_rate = -log(1-dminc_10)/10
  dminc_life = 1-exp(-dminc_rate*lifeexp)
  
  esrdrisk_rate = -log(1-esrdrisk_10)/10
  esrdrisk_life = 1-exp(-esrdrisk_rate*lifeexp)
  
  retinrisk_rate = -log(1-retinrisk_10)/10
  retinrisk_life = 1-exp(-retinrisk_rate*lifeexp)
  
  neurorisk_rate = -log(1-neurorisk_10)/10
  neurorisk_life = 1-exp(-neurorisk_rate*lifeexp)
  
  mortrisk_rate = -log(1-mortrisk_10)/10
  mortrisk_life = 1-exp(-mortrisk_rate*lifeexp)
  popl=length(ascvd_10)
  set.seed(iter*1)
  hr_ascvd = rnorm(popl, 1-0.91, (0.91-0.84)/1.96)/.2
  hr_ascvd[hr_ascvd<0]=0
  set.seed(iter*2)
  hr_dminc = rnorm(popl, 1-0.84, (0.84-0.78)/1.96)/.1
  hr_dminc[hr_dminc<0]=0
  set.seed(iter*3)
  delta_a1c = rnorm(popl, 0.32, (0.32-0.23)/1.96)/.1
  delta_a1c[delta_a1c<0]=0
  # appendix to doi:10.1001/jamainternmed.2014.2894
  # change in microvasc per 1% change in a1c
  hr_renal = 1-1/1.14/.9
  hr_neuro = 1-1/1.19/.9
  hr_retin = 1-1/1.29/.9
  set.seed(iter*4)
  hr_allmort = rnorm(popl, 1-0.83, (0.83-0.78)/1.96)/.2
  hr_allmort[hr_allmort<0]=0
  
  # rx effect 1
  set.seed(iter*5)
  rx1_heidelta = rnorm(popl,0.13,(.13-.09)/1.96)
  rx1_heidelta[rx1_heidelta<0]=0
  rx1_ascvd_10 = (1-rx1_heidelta*hr_ascvd)*ascvd_10
  rx1_ascvd_life = (1-rx1_heidelta*hr_ascvd)*ascvd_life
  rx1_ascvd_life[rx1_ascvd_life>1]=1
  rx1_dminc_10 = (1-rx1_heidelta*hr_dminc)*dminc_10
  rx1_dminc_life = (1-rx1_heidelta*hr_dminc)*dminc_life
  rx1_esrdrisk_10 = (1-rx1_heidelta*delta_a1c*hr_renal)*esrdrisk_10
  rx1_esrdrisk_life = (1-rx1_heidelta*delta_a1c*hr_renal)*esrdrisk_life
  rx1_retinrisk_10 = (1-rx1_heidelta*delta_a1c*hr_retin)*retinrisk_10
  rx1_retinrisk_life = (1-rx1_heidelta*delta_a1c*hr_retin)*retinrisk_life
  rx1_neurorisk_10 = (1-rx1_heidelta*delta_a1c*hr_neuro)*neurorisk_10
  rx1_neurorisk_life = (1-rx1_heidelta*delta_a1c*hr_neuro)* neurorisk_life
  rx1_allmort_10 = (1-rx1_heidelta*hr_allmort)*mortrisk_10
  rx1_allmort_life = (1-rx1_heidelta*hr_allmort)*mortrisk_life
  
  rx1_dminc_rate = -log(1-rx1_dminc_10)/10
  rx1_dminc_life = 1-exp(-rx1_dminc_rate*lifeexp)
  
  # rx effect 2
  set.seed(iter*6)
  rx2_heidelta = rnorm(popl,0.07,(.07-.03)/1.96) 
  rx2_heidelta[rx2_heidelta<0]=0
  rx2_ascvd_10 = (1-rx2_heidelta*hr_ascvd)*ascvd_10
  rx2_ascvd_life = (1-rx2_heidelta*hr_ascvd)*ascvd_life
  rx2_ascvd_life[rx2_ascvd_life>1]=1
  rx2_dminc_10 = (1-rx2_heidelta*hr_dminc)*dminc_10
  rx2_dminc_life = (1-rx2_heidelta*hr_dminc)*dminc_life
  rx2_esrdrisk_10 = (1-rx2_heidelta*delta_a1c*hr_renal)*esrdrisk_10
  rx2_esrdrisk_life = (1-rx2_heidelta*delta_a1c*hr_renal)*esrdrisk_life
  rx2_retinrisk_10 = (1-rx2_heidelta*delta_a1c*hr_retin)*retinrisk_10
  rx2_retinrisk_life = (1-rx2_heidelta*delta_a1c*hr_retin)*retinrisk_life
  rx2_neurorisk_10 = (1-rx2_heidelta*delta_a1c*hr_neuro)*neurorisk_10
  rx2_neurorisk_life = (1-rx2_heidelta*delta_a1c*hr_neuro)* neurorisk_life
  rx2_allmort_10 = (1-rx2_heidelta*hr_allmort)*mortrisk_10
  rx2_allmort_life = (1-rx2_heidelta*hr_allmort)*mortrisk_life
  
  rx2_dminc_rate = -log(1-rx2_dminc_10)/10
  rx2_dminc_life = 1-exp(-rx2_dminc_rate*lifeexp)
  
  set.seed(iter*7)
  ascvd_events_10 = rbinom(popl, 1, ascvd_10)
  set.seed(iter*8)
  ascvd_events_life = rbinom(popl, 1, ascvd_life)
  set.seed(iter*9)
  dminc_events_10 = rbinom(popl, 1, dminc_10)*(dm==0)
  set.seed(iter*10)
  dminc_events_life = rbinom(popl, 1, dminc_life)*(dm==0)
  set.seed(iter*11)
  dminter_10 = rbinom(popl,1,1-exp(-dminc_rate*10/2))
  set.seed(iter*12)
  dminter_life = rbinom(popl,1,1-exp(-dminc_rate*lifeexp/2))
  set.seed(iter*13)
  esrd_events_10 = rbinom(popl, 1, esrdrisk_10)*(dm==1 | dminter_10==1)
  set.seed(iter*14)
  esrd_events_life = rbinom(popl, 1, esrdrisk_life)*(dm==1 | dminter_life==1)
  set.seed(iter*15)
  retin_events_10 = rbinom(popl, 1, retinrisk_10)*(dm==1 | dminter_10==1)
  set.seed(iter*16)
  retin_events_life = rbinom(popl, 1, retinrisk_life)*(dm==1 | dminter_life==1)
  set.seed(iter*17)
  neuro_events_10 = rbinom(popl, 1, neurorisk_10)*(dm==1 | dminter_10==1)
  set.seed(iter*18)
  neuro_events_life = rbinom(popl, 1, neurorisk_life)*(dm==1 | dminter_life==1)
  set.seed(iter*19)
  mort_events_10 = rbinom(popl, 1, mortrisk_10)
  set.seed(iter*20)
  mort_events_life = rbinom(popl, 1, mortrisk_life)
  
   set.seed(iter*21)
  d_ascvd = rnorm(popl,.28,(.28-.06)/1.96)
  d_ascvd[d_ascvd<0]=0
  d_dminc = rnorm(popl,.01,.01/1.96)
  d_dminc[d_dminc<0]=0
  d_esrd = rnorm(popl,.57,(.57-.4)/1.96)
  d_esrd[d_esrd<0]=0
  d_retin = rnorm(popl,.2,(.2-.1)/1.96)
  d_retin[d_retin<0]=0
  d_neuro = rnorm(popl,.1,(.1-.05)/1.96)
  d_neuro[d_neuro<0]=0
  
  # QALYS lost
  q_ascvd_10 = d_ascvd*10/2*ascvd_events_10*(mort_events_10==0)+
    1*10/2*ascvd_events_10*(mort_events_10==1)
  q_ascvd_life = d_ascvd*lifeexp/2*ascvd_events_life*(mort_events_life==0)+
    1*lifeexp/2*ascvd_events_life*(mort_events_life==1)
  q_dminc_10 = d_dminc*10/2*dminc_events_10*(mort_events_10==0)+
    1*10/2*dminc_events_10*(mort_events_10==1)
  q_dminc_life = d_dminc*lifeexp/2*dminc_events_life*(mort_events_life==0)+
    1*lifeexp/2*dminc_events_life*(mort_events_life==1)
  q_esrd_10 = d_esrd*10/2*esrd_events_10*(mort_events_10==0)+
    1*10/2*esrd_events_10*(mort_events_10==1)
  q_esrd_life = d_esrd*lifeexp/2*esrd_events_life*(mort_events_life==0)+
    1*lifeexp/2*esrd_events_life*(mort_events_life==1)
  q_retin_10 = d_retin*10/2*retin_events_10*(mort_events_10==0)+
    1*10/2*retin_events_10*(mort_events_10==1)
  q_retin_life = d_retin*lifeexp/2*retin_events_life*(mort_events_life==0)+
    1*lifeexp/2*retin_events_life*(mort_events_life==1)
  q_neuro_10 = d_neuro*10/2*neuro_events_10*(mort_events_10==0)+
    1*10/2*neuro_events_10*(mort_events_10==1)
  q_neuro_life = d_neuro*lifeexp/2*neuro_events_life*(mort_events_life==0)+
    1*lifeexp/2*neuro_events_life*(mort_events_life==1)
  q_allmort_10 = 10/2*(ascvd_events_10==0)*(dminc_events_10==0)*(esrd_events_10==0)*(retin_events_10==0)*(neuro_events_10==0)*(mort_events_10==1)
  q_allmort_life = lifeexp/2*(ascvd_events_life==0)*(dminc_events_life==0)*(esrd_events_life==0)*(retin_events_life==0)*(neuro_events_life==0)*(mort_events_life==1)
  
  library(data.table)
  qalys_10 = data.table(q_ascvd_10, q_dminc_10, q_esrd_10, q_retin_10, q_neuro_10, q_allmort_10)/(1.03^5)
  qalys_life = data.table(q_ascvd_life, q_dminc_life, q_esrd_life, q_retin_life, q_neuro_life, q_allmort_life)/(1.03^(lifeexp/2))
  
  
  # QALYS lost under treatment 1
  
  set.seed(iter*7)
  rx1_ascvd_events_10 = rbinom(popl, 1, rx1_ascvd_10)
  set.seed(iter*8)
  rx1_ascvd_events_life = rbinom(popl, 1, rx1_ascvd_life)
  set.seed(iter*9)
  rx1_dminc_events_10 = rbinom(popl, 1, rx1_dminc_10)*(dm==0)
  set.seed(iter*10)
  rx1_dminc_events_life = rbinom(popl, 1, rx1_dminc_life)*(dm==0)
  set.seed(iter*11)
  rx1_dminter_10 = rbinom(popl,1,1-exp(-rx1_dminc_rate*10/2))
  set.seed(iter*12)
  rx1_dminter_life = rbinom(popl,1,1-exp(-rx1_dminc_rate*lifeexp/2))
  set.seed(iter*13)
  rx1_esrd_events_10 = rbinom(popl, 1, rx1_esrdrisk_10)*(dm==1 | rx1_dminter_10==1)
  set.seed(iter*14)
  rx1_esrd_events_life = rbinom(popl, 1, rx1_esrdrisk_life)*(dm==1 | rx1_dminter_life==1)
  set.seed(iter*15)
  rx1_retin_events_10 = rbinom(popl, 1, rx1_retinrisk_10)*(dm==1 | rx1_dminter_10==1)
  set.seed(iter*16)
  rx1_retin_events_life = rbinom(popl, 1, rx1_retinrisk_life)*(dm==1 | rx1_dminter_life==1)
  set.seed(iter*17)
  rx1_neuro_events_10 = rbinom(popl, 1, rx1_neurorisk_10)*(dm==1 | rx1_dminter_10==1)
  set.seed(iter*18)
  rx1_neuro_events_life = rbinom(popl, 1, rx1_neurorisk_life)*(dm==1 | rx1_dminter_life==1)
  set.seed(iter*19)
  rx1_mort_events_10 = rbinom(popl, 1, rx1_allmort_10)
  set.seed(iter*20)
  rx1_mort_events_life = rbinom(popl, 1, rx1_allmort_life)
  
  
  # QALYS lost
  rx1_q_ascvd_10 = d_ascvd*10/2*rx1_ascvd_events_10*(rx1_mort_events_10==0)+
    1*10/2*rx1_ascvd_events_10*(rx1_mort_events_10==1)
  rx1_q_ascvd_life = d_ascvd*lifeexp/2*rx1_ascvd_events_life*(rx1_mort_events_life==0)+
    1*lifeexp/2*rx1_ascvd_events_life*(rx1_mort_events_life==1)
  rx1_q_dminc_10 = d_dminc*10/2*rx1_dminc_events_10*(rx1_mort_events_10==0)+
    1*10/2*rx1_dminc_events_10*(rx1_mort_events_10==1)
  rx1_q_dminc_life = d_dminc*lifeexp/2*rx1_dminc_events_life*(rx1_mort_events_life==0)+
    1*lifeexp/2*rx1_dminc_events_life*(rx1_mort_events_life==1)
  rx1_q_esrd_10 = d_esrd*10/2*rx1_esrd_events_10*(rx1_mort_events_10==0)+
    1*10/2*rx1_esrd_events_10*(rx1_mort_events_10==1)
  rx1_q_esrd_life = d_esrd*lifeexp/2*rx1_esrd_events_life*(rx1_mort_events_life==0)+
    1*lifeexp/2*rx1_esrd_events_life*(rx1_mort_events_life==1)
  rx1_q_retin_10 = d_retin*10/2*rx1_retin_events_10*(rx1_mort_events_10==0)+
    1*10/2*rx1_retin_events_10*(rx1_mort_events_10==1)
  rx1_q_retin_life = d_retin*lifeexp/2*rx1_retin_events_life*(rx1_mort_events_life==0)+
    1*lifeexp/2*rx1_retin_events_life*(rx1_mort_events_life==1)
  rx1_q_neuro_10 = d_neuro*10/2*rx1_neuro_events_10*(rx1_mort_events_10==0)+
    1*10/2*rx1_neuro_events_10*(rx1_mort_events_10==1)
  rx1_q_neuro_life = d_neuro*lifeexp/2*rx1_neuro_events_life*(rx1_mort_events_life==0)+
    1*lifeexp/2*rx1_neuro_events_life*(rx1_mort_events_life==1)
  rx1_q_allmort_10 = 10/2*(rx1_ascvd_events_10==0)*(rx1_dminc_events_10==0)*(rx1_esrd_events_10==0)*(rx1_retin_events_10==0)*(rx1_neuro_events_10==0)*(rx1_mort_events_10==1)
  rx1_q_allmort_life = lifeexp/2*(rx1_ascvd_events_life==0)*(rx1_dminc_events_life==0)*(rx1_esrd_events_life==0)*(rx1_retin_events_life==0)*(rx1_neuro_events_life==0)*(rx1_mort_events_life==1)
  
  rx1_qalys_10 = data.table(rx1_q_ascvd_10, rx1_q_dminc_10, rx1_q_esrd_10, rx1_q_retin_10, rx1_q_neuro_10, rx1_q_allmort_10)/(1.03^5)
  rx1_qalys_life = data.table(rx1_q_ascvd_life, rx1_q_dminc_life, rx1_q_esrd_life, rx1_q_retin_life, rx1_q_neuro_life, rx1_q_allmort_life)/(1.03^(lifeexp/2))
  
  
  
  # QALYS lost under treatment 2
  
  set.seed(iter*7)
  rx2_ascvd_events_10 = rbinom(popl, 1, rx2_ascvd_10)
  set.seed(iter*8)
  rx2_ascvd_events_life = rbinom(popl, 1, rx2_ascvd_life)
  set.seed(iter*9)
  rx2_dminc_events_10 = rbinom(popl, 1, rx2_dminc_10)*(dm==0)
  set.seed(iter*10)
  rx2_dminc_events_life = rbinom(popl, 1, rx2_dminc_life)*(dm==0)
  set.seed(iter*11)
  rx2_dminter_10 = rbinom(popl,1,1-exp(-rx2_dminc_rate*10/2))
  set.seed(iter*12)
  rx2_dminter_life = rbinom(popl,1,1-exp(-rx2_dminc_rate*lifeexp/2))
  set.seed(iter*13)
  rx2_esrd_events_10 = rbinom(popl, 1, rx2_esrdrisk_10)*(dm==1 | rx2_dminter_10==1)
  set.seed(iter*14)
  rx2_esrd_events_life = rbinom(popl, 1, rx2_esrdrisk_life)*(dm==1 | rx2_dminter_life==1)
  set.seed(iter*15)
  rx2_retin_events_10 = rbinom(popl, 1, rx2_retinrisk_10)*(dm==1 | rx2_dminter_10==1)
  set.seed(iter*16)
  rx2_retin_events_life = rbinom(popl, 1, rx2_retinrisk_life)*(dm==1 | rx2_dminter_life==1)
  set.seed(iter*17)
  rx2_neuro_events_10 = rbinom(popl, 1, rx2_neurorisk_10)*(dm==1 | rx2_dminter_10==1)
  set.seed(iter*18)
  rx2_neuro_events_life = rbinom(popl, 1, rx2_neurorisk_life)*(dm==1 | rx2_dminter_life==1)
  set.seed(iter*19)
  rx2_mort_events_10 = rbinom(popl, 1, rx2_allmort_10)
  set.seed(iter*20)
  rx2_mort_events_life = rbinom(popl, 1, rx2_allmort_life)
  
  
  # QALYS lost
  rx2_q_ascvd_10 = d_ascvd*10/2*rx2_ascvd_events_10*(rx2_mort_events_10==0)+
    1*10/2*rx2_ascvd_events_10*(rx2_mort_events_10==1)
  rx2_q_ascvd_life = d_ascvd*lifeexp/2*rx2_ascvd_events_life*(rx2_mort_events_life==0)+
    1*lifeexp/2*rx2_ascvd_events_life*(rx2_mort_events_life==1)
  rx2_q_dminc_10 = d_dminc*10/2*rx2_dminc_events_10*(rx2_mort_events_10==0)+
    1*10/2*rx2_dminc_events_10*(rx2_mort_events_10==1)
  rx2_q_dminc_life = d_dminc*lifeexp/2*rx2_dminc_events_life*(rx2_mort_events_life==0)+
    1*lifeexp/2*rx2_dminc_events_life*(rx2_mort_events_life==1)
  rx2_q_esrd_10 = d_esrd*10/2*rx2_esrd_events_10*(rx2_mort_events_10==0)+
    1*10/2*rx2_esrd_events_10*(rx2_mort_events_10==1)
  rx2_q_esrd_life = d_esrd*lifeexp/2*rx2_esrd_events_life*(rx2_mort_events_life==0)+
    1*lifeexp/2*rx2_esrd_events_life*(rx2_mort_events_life==1)
  rx2_q_retin_10 = d_retin*10/2*rx2_retin_events_10*(rx2_mort_events_10==0)+
    1*10/2*rx2_retin_events_10*(rx2_mort_events_10==1)
  rx2_q_retin_life = d_retin*lifeexp/2*rx2_retin_events_life*(rx2_mort_events_life==0)+
    1*lifeexp/2*rx2_retin_events_life*(rx2_mort_events_life==1)
  rx2_q_neuro_10 = d_neuro*10/2*rx2_neuro_events_10*(rx2_mort_events_10==0)+
    1*10/2*rx2_neuro_events_10*(rx2_mort_events_10==1)
  rx2_q_neuro_life = d_neuro*lifeexp/2*rx2_neuro_events_life*(rx2_mort_events_life==0)+
    1*lifeexp/2*rx2_neuro_events_life*(rx2_mort_events_life==1)
  rx2_q_allmort_10 = 10/2*(rx2_ascvd_events_10==0)*(rx2_dminc_events_10==0)*(rx2_esrd_events_10==0)*(rx2_retin_events_10==0)*(rx2_neuro_events_10==0)*(rx2_mort_events_10==1)
  rx2_q_allmort_life = lifeexp/2*(rx2_ascvd_events_life==0)*(rx2_dminc_events_life==0)*(rx2_esrd_events_life==0)*(rx2_retin_events_life==0)*(rx2_neuro_events_life==0)*(rx2_mort_events_life==1)
  
  rx2_qalys_10 = data.table(rx2_q_ascvd_10, rx2_q_dminc_10, rx2_q_esrd_10, rx2_q_retin_10, rx2_q_neuro_10, rx2_q_allmort_10)/(1.03^5)
  rx2_qalys_life = data.table(rx2_q_ascvd_life, rx2_q_dminc_life, rx2_q_esrd_life, rx2_q_retin_life, rx2_q_neuro_life, rx2_q_allmort_life)/(1.03^(lifeexp/2))
   
  set.seed(iter*21)
  c_ascvd = rnorm(popl,61614,(61614-61095)/1.96)
  c_ascvd[c_ascvd<0]=0
  set.seed(iter*22)
  c_dminc = rnorm(popl,22526,(22526-21919)/1.96)
  c_dminc[c_dminc<0]=0
  set.seed(iter*23)
  c_esrd = rnorm(popl,293490,(293490-290584)/1.96)
  c_esrd[c_esrd<0]=0
  set.seed(iter*24)
  c_retin = rnorm(popl,26174,(26174-25578)/1.96)
  c_retin[c_retin<0]=0
  set.seed(iter*25)
  c_neuro = rnorm(popl,48432,(48432-46470)/1.96)
  c_neuro[c_neuro<0]=0
  
  
  # healthcare costs at baseline
  c_ascvd_10 = c_ascvd*ascvd_events_10
  c_ascvd_life = c_ascvd*ascvd_events_life
  c_dminc_10 = c_dminc*dminc_events_10
  c_dminc_life = c_dminc*dminc_events_life
  c_esrd_10 = c_esrd*esrd_events_10
  c_esrd_life = c_esrd*esrd_events_life
  c_retin_10 = c_retin*retin_events_10
  c_retin_life = c_retin*retin_events_life
  c_neuro_10 = c_neuro*neuro_events_10
  c_neuro_life = c_neuro*neuro_events_life
  
  costs_10 = data.table(c_ascvd_10, c_dminc_10, c_esrd_10, c_retin_10, c_neuro_10)/(1.03^5)
  costs_life = data.table(c_ascvd_life, c_dminc_life, c_esrd_life, c_retin_life, c_neuro_life)/(1.03^(lifeexp/2))
  
  # healthcare costs under rx 1
  rx1_c_ascvd_10 = c_ascvd*rx1_ascvd_events_10
  rx1_c_ascvd_life = c_ascvd*rx1_ascvd_events_life
  rx1_c_dminc_10 = c_dminc*rx1_dminc_events_10
  rx1_c_dminc_life = c_dminc*rx1_dminc_events_life
  rx1_c_esrd_10 = c_esrd*rx1_esrd_events_10
  rx1_c_esrd_life = c_esrd*rx1_esrd_events_life
  rx1_c_retin_10 = c_retin*rx1_retin_events_10
  rx1_c_retin_life = c_retin*rx1_retin_events_life
  rx1_c_neuro_10 = c_neuro*rx1_neuro_events_10
  rx1_c_neuro_life = c_neuro*rx1_neuro_events_life
  
  rx1_costs_10 = data.table(rx1_c_ascvd_10, rx1_c_dminc_10, rx1_c_esrd_10, rx1_c_retin_10, rx1_c_neuro_10)/(1.03^5)
  rx1_costs_life = data.table(rx1_c_ascvd_life, rx1_c_dminc_life, rx1_c_esrd_life, rx1_c_retin_life, rx1_c_neuro_life)/(1.03^(lifeexp/2))
  
  # healthcare costs under rx 2
  rx2_c_ascvd_10 = c_ascvd*rx2_ascvd_events_10
  rx2_c_ascvd_life = c_ascvd*rx2_ascvd_events_life
  rx2_c_dminc_10 = c_dminc*rx2_dminc_events_10
  rx2_c_dminc_life = c_dminc*rx2_dminc_events_life
  rx2_c_esrd_10 = c_esrd*rx2_esrd_events_10
  rx2_c_esrd_life = c_esrd*rx2_esrd_events_life
  rx2_c_retin_10 = c_retin*rx2_retin_events_10
  rx2_c_retin_life = c_retin*rx2_retin_events_life
  rx2_c_neuro_10 = c_neuro*rx2_neuro_events_10
  rx2_c_neuro_life = c_neuro*rx2_neuro_events_life
  
  rx2_costs_10 = data.table(rx2_c_ascvd_10, rx2_c_dminc_10, rx2_c_esrd_10, rx2_c_retin_10, rx2_c_neuro_10)/(1.03^5)
  rx2_costs_life = data.table(rx2_c_ascvd_life, rx2_c_dminc_life, rx2_c_esrd_life, rx2_c_retin_life, rx2_c_neuro_life)/(1.03^(lifeexp/2))
  
  
  set.seed(iter*21)
  p_ascvd = rnorm(popl,1085,(1085-521)/1.96)
  p_ascvd[p_ascvd<0]=0
  p_dminc = 0
  set.seed(iter*23)
  p_esrd = rnorm(popl,5811,(5811-2571)/1.96)
  p_esrd[p_esrd<0]=0
  p_retin = p_esrd
  p_neuro = p_esrd
  
  
  # productivity costs at baseline
  p_ascvd_10 = p_ascvd*10/2*ascvd_events_10
  p_ascvd_life = p_ascvd*lifeexp/2*ascvd_events_life
  p_dminc_10 = p_dminc*10/2*dminc_events_10
  p_dminc_life = p_dminc*lifeexp/2*dminc_events_life
  p_esrd_10 = p_esrd*10/2*esrd_events_10
  p_esrd_life = p_esrd*lifeexp/2*esrd_events_life
  p_retin_10 = p_retin*10/2*retin_events_10
  p_retin_life = p_retin*lifeexp/2*retin_events_life
  p_neuro_10 = p_neuro*10/2*neuro_events_10
  p_neuro_life = p_neuro*lifeexp/2*neuro_events_life
  
  prod_10 = data.table(p_ascvd_10, p_dminc_10, p_esrd_10, p_retin_10, p_neuro_10)
  prod_life = data.table(p_ascvd_life, p_dminc_life, p_esrd_life, p_retin_life, p_neuro_life)
  
  # productivity costs under rx 1
  rx1_p_ascvd_10 = p_ascvd*10/2*rx1_ascvd_events_10
  rx1_p_ascvd_life = p_ascvd*lifeexp/2*rx1_ascvd_events_life
  rx1_p_dminc_10 = p_dminc*10/2*rx1_dminc_events_10
  rx1_p_dminc_life = p_dminc*lifeexp/2*rx1_dminc_events_life
  rx1_p_esrd_10 = p_esrd*10/2*rx1_esrd_events_10
  rx1_p_esrd_life = p_esrd*lifeexp/2*rx1_esrd_events_life
  rx1_p_retin_10 = p_retin*10/2*rx1_retin_events_10
  rx1_p_retin_life = p_retin*lifeexp/2*rx1_retin_events_life
  rx1_p_neuro_10 = p_neuro*10/2*rx1_neuro_events_10
  rx1_p_neuro_life = p_neuro*lifeexp/2*rx1_neuro_events_life
  
  rx1_prod_10 = data.table(rx1_p_ascvd_10, rx1_p_dminc_10, rx1_p_esrd_10, rx1_p_retin_10, rx1_p_neuro_10)
  rx1_prod_life = data.table(rx1_p_ascvd_life, rx1_p_dminc_life, rx1_p_esrd_life, rx1_p_retin_life, rx1_p_neuro_life)
  
  # productivity costs under rx 2
  rx2_p_ascvd_10 = p_ascvd*10/2*rx2_ascvd_events_10
  rx2_p_ascvd_life = p_ascvd*lifeexp/2*rx2_ascvd_events_life
  rx2_p_dminc_10 = p_dminc*10/2*rx2_dminc_events_10
  rx2_p_dminc_life = p_dminc*lifeexp/2*rx2_dminc_events_life
  rx2_p_esrd_10 = p_esrd*10/2*rx2_esrd_events_10
  rx2_p_esrd_life = p_esrd*lifeexp/2*rx2_esrd_events_life
  rx2_p_retin_10 = p_retin*10/2*rx2_retin_events_10
  rx2_p_retin_life = p_retin*lifeexp/2*rx2_retin_events_life
  rx2_p_neuro_10 = p_neuro*10/2*rx2_neuro_events_10
  rx2_p_neuro_life = p_neuro*lifeexp/2*rx2_neuro_events_life
  
  rx2_prod_10 = data.table(rx2_p_ascvd_10, rx2_p_dminc_10, rx2_p_esrd_10, rx2_p_retin_10, rx2_p_neuro_10)
  rx2_prod_life = data.table(rx2_p_ascvd_life, rx2_p_dminc_life, rx2_p_esrd_life, rx2_p_retin_life, rx2_p_neuro_life)
  
  
  
  
  # intervention costs, with placeholder 15% overhead cost
  overhead = 0.15
  rx1_intcosts_10 = sum(300*10*(1+overhead)*(rx1_mort_events_10==0)+300*10/2*(1+overhead)*(rx1_mort_events_10==1))
  rx1_intcosts_life = sum(300*lifeexp*(1+overhead)*(rx1_mort_events_life==0)+300*lifeexp/2*(1+overhead)*(rx1_mort_events_life==1))
  
  rx2_intcosts_10 = sum(300*10*(1+overhead)*(rx2_mort_events_10==0)+300*10/2*(1+overhead)*(rx2_mort_events_10==1))
  rx2_intcosts_life = sum(300*lifeexp*(1+overhead)*(rx2_mort_events_life==0)+300*lifeexp/2*(1+overhead)*(rx2_mort_events_life==1))
  
  
  # economic benefits from intervention, placeholder $50 per person per year net increase benefit to farms
  econben = 50
  rx1_econben_10 =sum(econben*10*(rx1_mort_events_10==0)+econben*10/2*(rx1_mort_events_10==1))
  rx1_econben_life =sum(econben*lifeexp*(rx1_mort_events_life==0)+econben*lifeexp/2*(rx1_mort_events_life==1))
  
  rx2_econben_10 =sum(econben*10*(rx2_mort_events_10==0)+econben*10/2*(rx2_mort_events_10==1))
  rx2_econben_life =sum(econben*lifeexp*(rx2_mort_events_life==0)+econben*lifeexp/2*(rx2_mort_events_life==1))
  
  library(fBasics)
  library(matrixStats)
  out_baseline_10 = rowQuantiles(as.matrix(rbind(ascvd_10,dminc_10,esrdrisk_10,neurorisk_10,retinrisk_10)),probs=c(.025,.5,.975))
  out_rx1_10 = rowQuantiles(as.matrix(rbind(rx1_ascvd_10,rx1_dminc_10,rx1_esrdrisk_10,rx1_neurorisk_10,rx1_retinrisk_10)),probs=c(.025,.5,.975))
  out_rx2_10 = rowQuantiles(as.matrix(rbind(rx2_ascvd_10,rx2_dminc_10,rx2_esrdrisk_10,rx2_neurorisk_10,rx2_retinrisk_10)),probs=c(.025,.5,.975))

  out_baseline_life = rowQuantiles(as.matrix(rbind(ascvd_life,dminc_life,esrdrisk_life,neurorisk_life,retinrisk_life)),probs=c(.025,.5,.975))
  out_rx1_life = rowQuantiles(as.matrix(rbind(rx1_ascvd_life,rx1_dminc_life,rx1_esrdrisk_life,rx1_neurorisk_life,rx1_retinrisk_life)),probs=c(.025,.5,.975))
  out_rx2_life = rowQuantiles(as.matrix(rbind(rx2_ascvd_life,rx2_dminc_life,rx2_esrdrisk_life,rx2_neurorisk_life,rx2_retinrisk_life)),probs=c(.025,.5,.975))
    
  events_baseline_10 = rowSums(as.matrix(rbind(ascvd_events_10,dminc_events_10,esrd_events_10,neuro_events_10,retin_events_10)))/popl*10000
  events_rx1_10 = rowSums(as.matrix(rbind(rx1_ascvd_events_10,rx1_dminc_events_10,rx1_esrd_events_10,rx1_neuro_events_10,rx1_retin_events_10)))/popl*10000
  events_rx2_10 = rowSums(as.matrix(rbind(rx2_ascvd_events_10,rx2_dminc_events_10,rx2_esrd_events_10,rx2_neuro_events_10,rx2_retin_events_10)))/popl*10000
  
  events_baseline_life = rowSums(as.matrix(rbind(ascvd_events_life,dminc_events_life,esrd_events_life,neuro_events_life,retin_events_life)))/popl*10000
  events_rx1_life = rowSums(as.matrix(rbind(rx1_ascvd_events_life,rx1_dminc_events_life,rx1_esrd_events_life,rx1_neuro_events_life,rx1_retin_events_life)))/popl*10000
  events_rx2_life = rowSums(as.matrix(rbind(rx2_ascvd_events_life,rx2_dminc_events_life,rx2_esrd_events_life,rx2_neuro_events_life,rx2_retin_events_life)))/popl*10000
  
  
  # return(list(out_base_10 = out_baseline_10,
  #             out_rx1_10 = out_rx1_10,
  #             out_rx2_10 = out_rx2_10,
  #             out_base_life = out_baseline_life,
  #             out_rx1_life = out_rx1_life,
  #             out_rx2_life = out_rx2_life))
  
  
  
  # eventstab = data.table(events_baseline_10,events_rx1_10,events_rx2_10,events_baseline_life,events_rx1_life,events_rx2_life)
  # rx1_deltatab_10 = eventstab[,2]-eventstab[,1]
  # rx1_deltatab_life = eventstab[,5]-eventstab[,4]
  # rx2_deltatab_10 = eventstab[,3]-eventstab[,1]
  # rx2_deltatab_life = eventstab[,6]-eventstab[,4]
  # rx1_deltatab_10[rx1_deltatab_10>0]=mean(rx1_deltatab_10)
  # rx1_deltatab_life[rx1_deltatab_life>0]=mean(rx1_deltatab_life)
  # rx2_deltatab_10[rx2_deltatab_10>0]=mean(rx2_deltatab_10)
  # rx2_deltatab_life[rx2_deltatab_life>0]=mean(rx2_deltatab_life)
  # 
  # output = cbind(rx1_deltatab_10,rx1_deltatab_life,rx2_deltatab_10,rx2_deltatab_life)
  # output2 = matrix(unlist(output),ncol=4,byrow=F)
  # return(output2)
  
  out_q_10 = colSums(qalys_10)/popl*10000
  rx1_out_q_10 = colSums(rx1_qalys_10)/popl*10000
  rx2_out_q_10 = colSums(rx2_qalys_10)/popl*10000
  
  out_q_life = colSums(qalys_life)/popl*10000
  rx1_out_q_life = colSums(rx1_qalys_life)/popl*10000
  rx2_out_q_life = colSums(rx2_qalys_life)/popl*10000
  
  output = rbind(out_q_10,rx1_out_q_10,rx2_out_q_10,out_q_life,rx1_out_q_life,rx2_out_q_life)
  output2 = matrix(unlist(output),ncol=6,byrow=F)
  return(output2)
  
}
