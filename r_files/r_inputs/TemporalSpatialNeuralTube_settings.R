rm(list=ls())

## choose the directory where your inputs and outputs will be
workingdir="/Users/delasj/Documents/BriscoeLab/project_TempoSpace_reproduce/"

## Colors and shapes
sorted.DayGate <- c("D5_p1","D5_p2","D5_pM",
                    "D7_p1","D7_p2","D7_pM",
                    "D9_p1","D9_p2","D9_pM",
                    "D11_p1","D11_p2","D11_pM")
sorted.day <- c("D5","D7","D9","D11")
sorted.dayNfia <- c("D5_NFIAn_WT","D5_NFIAn_MUT",
                    "D7_NFIAn_WT","D7_NFIAn_MUT",
                    "D9_NFIAn_WT","D9_NFIAp_WT","D9_NFIAn_MUT",
                    "D11_NFIAp_WT","D11_NFIAn_MUT")


sorted.sample <- c("WT_D5_p1_NFIAn","WT_D5_p2_NFIAn","WT_D5_pM_NFIAn","WT_D7_p1_NFIAn","WT_D7_p2_NFIAn","WT_D7_pM_NFIAn","WT_D9_p1_NFIAn","WT_D9_p2_NFIAn","WT_D9_p2_NFIAp","WT_D9_pM_NFIAn","WT_D9_pM_NFIAp","WT_D11_p1_NFIAp","WT_D11_p2_NFIAp","WT_D11_pM_NFIAp")

sorted_DayGate <- c("D5_p1","D5_p2","D5_pM",
                    "D7_p1","D7_p2","D7_pM",
                    "D9_p1","D9_p2","D9_pM",
                    "D11_p1","D11_p2","D11_pM")
sorted_WTDayGate <- paste("WT",sorted_DayGate, sep="_")

gates <- c("p1","p2","pM")

colorIZ <- c("#abdff4","#f1df9a","#f19aac",
             "#55bee8","#e6c444","#e64466",
             "#1a91c1","#c19e1a","#c11a3d",
             "#0e506b","#6b570e","#7c1127")

color_days <- c("#fadede","#f3aaaa","#e96666","#cf1e1e")  

colorgates <- c("#55bee8","#e6c444","#e64466")

color_screen <- c("#e2e2e2","#a5143b","#00b184","#FF8370")
color_circlesscreen <- c("#a7a7a7","black","black","black")

Nr6a1colors <- c("#a7a7a7","#9f0000")

shapes4_manual = c(18,15,16,17) # these are block
shapes5_manual = c(25,21,22,23,24) # these are filled
shapes4_fill_manual = c(23,21,22,24)

sorted.sample.wReps.all <- c("WT_D5_p1_NFIAn_R1","WT_D5_p1_NFIAn_R2",
                             "WT_D5_p2_NFIAn_R1","WT_D5_p2_NFIAn_R2",
                             "WT_D5_pM_NFIAn_R1","WT_D5_pM_NFIAn_R2",
                             "WT_D7_p1_NFIAn_R1","WT_D7_p1_NFIAn_R2","WT_D7_p1_NFIAn_R3",
                             "WT_D7_p2_NFIAn_R1","WT_D7_p2_NFIAn_R2","WT_D7_p2_NFIAn_R3",
                             "WT_D7_pM_NFIAn_R1","WT_D7_pM_NFIAn_R2","WT_D7_pM_NFIAn_R3",
                             "WT_D9_p1_NFIAn_R2","WT_D9_p1_NFIAn_R3","WT_D9_p1_NFIAn_R1",
                             "WT_D9_p2_NFIAn_R1","WT_D9_p2_NFIAn_R2","WT_D9_p2_NFIAn_R3",
                             "WT_D9_p2_NFIAp_R2","WT_D9_p2_NFIAp_R3","WT_D9_p2_NFIAp_R1",
                             "WT_D9_pM_NFIAp_R1","WT_D9_pM_NFIAp_R2","WT_D9_pM_NFIAp_R3",
                             "WT_D9_pM_NFIAn_R2","WT_D9_pM_NFIAn_R1",
                             "WT_D11_p1_NFIAp_R1","WT_D11_p1_NFIAp_R2","WT_D11_p1_NFIAp_R3",
                             "WT_D11_p2_NFIAp_R1","WT_D11_p2_NFIAp_R2","WT_D11_p2_NFIAp_R3",
                             "WT_D11_pM_NFIAp_R1","WT_D11_pM_NFIAp_R2","WT_D11_pM_NFIAp_R3",
                             "MUT_D5_p1_NFIAn_R1","MUT_D5_p1_NFIAn_R2","MUT_D5_p1_NFIAn_R3",
                             "MUT_D5_p2_NFIAn_R1","MUT_D5_p2_NFIAn_R2","MUT_D5_p2_NFIAn_R3",
                             "MUT_D5_pM_NFIAn_R1","MUT_D5_pM_NFIAn_R2","MUT_D5_pM_NFIAn_R3",
                             "MUT_D7_p1_NFIAn_R1","MUT_D7_p1_NFIAn_R2","MUT_D7_p1_NFIAn_R3",
                             "MUT_D7_p2_NFIAn_R1","MUT_D7_p2_NFIAn_R2","MUT_D7_p2_NFIAn_R3",
                             "MUT_D7_pM_NFIAn_R1","MUT_D7_pM_NFIAn_R2","MUT_D7_pM_NFIAn_R3",
                             "MUT_D9_p1_NFIAn_R2","MUT_D9_p1_NFIAn_R1",
                             "MUT_D9_p2_NFIAn_R1","MUT_D9_p2_NFIAn_R2",
                             "MUT_D9_pM_NFIAn_R2","MUT_D9_pM_NFIAn_R1",
                             "MUT_D11_p1_NFIAn_R1","MUT_D11_p1_NFIAn_R2","MUT_D11_p1_NFIAn_R3",
                             "MUT_D11_p2_NFIAn_R1","MUT_D11_p2_NFIAn_R2","MUT_D11_p2_NFIAn_R3",
                             "MUT_D11_pM_NFIAn_R1","MUT_D11_pM_NFIAn_R2","MUT_D11_pM_NFIAn_R3")

sorted.sample.wReps.bycelltype <-c("WT_D5_p1_NFIAn_R1","WT_D5_p1_NFIAn_R2","WT_D5_p1_NFIAn_R3",
                                   "WT_D7_p1_NFIAn_R1","WT_D7_p1_NFIAn_R2","WT_D7_p1_NFIAn_R3",
                                   "WT_D9_p1_NFIAn_R1","WT_D9_p1_NFIAn_R2","WT_D9_p1_NFIAn_R3",
                                   "WT_D11_p1_NFIAp_R1","WT_D11_p1_NFIAp_R2","WT_D11_p1_NFIAp_R3",
                                   "WT_D5_p2_NFIAn_R1","WT_D5_p2_NFIAn_R2","WT_D5_p2_NFIAn_R3",
                                   "WT_D7_p2_NFIAn_R1","WT_D7_p2_NFIAn_R2","WT_D7_p2_NFIAn_R3",
                                   "WT_D9_p2_NFIAn_R1","WT_D9_p2_NFIAn_R2","WT_D9_p2_NFIAn_R3",
                                   "WT_D9_p2_NFIAp_R1","WT_D9_p2_NFIAp_R2","WT_D9_p2_NFIAp_R3",
                                   "WT_D11_p2_NFIAp_R1","WT_D11_p2_NFIAp_R2","WT_D11_p2_NFIAp_R3",
                                   "WT_D5_pM_NFIAn_R1","WT_D5_pM_NFIAn_R2","WT_D5_pM_NFIAn_R3",
                                   "WT_D7_pM_NFIAn_R1","WT_D7_pM_NFIAn_R2","WT_D7_pM_NFIAn_R3",
                                   "WT_D9_pM_NFIAn_R1","WT_D9_pM_NFIAn_R2","WT_D9_pM_NFIAn_R3",
                                   "WT_D9_pM_NFIAp_R1","WT_D9_pM_NFIAp_R2","WT_D9_pM_NFIAp_R3",
                                   "WT_D11_pM_NFIAp_R1","WT_D11_pM_NFIAp_R2","WT_D11_pM_NFIAp_R3")
