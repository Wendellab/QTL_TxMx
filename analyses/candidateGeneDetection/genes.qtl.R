AA <- scan("genes.with.aa.variants", what="")
eff <- scan("genes.with.eff.variants", what="")
DGE <- scan("new.genes.expression", what="")
TF <- scan("TF.gene.names", what="")
other <- scan("genes.of.interest", what="")

Fiber_color_Ca               <- scan("Fiber_color_Ca.genes               ", what="")
Fiber_length_LnCV            <- scan("Fiber_length_LnCV.genes            ", what="")
Fiber_quality_NS             <- scan("Fiber_quality_NS.genes             ", what="")
Flower_SD                    <- scan("Flower_SD.genes                    ", what="")
Phenology_FBFF               <- scan("Phenology_FBFF.genes               ", what="")
Plant_architecture_PHFB1     <- scan("Plant_architecture_PHFB1.genes     ", what="")
Seed_FSW                     <- scan("Seed_FSW.genes                     ", what="")
Fiber_color.category         <- scan("Fiber_color.category.genes         ", what="")
Fiber_length_Ln              <- scan("Fiber_length_Ln.genes              ", what="")
Fiber_quality_SFCn           <- scan("Fiber_quality_SFCn.genes           ", what="")
Fruiting_habit.category      <- scan("Fruiting_habit.category.genes      ", what="")
Phenology_GB                 <- scan("Phenology_GB.genes                 ", what="")
Plant_architecture_PHFB2     <- scan("Plant_architecture_PHFB2.genes     ", what="")
Seed_SCW                     <- scan("Seed_SCW.genes                     ", what="")
Fiber_color_Cb               <- scan("Fiber_color_Cb.genes               ", what="")
Fiber_length_Lw              <- scan("Fiber_length_Lw.genes              ", what="")
Fiber_quality_TrS            <- scan("Fiber_quality_TrS.genes            ", what="")
Fruiting_habit_NF            <- scan("Fruiting_habit_NF.genes            ", what="")
Phenology_TNFF               <- scan("Phenology_TNFF.genes               ", what="")
Plant_architecture_PH        <- scan("Plant_architecture_PH.genes        ", what="")
Seed_SW                      <- scan("Seed_SW.genes                      ", what="")
Fiber_color_CL               <- scan("Fiber_color_CL.genes               ", what="")
Fiber_length_UQLw            <- scan("Fiber_length_UQLw.genes            ", what="")
Flower.category              <- scan("Flower.category.genes              ", what="")
Fruiting_habit_PHTN          <- scan("Fruiting_habit_PHTN.genes          ", what="")
Plant_architecture_BA        <- scan("Plant_architecture_BA.genes        ", what="")
Plant_architecture_SP        <- scan("Plant_architecture_SP.genes        ", what="")
Fiber_length.category        <- scan("Fiber_length.category.genes        ", what="")
Fiber_quality.category       <- scan("Fiber_quality.category.genes       ", what="")
Flower_CS                    <- scan("Flower_CS.genes                    ", what="")
Fruiting_habit_TNFB          <- scan("Fruiting_habit_TNFB.genes          ", what="")
Plant_architecture.category  <- scan("Plant_architecture.category.genes  ", what="")
Seed_AL                      <- scan("Seed_AL.genes                      ", what="")
Fiber_length_L25n            <- scan("Fiber_length_L25n.genes            ", what="")
Fiber_quality_Fine           <- scan("Fiber_quality_Fine.genes           ", what="")
Flower_PC                    <- scan("Flower_PC.genes                    ", what="")
Fruiting_habit_TN            <- scan("Fruiting_habit_TN.genes            ", what="")
Plant_architecture_FB1       <- scan("Plant_architecture_FB1.genes       ", what="")
Seed_BW                      <- scan("Seed_BW.genes                      ", what="")
Fiber_length_L5n             <- scan("Fiber_length_L5n.genes             ", what="")
Fiber_quality_MR             <- scan("Fiber_quality_MR.genes             ", what="")
Flower_PS                    <- scan("Flower_PS.genes                    ", what="")
Phenology.category           <- scan("Phenology.category.genes           ", what="")
Plant_architecture_FB2       <- scan("Plant_architecture_FB2.genes       ", what="")
Seed.category                <- scan("Seed.category.genes                ", what="")


master_list <- sort(unique(c(AA,eff,DGE,TF,other,Fiber_color.category,Fiber_color_Ca, Fiber_color_Cb, Fiber_color_CL, Fiber_length.category, 
Fiber_length_L25n, Fiber_length_L5n, Fiber_length_Ln, Fiber_length_LnCV, Fiber_length_Lw, Fiber_length_UQLw, Fiber_quality.category, Fiber_quality_Fine,Fiber_quality_MR, 
Fiber_quality_NS, Fiber_quality_SFCn, Fiber_quality_TrS, Flower.category, Flower_CS, Flower_PC, Flower_PS, Flower_SD, Fruiting_habit.category, Fruiting_habit_NF, 
Fruiting_habit_PHTN, Fruiting_habit_TN, Fruiting_habit_TNFB, Phenology.category, Phenology_FBFF, Phenology_GB, Phenology_TNFF, Plant_architecture.category, Plant_architecture_BA, 
Plant_architecture_FB1, Plant_architecture_FB2, Plant_architecture_PH, Plant_architecture_PHFB1, Plant_architecture_PHFB2, Plant_architecture_SP, Seed.category, Seed_AL, Seed_BW, 
Seed_FSW, Seed_SCW, Seed_SW)))

df <- as.data.frame(master_list)
row.names(df) <- master_list
df$AA <- df$master_list %in% AA
df$eff                         <- df$master_list %in% eff
df$DGE                         <- df$master_list %in% DGE
df$TF                          <- df$master_list %in% TF
df$other                       <- df$master_list %in% other
df$Fiber_color.category        <- df$master_list %in% Fiber_color.category
df$Fiber_color_Ca              <- df$master_list %in% Fiber_color_Ca 
df$Fiber_color_Cb              <- df$master_list %in% Fiber_color_Cb 
df$Fiber_color_CL              <- df$master_list %in% Fiber_color_CL 
df$Fiber_length.category       <- df$master_list %in% Fiber_length.category 
df$Fiber_length_L25n           <- df$master_list %in% Fiber_length_L25n 
df$Fiber_length_L5n            <- df$master_list %in% Fiber_length_L5n 
df$Fiber_length_Ln             <- df$master_list %in% Fiber_length_Ln 
df$Fiber_length_LnCV           <- df$master_list %in% Fiber_length_LnCV 
df$Fiber_length_Lw             <- df$master_list %in% Fiber_length_Lw 
df$Fiber_length_UQLw           <- df$master_list %in% Fiber_length_UQLw 
df$Fiber_quality.category      <- df$master_list %in% Fiber_quality.category 
df$Fiber_quality_Fine          <- df$master_list %in% Fiber_quality_Fine 
df$Fiber_quality_MR            <- df$master_list %in% Fiber_quality_MR 
df$Fiber_quality_NS            <- df$master_list %in% Fiber_quality_NS 
df$Fiber_quality_SFCn          <- df$master_list %in% Fiber_quality_SFCn 
df$Fiber_quality_TrS           <- df$master_list %in% Fiber_quality_TrS 
df$Flower.category             <- df$master_list %in% Flower.category 
df$Flower_CS                   <- df$master_list %in% Flower_CS 
df$Flower_PC                   <- df$master_list %in% Flower_PC 
df$Flower_PS                   <- df$master_list %in% Flower_PS 
df$Flower_SD                   <- df$master_list %in% Flower_SD 
df$Fruiting_habit.category     <- df$master_list %in% Fruiting_habit.category 
df$Fruiting_habit_NF           <- df$master_list %in% Fruiting_habit_NF 
df$Fruiting_habit_PHTN         <- df$master_list %in% Fruiting_habit_PHTN 
df$Fruiting_habit_TN           <- df$master_list %in% Fruiting_habit_TN 
df$Fruiting_habit_TNFB         <- df$master_list %in% Fruiting_habit_TNFB 
df$Phenology.category          <- df$master_list %in% Phenology.category 
df$Phenology_FBFF              <- df$master_list %in% Phenology_FBFF 
df$Phenology_GB                <- df$master_list %in% Phenology_GB 
df$Phenology_TNFF              <- df$master_list %in% Phenology_TNFF 
df$Plant_architecture.category <- df$master_list %in% Plant_architecture.category 
df$Plant_architecture_BA       <- df$master_list %in% Plant_architecture_BA 
df$Plant_architecture_FB1      <- df$master_list %in% Plant_architecture_FB1 
df$Plant_architecture_FB2      <- df$master_list %in% Plant_architecture_FB2 
df$Plant_architecture_PH       <- df$master_list %in% Plant_architecture_PH 
df$Plant_architecture_PHFB1    <- df$master_list %in% Plant_architecture_PHFB1 
df$Plant_architecture_PHFB2    <- df$master_list %in% Plant_architecture_PHFB2 
df$Plant_architecture_SP       <- df$master_list %in% Plant_architecture_SP 
df$Seed.category               <- df$master_list %in% Seed.category 
df$Seed_AL                     <- df$master_list %in% Seed_AL 
df$Seed_BW                     <- df$master_list %in% Seed_BW 
df$Seed_FSW                    <- df$master_list %in% Seed_FSW 
df$Seed_SCW                    <- df$master_list %in% Seed_SCW 
df$Seed_SW                     <- df$master_list %in% Seed_SW

df$master_list <- NULL

funmat <- read.table("annotated.Gohir.match", row.names=1, header=T, sep="\t")

bothmat <- merge(df,funmat,by="row.names",all.x=T)
bothmat$total <- rowSums(subset(bothmat,select=c(AA, eff, DGE, TF, other, relevance)))

write.table(bothmat, file="annotated.Gohir.tbl", sep="\t", quote=F, row.names=F)

#grep("category",names(bothmat), value=T) 

#[1] "Fiber_color.category"        "Fiber_length.category"       "Fiber_quality.category"      "Flower.category"             "Fruiting_habit.category"     "Phenology.category"          "Plant_architecture.category"
#[8] "Seed.category" 

FCmat <- subset(bothmat,Fiber_color.category==1,select=c(Row.names,AA, eff, DGE, TF, other, relevance,total,Best.hit.arabi.name,arabi.symbol,arabi.defline, description))
FCmat <- FCmat[order(-FCmat$total),]

FLmat <- subset(bothmat,Fiber_length.category==1,select=c(Row.names,AA, eff, DGE, TF, other, relevance,total,Best.hit.arabi.name,arabi.symbol,arabi.defline, description))
FLmat <- FLmat[order(-FLmat$total),]

FQmat <- subset(bothmat,Fiber_quality.category==1,select=c(Row.names,AA, eff, DGE, TF, other, relevance,total,Best.hit.arabi.name,arabi.symbol,arabi.defline, description))
FQmat <- FQmat[order(-FQmat$total),]

Flowermat <- subset(bothmat,Flower.category==1,select=c(Row.names,AA, eff, DGE, TF, other, relevance,total,Best.hit.arabi.name,arabi.symbol,arabi.defline, description))
Flowermat <- Flowermat[order(-Flowermat$total),]

Fruitmat <- subset(bothmat,Fruiting_habit.category==1,select=c(Row.names,AA, eff, DGE, TF, other, relevance,total,Best.hit.arabi.name,arabi.symbol,arabi.defline, description))
Fruitmat <- Fruitmat[order(-Fruitmat$total),]

Phenmat <- subset(bothmat,Phenology.category==1,select=c(Row.names,AA, eff, DGE, TF, other, relevance,total,Best.hit.arabi.name,arabi.symbol,arabi.defline, description))
Phenmat <- Phenmat[order(-Phenmat$total),]

Seedmat <- subset(bothmat,Seed.category==1,select=c(Row.names,AA, eff, DGE, TF, other, relevance,total,Best.hit.arabi.name,arabi.symbol,arabi.defline, description))
Seedmat <- Seedmat[order(-Seedmat$total),]

PAmat <- subset(bothmat,Plant_architecture.category==1,select=c(Row.names,AA, eff, DGE, TF, other, relevance,total,Best.hit.arabi.name,arabi.symbol,arabi.defline, description))
PAmat <- PAmat[order(-PAmat$total),]

write.table(FCmat, file="annotated.Gohir.FCmat.tbl", sep="\t", quote=F, row.names=F)
write.table(FLmat, file="annotated.Gohir.FLmat.tbl", sep="\t", quote=F, row.names=F)
write.table(FQmat, file="annotated.Gohir.FQmat.tbl", sep="\t", quote=F, row.names=F)
write.table(Flowermat, file="annotated.Gohir.Flowermat.tbl", sep="\t", quote=F, row.names=F)
write.table(Fruitmat, file="annotated.Gohir.Fruitmat.tbl", sep="\t", quote=F, row.names=F)
write.table(Phenmat, file="annotated.Gohir.Phenmat.tbl", sep="\t", quote=F, row.names=F)
write.table(Seedmat, file="annotated.Gohir.Seedmat.tbl", sep="\t", quote=F, row.names=F)
write.table(PAmat, file="annotated.Gohir.PAmat.tbl", sep="\t", quote=F, row.names=F)








