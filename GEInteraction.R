library(tidyverse)
library(janitor)
library(lme4)

## original data
betaValue <- read_delim("00.MZ_355pair_vmeQTL/00.Original_info_file/450k_vCpGs_twins_beta_enmix_probeqced_990ids_ready.csv",
                        col_names = T, delim=",")
SNP <- read_delim("00.MZ_355pair_vmeQTL/02.Genotype_MZ_355pair/allvmeQTLs.Genotype",
                  col_names=T, delim="\t")
info <- read_delim("00.MZ_355pair_vmeQTL/00.Original_info_file/info_450k_blood_355pair_MZ_female.csv",
                   col_names=T, delim="\t") %>%
              select(Barcode, Sentrix_ID, KCLfam, Age, BMI_Nearest, Smoking_Nearest, CD8.naive, CD8T, Mono, Gran, MLR) %>% 
              mutate(Smoking_Grouped = ifelse(Smoking_Nearest %in% c(0,1),0,1)) %>%
              mutate(Smoking_Nearest = Smoking_Grouped) %>% select(-Smoking_Grouped)
#vmeQTLcis <- read_delim("03.TensorQTL.cis/cis.1027vCpG",
#                          col_names=T, delim="\t") %>% mutate(SNP=variant_id, CpG=phenotype_id) %>% select(SNP, CpG)
vmeQTL <- read_delim("03.runLDforall/vCpG.1017cis.656trans.removeMeanvar")

# roughly select the top SNPs for each CpG
out <- list()
for (i in c(1:nrow(vmeQTL))){
  print(i)
  temp <- vmeQTL[i,] 
  beta_temp <- betaValue %>% filter(cg==temp$CpG) %>% t() %>% data.frame() %>% row_to_names(row_number = 1)
  beta_temp$Barcode <- rownames(beta_temp)
  colnames(beta_temp) <- c("cg","Barcode")
  beta_temp1 <- merge(beta_temp, info, by.x="Barcode") %>% mutate(cg=as.numeric(cg), KCLfam=as.numeric(KCLfam))
  
  SNP_temp <- SNP %>% filter(info==temp$SNP) %>% select(-CHR,-SNP,-`(C)M`,-POS,-COUNTED,-ALT) %>%
    t() %>% data.frame() %>% row_to_names(row_number = 1)
  SNP_temp <- SNP_temp %>% mutate(KCLfam=rownames(SNP_temp))
  
  df <- merge(SNP_temp, beta_temp1, by.x="KCLfam") %>% select(-Barcode) %>% mutate_if(is.character, as.numeric)
  colnames(df) <- c("KCLfam","SNP",colnames(df)[3:12])
  
  df$resi.BMI <- residuals(lmer(cg ~ (1|KCLfam) + (1|Sentrix_ID) + Age + Smoking_Nearest + CD8.naive + CD8T + Mono + Gran, data=df))
  df$resi.Smoke <- residuals(lmer(cg ~ (1|KCLfam) + (1|Sentrix_ID) + Age + BMI_Nearest + CD8.naive + CD8T + Mono + Gran, data=df))
  df$resi.CD8naive <- residuals(lmer(cg ~ (1|KCLfam) + (1|Sentrix_ID) + Age + BMI_Nearest + Smoking_Nearest + CD8T + Mono + Gran, data=df))
  df$resi.CD8T <- residuals(lmer(cg ~ (1|KCLfam) + (1|Sentrix_ID) + Age + BMI_Nearest + Smoking_Nearest + CD8.naive + Mono + Gran, data=df))
  df$resi.Mono <- residuals(lmer(cg ~ (1|KCLfam) + (1|Sentrix_ID) + Age + BMI_Nearest + Smoking_Nearest + CD8.naive + CD8T + Gran, data=df))
  df$resi.Gran <- residuals(lmer(cg ~ (1|KCLfam) + (1|Sentrix_ID) + Age + BMI_Nearest + Smoking_Nearest + CD8.naive + CD8T + Mono, data=df))
  df$resi.MLR <- residuals(lmer(cg ~ (1|KCLfam) + (1|Sentrix_ID) + Age + BMI_Nearest + Smoking_Nearest + CD8.naive + CD8T + Mono + Gran, data=df))
  
  EffectSize <- data.frame(vCpG = temp$CpG,
                           vmeQTL = temp$SNP,
                           Measure = "effectsize",
                           BMI = summary(lm(resi.BMI ~ SNP*BMI_Nearest, data=df))$coef[[4]],
                           Smoking = summary(lm(resi.Smoke ~ SNP*Smoking_Nearest, data=df))$coef[[4]],
                           CD8.naive = summary(lm(resi.CD8naive ~ SNP*CD8.naive, data=df))$coef[[4]],
                           CD8T = summary(lm(resi.CD8T ~ SNP*CD8T, data=df))$coef[[4]],
                           Mono = summary(lm(resi.Mono ~ SNP*Mono, data=df))$coef[[4]],
                           Gran = summary(lm(resi.Gran ~ SNP*Gran, data=df))$coef[[4]],
                           MLR = summary(lm(resi.MLR ~ SNP*MLR, data=df))$coef[[4]])
  Pvalue <- data.frame(vCpG = temp$CpG,
                       vmeQTL = temp$SNP,
                       Measure = "pvalue",
                       BMI = summary(lm(resi.BMI ~ SNP*BMI_Nearest, data=df))$coef[[16]],
                       Smoking = summary(lm(resi.Smoke ~ SNP*Smoking_Nearest, data=df))$coef[[16]],
                       CD8.naive = summary(lm(resi.CD8naive ~ SNP*CD8.naive, data=df))$coef[[16]],
                       CD8T = summary(lm(resi.CD8T ~ SNP*CD8T, data=df))$coef[[16]],
                       Mono = summary(lm(resi.Mono ~ SNP*Mono, data=df))$coef[[16]],
                       Gran = summary(lm(resi.Gran ~ SNP*Gran, data=df))$coef[[16]],
                       MLR = summary(lm(resi.MLR ~ SNP*MLR, data=df))$coef[[16]])
  out[[i]] <- rbind(EffectSize, Pvalue)
}
output <- bind_rows(out)
write.table(output, "05.GEI/GEI.vmeQTL.230512.1698assoc", quote=F, row.names=F, col.names = T, sep="\t")

