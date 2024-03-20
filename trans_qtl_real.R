#### matrix eQTL
library("MatrixEQTL")
args <- commandArgs(trailingOnly = TRUE)
perm <- as.numeric(args[1])*100 + as.numeric(args[2]) - 1
base.dir = "/01.MZ_vmeQTL/00.MZ_355pair_vmeQTL/"
useModel = modelLINEAR
SNP_file_name = paste(base.dir, "02.Genotype_MZ_355pair/allSNPs.additive.MHF5", sep="")
CpG_file_name = paste("03.MatrixQTL.trans/meQTL/input/450k_blood_beta_noShuffledinput_mean_", perm, sep="")
snps_location_file_name = paste(base.dir, "02.Genotype_MZ_355pair/SNPlocs.matrixEQTL", sep="")
cpgs_location_file_name = paste(base.dir, "01.DNAm_MZ_355pair/CpGlocs.matrixEQTL", sep="");
output_file_name = paste("03.MatrixQTL.perm.trans/meQTL/output/Trans_out_real", perm, sep="")
pvOutputThreshold = 5e-2
pvOutputThreshold.cis = 1e-1

snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipColumns = 7
snps$fileSliceSize = 2000
snps$LoadFile(SNP_file_name)

probes = SlicedData$new()
probes$fileDelimiter = ","
probes$fileOmitCharacters = "NA"
probes$fileSkipColumns = 1
probes$fileSliceSize = 2000
probes$LoadFile(CpG_file_name)

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
cpgspos = read.table(cpgs_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = probes,
  cvrt = SlicedData$new(),
  output_file_name = output_file_name,
  # output_file_name.cis = "",
  pvOutputThreshold = pvOutputThreshold,
  pvOutputThreshold.cis = pvOutputThreshold.cis,
  useModel = useModel,
  snpspos = snpspos,
  genepos = cpgspos,
  cisDist = 1e6,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE)

saveRDS(me, paste("03.MatrixQTL.perm.trans/meQTL/output/transmeQTL.real", perm, ".rds", sep=""))
