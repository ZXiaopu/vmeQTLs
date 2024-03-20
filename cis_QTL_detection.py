import pandas as pd
import torch
import tensorqtl
import sys
import pickle
from tensorqtl import genotypeio, cis, trans

chr = str(sys.argv[1])
plink_prefix_path = './TensorQTLinput/chr'+chr+'.additive.MHF5'
expression_bed= './TensorQTLinput/allprobe_noNA.variance.bed.gz'
covariates_file = './TensorQTLinput/covariates.variance'
prefix='cis.vmeQTL.chr'+chr+'.association.1103'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df['chr']=='chr'+chr], 
                phenotype_pos_df.loc[phenotype_pos_df['chr']=='chr'+chr], 
                prefix,window=1000000)

cis_df = cis.map_cis(genotype_df, variant_df, 
                     phenotype_df.loc[phenotype_pos_df['chr']=='chr'+chr],
                     phenotype_pos_df.loc[phenotype_pos_df['chr']=='chr'+chr],
                     window=1000000, nperm=10000)
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)
cis_df = pd.read_csv('cis.vmeQTL.chr'+chr+'.perm10000.csv', index_col=1)

indep_df = cis.map_independent(genotype_df, variant_df, cis_df,
                               phenotype_df.loc[phenotype_pos_df['chr']=='chr'+chr], phenotype_pos_df.loc[phenotype_pos_df['chr']=='chr'+chr], covariates_df)
indep_df.to_csv('cis.vmeQTL.indep.chr'+chr+'.csv')
