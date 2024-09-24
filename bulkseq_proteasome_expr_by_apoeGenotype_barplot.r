library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)

msbb_dir = '/home/shann/Documents/Natura/proteosome/MSBB/transcriptomics/'
mayo_dir = '/home/shann/Documents/Natura/proteosome/MAYO/transcriptomics/'
rosmap_dir = '/home/shann/Documents/Natura/proteosome/ROSMAP/transcriptomics/'
out_dir = '/home/shann/Documents/Natura/proteasome_apoe/bulkseq/'
proteasome_gene_id_df = read.table(file = paste0('/home/shann/Documents/Natura/proteosome/gene_id.txt'), header = T, sep = '\t')
proteasome_gene_id_df = proteasome_gene_id_df[!proteasome_gene_id_df$OGS %in% c('ADRM1','MAFF','MAFG','MAFK'),]

#helper
##rename gene to OGS
gene_rename = function(df, gene_id_df){
  new_id_df = gene_id_df[gene_id_df$ENSG_ID %in% rownames(df),]
  new_df = df[rownames(df) %in% gene_id_df$ENSG_ID,]
  rownames(new_id_df) = new_id_df$ENSG_ID
  new_id_df = new_id_df[rownames(new_df),]
  rownames(new_df) = new_id_df$OGS
  return(new_df)
}

br_meta_gen = function(vst = msbb_bm10_vst, meta = msbb_meta){
  br_meta = meta[rownames(meta) %in% colnames(vst),]
  br_meta = br_meta[colnames(vst),]
  return(br_meta)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

sem_cal = function(vec){
  return(sd(vec)/sqrt(length(vec)))
}

#proteasome_gene_annotation
proteasome_19s_genes = c(paste0('PSMC',c(1:6)), paste0('PSMD',c(1:14)))
proteasome_20s_genes = c(paste0('PSMA',c(1:7)),paste0('PSMB',c(1:7)))
immunoproteasome_genes = c(paste0('PSMB',c(8:10)))
proteasome_chaperon_genes = c(paste0('PSMG',c(1:4)))
proteasome_activator_genes = c(paste0('PSME',c(1:4)))
#proteasome_tf_genes = c('NFE2L1','MAFF','MAFG','MAFK')
#proteasome_tf_genes = c('NFE2L1','IRF1','STAT1','EGR1')
proteasome_tf_genes = c('NFE2L1','NFE2L2')
proteasome_receptor_genes = c('ADRM1')
housekeeper_genes = c('ACTB','GAPDH')

substructure = rep(c('19S','20S','Immunoproteasome','Chaperon','Activator','TF','Receptor','housekeeper'),c(length(proteasome_19s_genes),length(proteasome_20s_genes),length(immunoproteasome_genes),length(proteasome_chaperon_genes),length(proteasome_activator_genes),length(proteasome_tf_genes),length(proteasome_receptor_genes),length(housekeeper_genes)))
gene = c(proteasome_19s_genes,proteasome_20s_genes,immunoproteasome_genes,proteasome_chaperon_genes,proteasome_activator_genes,proteasome_tf_genes,proteasome_receptor_genes,housekeeper_genes)
gene_substructure_df = data.frame(gene = gene,substructure = substructure)
rownames(gene_substructure_df) = gene
#gene_substructure_df = gene_substructure_df[!gene_substructure_df$substructure %in% c('Receptor','housekeeper'),]
gene_substructure_df = gene_substructure_df[!gene_substructure_df$substructure %in% c('Receptor','housekeeper'),]
gene_substructure_df$broadSubstr = with(gene_substructure_df,ifelse(substructure == '19S','CP',ifelse(substructure == '20S','CP',ifelse(substructure == 'Immunoproteasome','IP',ifelse(substructure == 'Chaperon','Chaperon',ifelse(substructure == 'Activator','IP','TF'))))))


##vst
#MSBB vst
msbb_bm10_vst = read.csv(file = paste0(msbb_dir,'MSBB_BM_10_expr_vst.csv'), header = T, row.names = 1)
msbb_bm10_prot_vst = msbb_bm10_vst[rownames(msbb_bm10_vst) %in% proteasome_gene_id_df$ENSG_ID,]
msbb_bm22_vst = read.csv(file = paste0(msbb_dir,'MSBB_BM_22_expr_vst.csv'), header = T, row.names = 1)
msbb_bm22_prot_vst = msbb_bm22_vst[rownames(msbb_bm22_vst) %in% proteasome_gene_id_df$ENSG_ID,]
msbb_bm36_vst = read.csv(file = paste0(msbb_dir,'MSBB_BM_36_expr_vst.csv'), header = T, row.names = 1)
msbb_bm36_prot_vst = msbb_bm36_vst[rownames(msbb_bm36_vst) %in% proteasome_gene_id_df$ENSG_ID,]
msbb_bm44_vst = read.csv(file = paste0(msbb_dir,'MSBB_BM_44_expr_vst.csv'), header = T, row.names = 1)
msbb_bm44_prot_vst = msbb_bm44_vst[rownames(msbb_bm44_vst) %in% proteasome_gene_id_df$ENSG_ID,]
#MAYO vst
mayo_cer_vst = read.csv(file = paste0(mayo_dir,'MAYO_CER_expr_vst.csv'), header = T, row.names = 1)
mayo_tcx_vst = read.csv(file = paste0(mayo_dir,'MAYO_TCX_expr_vst.csv'), header = T, row.names = 1)
colnames(mayo_cer_vst) = gsub('^X','',colnames(mayo_cer_vst))
colnames(mayo_tcx_vst) = gsub('^X','',colnames(mayo_tcx_vst))
mayo_cer_prot_vst = mayo_cer_vst[rownames(mayo_cer_vst) %in% proteasome_gene_id_df$ENSG_ID,]
mayo_tcx_prot_vst = mayo_tcx_vst[rownames(mayo_tcx_vst) %in% proteasome_gene_id_df$ENSG_ID,]
#ROSMAP vst
rosmap_vst = read.csv(file = paste0(rosmap_dir,'ROSMAP_expr_vst.csv'), header = T, row.names = 1)
colnames(rosmap_vst) = gsub('^X','',colnames(rosmap_vst))
rosmap_prot_vst = rosmap_vst[rownames(rosmap_vst) %in% proteasome_gene_id_df$ENSG_ID,]

##rename gene to official gene symbol
#msbb
msbb_bm10_prot_vst = gene_rename(msbb_bm10_prot_vst, proteasome_gene_id_df)
msbb_bm22_prot_vst = gene_rename(msbb_bm22_prot_vst, proteasome_gene_id_df)
msbb_bm36_prot_vst = gene_rename(msbb_bm36_prot_vst, proteasome_gene_id_df)
msbb_bm44_prot_vst = gene_rename(msbb_bm44_prot_vst, proteasome_gene_id_df)
#mayo
mayo_cer_prot_vst = gene_rename(mayo_cer_prot_vst, proteasome_gene_id_df)
mayo_tcx_prot_vst = gene_rename(mayo_tcx_prot_vst, proteasome_gene_id_df)
#rosmap
rosmap_prot_vst = gene_rename(rosmap_prot_vst, proteasome_gene_id_df)

##meta
msbb_meta = read.csv(file = paste0(msbb_dir,'MSBB_meta.csv'), header = T, row.names = 1)
colnames(msbb_meta) = gsub('Braak','braaksc',colnames(msbb_meta))
mayo_meta = read.csv(file = paste0(mayo_dir,'MAYO_meta.csv'), header = T, row.names = 1)
rosmap_meta = read.csv(file = paste0(rosmap_dir, 'ROSMAP_meta.csv'), header = T, row.names = 1)

msbb_meta$braak = paste0('braak',msbb_meta$braaksc)
mayo_meta$braak = paste0('braak',mayo_meta$Braak)
mayo_meta$braak = gsub('braakNA',NA,mayo_meta$braak)
rosmap_meta$braak = paste0('braak',rosmap_meta$braaksc)
msbb_meta$sampleID = rownames(msbb_meta)
mayo_meta$sampleID = rownames(mayo_meta)
rosmap_meta$sampleID = rownames(rosmap_meta)
msbb_meta$dataset = 'msbb'
mayo_meta$dataset = 'mayo'
rosmap_meta$dataset = 'rosmap'

#apoe4status
msbb_meta$apoe4status = with(msbb_meta,ifelse(apoeGenotype==22,'APOE4noncarrier',ifelse(apoeGenotype==23,'APOE4noncarrier',ifelse(apoeGenotype==24,'APOE4carrier',ifelse(apoeGenotype==33,'APOE4noncarrier',ifelse(apoeGenotype==34,'APOE4carrier',ifelse(apoeGenotype==44,'APOE4carrier','NA')))))))
rosmap_meta$apoe4status = with(rosmap_meta,ifelse(apoe_genotype==22,'APOE4noncarrier',ifelse(apoe_genotype==23,'APOE4noncarrier',ifelse(apoe_genotype==24,'APOE4carrier',ifelse(apoe_genotype==33,'APOE4noncarrier',ifelse(apoe_genotype==34,'APOE4carrier',ifelse(apoe_genotype==44,'APOE4carrier','NA')))))))
mayo_meta$apoe4status = with(mayo_meta,ifelse(apoeGenotype==22,'APOE4noncarrier',ifelse(apoeGenotype==23,'APOE4noncarrier',ifelse(apoeGenotype==24,'APOE4carrier',ifelse(apoeGenotype==33,'APOE4noncarrier',ifelse(apoeGenotype==34,'APOE4carrier',ifelse(apoeGenotype==44,'APOE4carrier','NA')))))))
#apoeGenotype
rosmap_meta$apoeGenotype = rosmap_meta$apoe_genotype
#msbb_meta$apoeGenotype = factor(msbb_meta$apoeGenotype,levels = c(22,23,33,24,34,44))
#rosmap_meta$apoeGenotype = factor(rosmap_meta$apoeGenotype,levels = c(22,23,33,24,34,44))
#mayo_meta$apoeGenotype = factor(mayo_meta$apoeGenotype,levels = c(22,23,33,24,34,44))
#remove AD genotype 24 in msbb
msbb_meta = msbb_meta %>% filter(!(diagnosis %in% 'AD' & apoeGenotype %in% c(24)))

msbb_bm10_meta = br_meta_gen(msbb_bm10_vst,msbb_meta)
msbb_bm22_meta = br_meta_gen(msbb_bm22_vst,msbb_meta)
msbb_bm36_meta = br_meta_gen(msbb_bm36_vst,msbb_meta)
msbb_bm44_meta = br_meta_gen(msbb_bm44_vst,msbb_meta)
mayo_tcx_meta = br_meta_gen(mayo_tcx_vst,mayo_meta)
mayo_cer_meta = br_meta_gen(mayo_cer_vst,mayo_meta)
rosmap_meta = rosmap_meta[colnames(rosmap_prot_vst),]

msbb_bm10_meta$diagnosis = factor(msbb_bm10_meta$diagnosis,levels = c('Control','AD'))
msbb_bm22_meta$diagnosis = factor(msbb_bm22_meta$diagnosis,levels = c('Control','AD'))
msbb_bm36_meta$diagnosis = factor(msbb_bm36_meta$diagnosis,levels = c('Control','AD'))
msbb_bm44_meta$diagnosis = factor(msbb_bm44_meta$diagnosis,levels = c('Control','AD'))
mayo_tcx_meta$diagnosis = factor(mayo_tcx_meta$diagnosis,levels = c('Control','AD'))
mayo_cer_meta$diagnosis = factor(mayo_cer_meta$diagnosis,levels = c('Control','AD'))
rosmap_meta$diagnosis = factor(rosmap_meta$diagnosis,levels = c('Control','AD'))

#boxplot helper
barplot_fig_gen = function(vst = msbb_bm10_prot_vst, meta = msbb_bm10_meta, group = c('Control','AD'), gene_substructure_df = gene_substructure_df, substructures = c('CP','IP')){
  commonSamples = intersect(colnames(vst),rownames(meta))
  vst = vst[,colnames(vst) %in% commonSamples]
  meta = meta[rownames(meta) %in% commonSamples,]
  meta = meta[colnames(vst),]
  gene_substructure_df2 = gene_substructure_df[gene_substructure_df$broadSubstr %in% substructures,]
  commonGenes = intersect(rownames(vst),rownames(gene_substructure_df2))
  vst = vst[rownames(vst) %in% commonGenes,]
  gene_substructure_df2 = gene_substructure_df2[rownames(gene_substructure_df2) %in% commonGenes,]
  vst = vst[rownames(gene_substructure_df2),]
  vst$broadSubstr = gene_substructure_df2$broadSubstr
  vst_agg = aggregate(.~broadSubstr,vst,mean)
  rownames(vst_agg) = vst_agg$broadSubstr
  vst_agg = vst_agg[,-1]
  vst_agg = as.data.frame(t(vst_agg))
  vst_agg$diagnosis = meta$diagnosis
  vst_agg$apoeGenotype = meta$apoeGenotype
  #vst_agg <- vst_agg %>% filter(!(diagnosis %in% 'AD' & apoeGenotype %in% c(24)))
  vst_agg = vst_agg[!is.na(vst_agg$apoeGenotype),]
  vst_agg$apoeGenotype = factor(vst_agg$apoeGenotype,levels = sort(unique(vst_agg$apoeGenotype),decreasing = F))
  vst_agg$diagnosis = factor(vst_agg$diagnosis,levels = group)

  vst_agg_mean = aggregate(.~diagnosis+apoeGenotype,vst_agg,mean)
  vst_agg_sem = aggregate(.~diagnosis+apoeGenotype,vst_agg,sem_cal)
  vst_aggl_mean = melt(vst_agg_mean,id.vars = c('diagnosis','apoeGenotype'),value.name = 'mean')
  vst_aggl_sem = melt(vst_agg_sem,id.vars = c('diagnosis','apoeGenotype'),value.name = 'sem')
  vst_aggl_mean$sem = vst_aggl_sem$sem
  #vst_agg2 = aggregate(.~diagnosis+apoeGenotype,vst_agg,mean)
  colnames(vst_aggl_mean)[3] = 'proteasome_type'
  p = ggplot(vst_aggl_mean, aes(x=apoeGenotype, y=mean, fill=proteasome_type, color = diagnosis)) + geom_bar(stat = 'identity', position = 'dodge') + geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .2, position = position_dodge(.9)) + scale_fill_manual(values = c('grey90','grey50'))
  return(p)
}

#proteasome type and diagnosis together
#msbb bm10
p = barplot_fig_gen(msbb_bm10_prot_vst,msbb_bm10_meta,c('Control','AD'),gene_substructure_df,c('CP','IP'))
pdf(paste0(out_dir,'msbb_bm10_proteasome_by_apoeGenotype_barplot.pdf'),width = 7, height = 7)
print(p)
dev.off()
#msbb bm22
p = barplot_fig_gen(msbb_bm22_prot_vst,msbb_bm22_meta,c('Control','AD'),gene_substructure_df,c('CP','IP'))
pdf(paste0(out_dir,'msbb_bm22_proteasome_by_apoeGenotype_barplot.pdf'),width = 7, height = 7)
print(p)
dev.off()
#msbb bm36
p = barplot_fig_gen(msbb_bm36_prot_vst,msbb_bm36_meta,c('Control','AD'),gene_substructure_df,c('CP','IP'))
pdf(paste0(out_dir,'msbb_bm36_proteasome_by_apoeGenotype_barplot.pdf'),width = 7, height = 7)
print(p)
dev.off()
#msbb bm44
p = barplot_fig_gen(msbb_bm44_prot_vst,msbb_bm44_meta,c('Control','AD'),gene_substructure_df,c('CP','IP'))
pdf(paste0(out_dir,'msbb_bm44_proteasome_by_apoeGenotype_barplot.pdf'),width = 7, height = 7)
print(p)
dev.off()

#mayo tcx
mayo_tcx_meta = br_meta_gen(mayo_tcx_vst,mayo_meta)
mayo_tcx_meta$diagnosis2 = with(mayo_tcx_meta,ifelse(diagnosis == 'control','Control',ifelse(diagnosis == 'Alzheimer Disease','AD',ifelse(diagnosis == 'pathological aging','PA','PSP'))))
mayo_tcx_meta$diagnosis = mayo_tcx_meta$diagnosis2
mayo_tcx_meta = mayo_tcx_meta[mayo_tcx_meta$diagnosis %in% c('Control','AD'),]
mayo_tcx_prot_vst = mayo_tcx_prot_vst[,colnames(mayo_tcx_prot_vst) %in% rownames(mayo_tcx_meta)]
p = barplot_fig_gen(mayo_tcx_prot_vst,mayo_tcx_meta,c('Control','AD'),gene_substructure_df,c('CP','IP'))
pdf(paste0(out_dir,'mayo_tcx_proteasome_by_apoeGenotype_barplot.pdf'),width = 7, height = 7)
print(p)
dev.off()
#mayo cer
mayo_cer_meta = br_meta_gen(mayo_cer_vst,mayo_meta)
mayo_cer_meta$diagnosis2 = with(mayo_cer_meta,ifelse(diagnosis == 'control','Control',ifelse(diagnosis == 'Alzheimer Disease','AD',ifelse(diagnosis == 'pathological aging','PA','PSP'))))
mayo_cer_meta$diagnosis = mayo_cer_meta$diagnosis2
mayo_cer_meta = mayo_cer_meta[mayo_cer_meta$diagnosis %in% c('Control','AD'),]
mayo_cer_prot_vst = mayo_cer_prot_vst[,colnames(mayo_cer_prot_vst) %in% rownames(mayo_cer_meta)]
p = barplot_fig_gen(mayo_cer_prot_vst,mayo_cer_meta,c('Control','AD'),gene_substructure_df,c('CP','IP'))
pdf(paste0(out_dir,'mayo_cer_proteasome_by_apoeGenotype_barplot.pdf'),width = 7, height = 7)
print(p)
dev.off()
#rosmap dlpfc
p = barplot_fig_gen(rosmap_prot_vst,rosmap_meta,c('Control','AD'),gene_substructure_df,c('CP','IP'))
pdf(paste0(out_dir,'rosmap_dlpfc_proteasome_by_apoeGenotype_barplot.pdf'),width = 7, height = 7)
print(p)
dev.off()
