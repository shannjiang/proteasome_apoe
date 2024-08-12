library(Seurat)
#library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(reshape2)

#work_dir = '/mnt/mfs/ctcn/datasets/rosmap/rnaseq/dlpfcTissue/snrnaseq/annotation.2024-02-26/Seurat/'
#work_dir = '/mnt/mfs/ctcn/datasets/rosmap/rnaseq/dlpfcTissue/snrnaseq/annotation.2022-05-10/Seurat/'
work_dir = '/home/shann/Documents/synapse/ROSMAP/snrnaseq/DeJager/seurat_files/'
out_dir = '/home/shann/Documents/Natura/proteasome_apoe/DeJager_rosmap_snRNAseq/'
if(!dir.exists(out_dir)){dir.create(out_dir,recursive=T)}

#celltypes = c('cux2+','astrocytes','cux2-','inhibitory','microglia','oligodendrocytes','opcs')
celltypes = 'endo'
celltype = 'endo'

#rosmap_meta = read.csv(file = paste0(work_dir,'ROSMAP_clinical.csv'),header = T)

#proteasome_gene_annotation
proteasome_19s_genes = c(paste0('PSMC',c(1:6)), paste0('PSMD',c(1:14)))
proteasome_20s_genes = c(paste0('PSMA',c(1:7)),paste0('PSMB',c(1:7)))
immunoproteasome_genes = c(paste0('PSMB',c(8:10)))
proteasome_chaperon_genes = c(paste0('PSMG',c(1:4)))
proteasome_activator_genes = c(paste0('PSME',c(1:4)))
#proteasome_tf_genes = c('NFE2L1','MAFF','MAFG','MAFK')
#proteasome_tf_genes = c('NFE2L1','STAT1','IRF1','EGR1')
proteasome_tf_genes = c('NFE2L1','NFE2L2')
proteasome_receptor_genes = c('ADRM1')
housekeeper_genes = c('ACTB','GAPDH')

#substructure = rep(c('19S','20S','Immunoproteasome','Chaperon','Activator','TF','Receptor','housekeeper'),c(length(proteasome_19s_genes),length(proteasome_20s_genes),length(immunoproteasome_genes),length(proteasome_chaperon_genes),length(proteasome_activator_genes),length(proteasome_tf_genes),length(proteasome_receptor_genes),length(housekeeper_genes)))
substructure = rep(c('19S','20S','IP','Assembly chaperones','IP activator complex','TF','Receptor','housekeeper'),c(length(proteasome_19s_genes),length(proteasome_20s_genes),length(immunoproteasome_genes),length(proteasome_chaperon_genes),length(proteasome_activator_genes),length(proteasome_tf_genes),length(proteasome_receptor_genes),length(housekeeper_genes)))
gene = c(proteasome_19s_genes,proteasome_20s_genes,immunoproteasome_genes,proteasome_chaperon_genes,proteasome_activator_genes,proteasome_tf_genes,proteasome_receptor_genes,housekeeper_genes)
gene_substructure_df = data.frame(gene = gene,substructure = substructure)
rownames(gene_substructure_df) = gene
#gene_substructure_df = gene_substructure_df[!gene_substructure_df$substructure %in% c('Receptor','housekeeper'),]
gene_substructure_df = gene_substructure_df[!gene_substructure_df$substructure %in% c('Receptor','housekeeper'),]
gene_substructure_df$substructure[46] = gene_substructure_df$gene[46]
gene_substructure_df$substructure[47] = gene_substructure_df$gene[47]
#gene_substructure_df2 = gene_substructure_df[gene_substructure_df$substructure %in% c('20S','19S','IP'),]

##helper
cel_exp = function(x){
  x=log2(x+1)
  return(mean(x))
}

cel_exp_pct = function(x){
  return(sum(x != 0)/length(x))
}

cel_exp_ct = function(x){
  return(sum(x != 0))
}

for(celltype in celltypes){
#seu <- LoadH5Seurat(paste0(work_dir,celltype,'/',celltype,'.h5Seurat'))
#seuObj = load(paste0(work_dir,celltype,'_seu.rda'))
seuObj = load('/home/shann/Documents/synapse/ROSMAP/snrnaseq/DeJager/seurat_files/endo_seu.rda')
cts = seu@assays$RNA@counts
gene_substructure_df2 = gene_substructure_df[gene_substructure_df$substructure %in% c('20S','19S','IP'),]
cts = cts[rownames(cts) %in% gene_substructure_df2$gene,]
cts = as.data.frame(as.matrix(cts))
cts2 = as.data.frame(apply(cts,2,as.integer))
rownames(cts2) = rownames(cts)
cts = cts2
gene_substructure_df2 = gene_substructure_df2[rownames(gene_substructure_df2) %in% rownames(cts),]
gene_substructure_df2 = gene_substructure_df2[rownames(cts),]
cts$substructure = gene_substructure_df2$substructure
substr_cts = aggregate(.~substructure,cts,sum)
rownames(substr_cts) = substr_cts$substructure
substr_cts = substr_cts[,!colnames(substr_cts) %in% 'substructure']
substr_cts = as.data.frame(t(substr_cts))

meta = seu@meta.data
meta$diagnosis = with(meta,ifelse(ceradsc==1,'AD',ifelse(ceradsc==2,'AD','Control')))
substr_cts$subtype = with(meta,paste0(cell.type,'_',diagnosis))
cel_exp_df = aggregate(.~subtype,substr_cts,cel_exp)
cel_exp_long = melt(cel_exp_df)
colnames(cel_exp_long) = c('celltype','substructure','cell_exp')
cel_exp_pct_df = aggregate(.~subtype,substr_cts,cel_exp_pct)
cel_exp_pct_long = melt(cel_exp_pct_df)
colnames(cel_exp_pct_long) = c('celltype','substructure','cell_exp_pct')
cel_exp_ct_df = aggregate(.~subtype,substr_cts,cel_exp_ct)
cel_exp_ct_long = melt(cel_exp_ct_df)
colnames(cel_exp_ct_long) = c('celltype','substructure','cell_exp_ct')
dot_df = cel_exp_long
dot_df$cell_exp_pct = cel_exp_pct_long$cell_exp_pct
dot_df$cell_exp_ct = cel_exp_ct_long$cell_exp_ct
midpoint = median(dot_df$cell_exp)
#celltype order
celltypes = unique(meta$cell.type)
celltypes = celltypes[order(celltypes)]
celltypes = rep(celltypes,each=2)
celltypes = paste0(celltypes,'_',c('Control','AD'))
dot_df$celltype = factor(dot_df$celltype,levels = celltypes)
#graph
p = ggplot(dot_df, aes(x=celltype, y = substructure, color = cell_exp, size = cell_exp_pct)) +
  geom_point() + geom_text(label = dot_df$cell_exp_ct, nudge_x = 0, nudge_y = 0.25, check_overlap = T, color = 'black', size = 2) +
  scale_colour_gradient2(name = 'cell exp', low = 'blue', mid = 'white', high = 'red', midpoint = midpoint) +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
pdf(file = paste0(out_dir, 'DeJager_rosmap_snRNAseq_',celltype,'_subtype_by_diagnosis_dotplot_BWR.pdf'), width = 9, height = 5)
print(p)
dev.off()
}