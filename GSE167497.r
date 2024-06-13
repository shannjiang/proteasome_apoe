library(rhdf5)
library(Matrix)
library(Seurat)

work_dir = '/home/shann/Documents/GEO/GSE167497/samples/'

samples = list.files(work_dir)
samples = gsub('_raw_gene_bc_matrices_h5\\.h5','',samples)
samples2 = samples[1:3]
sample = samples2[1]

h5_to_ctMtx = function(sample = sample){
  gsmid = sub('_.*','',sample)
  h5 = h5read(paste0(work_dir,sample,'_raw_gene_bc_matrices_h5.h5'),'mm10-1.2.0_premrna')
  barcodes = paste0(gsmid,'_',h5$barcodes)
  counts <- sparseMatrix(
    dims = h5$shape,
    i = as.numeric(h5$indices),
    p = as.numeric(h5$indptr),
    x = as.numeric(h5$data),
    index1 = FALSE
  )
  #denseCts = as.matrix(counts)
  colnames(counts) <- barcodes
  rownames(counts) <- h5$gene_names
  denseCts = t(sapply(by(counts,rownames(counts),colSums),identity))
  counts = Matrix(denseCts,sparse = TRUE)
  remove(denseCts)
  return(counts)
}

cts_lst = list()
for(sample in samples2){
  #h5 = H5Fopen(paste0(work_dir,sample,'_raw_gene_bc_matrices_h5.h5'))
  print(paste0('process ',sample))
  cts_lst[length(cts_lst)+1] = h5_to_ctMtx(sample)
}

cts = do.call(cbind,cts_lst)
seuObj = CreateSeuratObject(cts, project = "GSE167497", assay = "RNA")
