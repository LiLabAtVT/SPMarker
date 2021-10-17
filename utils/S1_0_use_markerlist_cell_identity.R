##this script is to make correlation between two matrix

args <- commandArgs(T)

ipt_gene_Cell_exp_mtx_fl   <- as.character(args[1]) 
ipt_cellmarker_mtx_fl  <- as.character(args[2]) 
ipt_output_dir <- as.character(args[3])

ipt_gene_Cell_exp_mtx_dt <- read.csv(ipt_gene_Cell_exp_mtx_fl,row.names = 1)
ipt_cellmarker_dt <- read.delim(ipt_cellmarker_mtx_fl,row.names = 1)

##put the gene cell exp ahead of cellmarker fl
corr_file <- cor(ipt_gene_Cell_exp_mtx_dt,ipt_cellmarker_dt)

##remove na
corr_file <- corr_file[complete.cases(corr_file),]

##select the max of the corr_file
max_order <- apply(corr_file,1,which.max)
max_cell_dt <- corr_file[names(max_order),max_order]

results <- data.frame(matrix(nrow=length(max_order),ncol=2))
rownames(results) <- names(max_order)
colnames(results) <- c('cell_type','prob')

for (i in 1:length(max_order)){
  col_order = max_order[i]
  cellname = names(max_order)[i]
  corr_val = corr_file[cellname,col_order]
  celltype = colnames(corr_file)[col_order]
  results[cellname,1] <- celltype
  results[cellname,2] <- corr_val
}

write.csv(results,paste0(ipt_output_dir,'/opt_meta.csv'))




