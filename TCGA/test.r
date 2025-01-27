#** First back after Christmas. Checked access and viewing TCGA data **
data <- read.csv("/mainfs/ddnb/TCGA_data/full_TCGA-DLBC_matrix.csv")
View(data)

dim(data) # 60660 genes, 49 samples
