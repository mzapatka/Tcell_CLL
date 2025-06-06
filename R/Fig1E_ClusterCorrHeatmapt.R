
library("heatmaply")
library("plotly")
my_data <- read.delim("CyTOF2021_SampleClusterCells_CLL_LNs.txt")

my_data2 <- my_data[,-1]
rownames(my_data2) <- my_data[,1]	#name the rows by a column value

# Pearson correlation
res <- cor(my_data2)
round(res, 2)

# Spearman correlation (change res to res2!)
res2 <- cor(my_data2, method = "spearman")
round(res2, 2)


## We use this function to calculate a matrix of p-values from correlation tests
## https://stackoverflow.com/a/13112337/4747043
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}
p <- cor.test.p(my_data2)

#bootstraping
library(reshape2)
p_list <- c()
r_list <- c()
for (i in 1:nrow(my_data2)){	#go through each row of my_data2
  subset_data <- my_data2[-i,]	#remove one row at a time
  print("Size of subset matrix")
  print(dim(subset_data))
  print(rownames(subset_data))
  p_melt <- cbind(i,melt(cor.test.p(subset_data)))	#make a matrix of i-cell-cell-value pairs instead of a matrix
  r_melt <- cbind(i,melt(cor(subset_data)))
  p_list <- rbind(p_list, p_melt) 	#add the p-values result of this -1 sample "melted" list to the combined list
  r_list <- rbind(r_list, r_melt) 	#add the r result of this -1 sample "melted" list to a combined list
}
#now the list is complete, find which value is the worst-p and the r
colnames(p_list) = c("DropSample","Cell1","Cell2","value")
colnames(r_list) = c("DropSample","Cell1","Cell2","value")
cell_cell_combinations = unique(p_list[,c("Cell1","Cell2")])	#all possible combination of cell1-cell2

#this routine find the max-p and the corresponding r in all the values computed above
p_boostraped <- c()
r_boostraped <- c()
for (i in 1:nrow(cell_cell_combinations)){	#go through each row of my_data2
  p_combination <- p_list[which(p_list$Cell1 == cell_cell_combinations[i,1] & p_list$Cell2 == cell_cell_combinations[i,2]),]
  p_combination_max <- max(p_combination$value)
  max_drop_sample_index <- p_combination[which(p_combination$value == max(p_combination_max)),1][1]	#get the first drop-sample that have p-value = max(p-value)
  r_combination_max <- r_list[which(r_list$DropSample == max_drop_sample_index & r_list$Cell1 == cell_cell_combinations[i,1] & r_list$Cell2 == cell_cell_combinations[i,2]),]$value	#get the r corresponding = the p-value
  p_boostraped <- rbind(p_boostraped,data.frame(cell_cell_combinations[i,],value=p_combination_max,stringsAsFactors=F))
  r_boostraped <- rbind(r_boostraped,data.frame(cell_cell_combinations[i,],value=r_combination_max,stringsAsFactors=F))
}

#from melted format back to matrix
p_boostraped_df = acast(Cell1~Cell2,data=p_boostraped)
r_boostraped_df = acast(Cell1~Cell2,data=r_boostraped)

# Getting the dot size (p-val) right
testpval <- (-log10(p_boostraped_df))

test1 <- ifelse(testpval > 1.3,
                yes = 1.3,
                no = testpval)

# Plotting the figure
heatmaply(
  r_boostraped_df,
  node_type = "scatter",
  point_size_mat = (test1), 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation"),
  plot_method = "ggplot", limits = c(-1, 1), colors = cool_warm, grid_size = 0.1,
  dir.create("folder"))

pdf("~/Desktop/plot.pdf", width=10, height=8)
heatmaply(
  r_boostraped_df,
  node_type = "scatter",
  point_size_mat = (test1), 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation"),
  plot_method = "ggplot", limits = c(-1, 1), colors = cool_warm, grid_size = 1,
  return_ppxpy = TRUE,
  browseURL("heatmaply_plot.pdf"))
dev.off()



