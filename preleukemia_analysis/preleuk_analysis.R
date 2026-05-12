#Preleukemia scRNAseq Sample Analysis (Part 1 in R)
#Naraen Palanikumar
#This is part 1 of the analysis and dashboard completed of publicly available scRNAseq/gene signature/reference atlas data from this paper (https://doi.org/10.1016/j.xgen.2023.100426) and AML clinical data from NCI.

# --- Data Summary ---
#Raw: Contains scRNAseq data from 38 samples as well as a metadata file
#Reference: Contains reference file for cell annotation
#Seurat_Inter: Contains pre-normalized data to integrate
#R_Final: Contains the final integrated seurat file for further analysis

# --- Libraries ---
library(Seurat)     #core scRNAseq data manipulation package
library(Matrix)     #handling sparse count matrices
library(hdf5r)      #reading 10X seq files with .h5 extension
library(patchwork)  #combining plots
library(ggplot2)    #further plotting
library(SeuratDisk) #converting file formats

# --- Load Data ---
# Use a project-relative working directory so collaborators can run this script
# without changing hardcoded user-specific paths.
resolve_project_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path)))
  }
  return(getwd())
}

setwd(resolve_project_dir())
samples_df <- read.csv('data/raw/sample_metadata.csv') #reads sample metadata

for (i in 1:nrow(samples_df)) { #loops through each of the samples
  filepath <- samples_df$file_path[i] #path of sample data itself
  sample_id <- samples_df$sample_id[i] #sample id name
  condition <- samples_df$condition[i] #condition (i.e. WT vs Mutant)
  
  print(paste("Stage 1 - Processing sample: ", i, "/", nrow(samples_df), " (", sample_id, ")")) #keeps track of progress
  
  counts <- Read10X_h5(filepath) #reads the h5 file
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200) #creates seurat obj
  seurat_obj$sample_id <- sample_id #appends sample_id
  seurat_obj$condition <- condition #appends condition
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-") #adds mtdna column
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10) #QC cells based on filter
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE) #normalizes the data
  
  saveRDS(seurat_obj, file = paste0("data/seurat_inter/stage1_", sample_id, ".rds")) #saves to file to seurat_inter
  
  rm(seurat_obj, counts) #removes unwanted obj from memory
  gc() #garbage collection
}

# --- Read Intermediate Data ---
rds_files <- list.files("data/seurat_inter", pattern = "*.rds", full.names = TRUE) #lists the intermediate files
seurat_list <- lapply(rds_files, readRDS) #reads all the rds files
seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, reference = 1) #uses the first file as a reference anchor
seurat_integrated <- IntegrateData(anchorset = seurat_anchors) #integrates the data using anchors

rm(rds_files, seurat_list, seurat_anchors) #removes unnecessary files
gc() #garbage collection

# --- Finding Transfer Anchors ---
reference_atlas <- readRDS("data/reference/reference_atlas.rds") #reads the reference atlas
reference_atlas <- FindVariableFeatures(reference_atlas, selection.method = "vst", nfeatures = 2000) #finds variable features
reference_atlas <- NormalizeData(reference_atlas, verbose = FALSE) #normalizes the data
reference_atlas <- ScaleData(reference_atlas, verbose = FALSE) #scales the data
reference_atlas <- RunPCA(reference_atlas, npcs = 30, verbose = FALSE) #runs a principal component analysis

transfer_anchors <- FindTransferAnchors( #creates reference anchors
  reference = reference_atlas,
  query = seurat_integrated,
  dims = 1:30,
  reference.reduction = "pca",
  features = VariableFeatures(reference_atlas)
)

gc() #garbage collection

# --- Process Integrated Data ---
DefaultAssay(seurat_integrated) <- "integrated" #selects that we want to run a integrated workflow
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE) #scales the data
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30, verbose = FALSE) #runs a principal component analysis

predictions <- TransferData( #transfers cell type labels from ref to dataset
  anchorset = transfer_anchors,
  refdata = reference_atlas$celltype,
  weight.reduction = seurat_integrated[["pca"]],
  dims = 1:30
)

seurat_integrated <- AddMetaData(seurat_integrated, metadata = predictions) #adds this metadata to seurat integrated

seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:30) #runs a umap
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", dims = 1:30) #finds nearest neighbors
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5) #finds similar clusters

rm(predictions, transfer_anchors, reference_atlas) #removes unnecessary files
gc() #garbage collection

# --- Plotting Results ---
p1 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE) #predicted UMAP
p2 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "condition", cols = c("Mutant" = "purple", "WT" = "yellow")) #WT vs Mutant
png("preleuk_dashboard/plots/UMAP_Integrated.png", width = 12, height = 6, units = "in", res = 300) #gives a path for the plot
print(p1 + p2) #saves the plot
dev.off() #turns off graphical device

# --- Cell Abundance Statistical Analysis ---
cell_counts <- table(seurat_integrated$predicted.id, seurat_integrated$condition) #puts the cell counts from each group into a table
print(cell_counts) #prints this table
chisq.test(cell_counts) #prints the chi-sq test for difference between these 

cell_counts_df <- as.data.frame(cell_counts) #converts table to data frame
colnames(cell_counts_df) <- c("cell_type", "group", "count") #renames columns for clarity
p3 <- ggplot(cell_counts_df, aes(x = cell_type, y = count, fill = group)) +
             geom_bar(stat = "identity", position = "dodge") +
             labs(title = "Cell Type Abundance in WT vs. Mutant Groups",
             x = "Cell Type",
             y = "Number of Cells",
             fill = "Group") +
             theme_minimal() +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)) #creates a bar plots
png("preleuk_dashboard/plots/cellabundance_plot.png", width = 12, height = 6, units = "in", res = 300) #gives a path for the plot
print(p3) #saves the plot to file
dev.off() #turns off graphical device

gc() #garbage collection

SaveH5Seurat(seurat_integrated, filename = "data/seurat_inter/seurat_integrated.h5Seurat") #saves our seurat object
Convert("data/seurat_inter/seurat_integrated.h5Seurat", dest = "h5ad", overwrite = TRUE) #converts it to format readable by anndata

print("Script Successfully Completed!") #completion print

