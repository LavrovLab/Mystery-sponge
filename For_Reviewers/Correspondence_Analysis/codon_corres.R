library(ade4)
library(coRdon)
library(Biostrings)
library(seqinr)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# Make codon frequency matrix for cds in MysSpon Mt DNA ----
Myscds = readSet(file="MysSponGenes.txt")
Hplakcds = readSet(file="HplakGenes.txt")
Xdcds = readSet(file="XdGenes.txt")
Psymcds = readSet(file="PsymGenes.txt")

Mys_cox2 = readSet(file="MysSpon_cox2.txt")
Hplak_cox2 = readSet(file="hplak_cox2.txt")
Xd_cox2 = readSet(file="xd_cox2.txt")
Psym_cox2 = readSet(file="Psym_cox2.txt")

Mysctable = codonTable(Myscds)
HpCtable = codonTable(Hplakcds)
XdCtable = codonTable(Xdcds)
PsymCtable = codonTable(Psymcds)

Mys_cox2ctable = codonTable(Mys_cox2)
Hplak_cox2ctable = codonTable(Hplak_cox2)
Xd_cox2ctable = codonTable(Xd_cox2)
Psym_cox2ctable = codonTable(Psym_cox2)

Mysfq_Mat = Mysctable@counts
HpFqMAt = HpCtable@counts
XdFqMat = XdCtable@counts
PsymFqMat = PsymCtable@counts

Myscox2fqmat = Mys_cox2ctable@counts
Hpcox2fqmat = Hplak_cox2ctable@counts
Xdcox2fqmat = Xd_cox2ctable@counts
Psymcox2fqmat = Psym_cox2ctable@counts

rownames(Mysfq_Mat) = Mysctable@ID
rownames(HpFqMAt) = HpCtable@ID
rownames(XdFqMat) = XdCtable@ID
rownames(PsymFqMat) = PsymCtable@ID
total_codons = c(sum(Mysfq_Mat), sum(HpFqMAt), sum(XdFqMat), sum(PsymFqMat))

rownames(Myscox2fqmat) = Mys_cox2ctable@ID
rownames(Hpcox2fqmat) = Hplak_cox2ctable@ID
rownames(Xdcox2fqmat) = Xd_cox2ctable@ID
rownames(Psymcox2fqmat) = Psym_cox2ctable@ID

# Codon position base frequency ------
# List of codon frequency tables (example: fq_Mat1, fq_Mat2, fq_Mat3)
codon_tables <- list(Mysfq_Mat, HpFqMAt, XdFqMat, PsymFqMat)  
codon_tables_cox2 = list(Myscox2fqmat,Hpcox2fqmat,Xdcox2fqmat,Psymcox2fqmat)
# Number of tables
num_tables <- length(codon_tables)
num_tables_cox2 = length(codon_tables_cox2)
# Initialize a list to store results for each table
total_usage_list <- vector("list", num_tables_cox2)
names(total_usage_list) = c("H. sp nov", "Hplak", "Xd", "Psym")

# Loop through each codon frequency table
for (t in 1:num_tables_cox2) {
  
  # Convert frequency matrix to summed codon counts
  cdnSum <- as.matrix(colSums(codon_tables_cox2[[t]]))
  
  # Initialize total_position_n matrix for this table
  total_position_n <- matrix(data = 0, nrow = 4, ncol = 3)
  rownames(total_position_n) <- c("A", "G", "C", "T")
  colnames(total_position_n) <- c("1", "2", "3")
  
  # Function to calculate total usage for a given base at a given position
  get_total_usage <- function(base, position) {
    pattern <- paste0("^", strrep(".", position - 1), base)  # Construct regex pattern
    rows <- grep(pattern, rownames(cdnSum))  # Find matching rows
    sum(cdnSum[rows, ])  # Sum the values
  }
  
  # Loop through bases and positions to fill the matrix
  bases <- c("A", "G", "C", "T")
  
  for (i in 1:4) {  # Loop over row indices (bases)
    for (j in 1:3) {  # Loop over column indices (positions)
      total_position_n[i, j] <- get_total_usage(bases[i], j)
    }
  }
  
  # Store the result in the list
  total_usage_list[[t]] <- total_position_n
}

# Print the results
print(total_usage_list)
total_usage_list <- lapply(total_usage_list, function(mat) {
  mat / (sum(mat)/3)  # Divide each cell by the total sum of the dataset
})
print(total_usage_list)


num_datasets <- length(total_usage_list)   # Count how many datasets
num_bases <- nrow(total_usage_list[[1]])  # Number of bases (A, G, C, T)
num_positions <- ncol(total_usage_list[[1]])  # Number of codon positions (1, 2, 3)

# Step 2: Set Up Plot Layout (2x2 Grid)

par(mar = c(2, 2, 2, 1))
par(mfrow = c(2, 2))  # 2x2 grid of plots

# Step 3: Define Labels and Colors

bases <- rownames(total_usage_list[[1]])  # ["A", "G", "C", "T"]
positions <- colnames(total_usage_list[[1]])  # ["1", "2", "3"]
colors <- c("#A6CEE3", "#B2DF8A","#FDBF6F", "#CAB2D6")  
# Use dataset names if available; otherwise, default to "Dataset 1, Dataset 2, ..."
dataset_names <- if (!is.null(names(total_usage_list))) names(total_usage_list) else paste("Dataset", 1:num_datasets)

# Find the maximum y-axis value across all datasets for consistent scaling
ymax <- max(sapply(total_usage_list, max)) * 1.2  

# Step 4: Loop Through Each Base (A, G, C, T) and Plot

for (i in 1:num_bases) {
  
  # Collect data for the current base across all datasets
  values <- sapply(total_usage_list, function(mat) mat[i, ])
  
  # Transpose to align bars correctly for positions
  values <- t(values)
  
  # Create the bar plot
  barplot(values, 
          beside=TRUE, 
          col=colors, 
          ylim=c(0, ymax),
          names.arg=positions,
          main=paste("Base", bases[i]),
          xlab="Codon Position",
          ylab="Frequency")
  
  # Add a legend only in the first plot
  if (i == 1) {
    legend("topright",
           legend = dataset_names,
           fill = colors,
           bty = "n",
           cex = 0.7,               # Shrinks text size
           ncol = 2)                # Arrange legend in 2 columns
  }
  
}

par(mfrow = c(1,1))
#_________________________________________________________________________________________________________________

# correspondence analysis for MysSpon mt cds ----
codon_ca = dudi.coa(Mysfq_Mat, scannf = FALSE, nf = 2)

s.label(
        codon_ca$li, 
        sub = "Gene Distribution in codon CA Space", 
        clabel = 0, 
        cpoint = 1.5, 
        xax = 1, 
        yax = 2, 
        xlim = c(-0.4,1), 
        ylim = c(-0.6,0.4), 
        addaxes = TRUE
        )

text(
      codon_ca$li[, 1], codon_ca$li[, 2] +0.01,
      label = rownames(codon_ca$li),
      cex = 0.8,
      col = "black",
     )

axis(1, at = seq(-1, 1, by = 0.1), labels = seq(-1, 1, by = 0.1))
axis(2, at = seq(-1, 1, by = 0.1), labels = seq(-1, 1, by = 0.1))

summary(codon_ca)
syncuo = t(t(SCUO(cTobject = Mysctable, id_or_name2 = "4", alt.init = TRUE)))
rownames(syncuo) = t(t(Mysctable@ID))

# ----

# CA Analysis of amino acid frequency  ----

mtAA_full = readAAMultipleAlignment(file="Haplo_aa.aln", "FASTA")
AAfqTable = alphabetFrequency(mtAA_full)
rownames(AAfqTable) = names(mtAA_full@unmasked)
AAfqTable = AAfqTable[1:46,1:20]
clades = read_excel("Clades.xlsx")

mtAA_ca = dudi.coa(AAfqTable,scannf = FALSE, nf = 2)

plot_coords = as.data.frame(mtAA_ca$li)
clade_fac = as.factor(clades$tags)
plot_coords <- plot_coords[match(clades$ids, row.names(plot_coords)), ]
point_colors <- setNames(
                          c(
                            "#E57373","red","#64B5F6","#81C784","#BA68C8",
                            "#4DB6AC","#FF8A65","#A1887F","#90A4AE"
                            ),
                          unique(clade_fac)
                         )  # Custom colors

s.class(dfxy = plot_coords, fac = t(t(clade_fac)),
        col = point_colors,
        cstar = 0,
        sub = "mtCDS in CA Space",
        clabel = 0,
        xax = 1, yax = 2,
        xlim = c(-.6,.3),
        cpoint=1
        
        
        )
title(xlab = "CA Axis 1", ylab = "CA Axis 2")
# Add numerical axis ticks and labels
axis(1, at = seq(-1, 1, by = 0.1), labels = seq(-1, 1, by = 0.1))
axis(2, at = seq(-1, 1, by = 0.1), labels = seq(-1, 1, by = 0.1))


legend("topleft",
       legend = levels(clade_fac),
       fill = point_colors,
       title = "Clades",
       ncol = 2,
       cex = 0.6,
       pt.cex = .2,
       inset = c(0.2,0.25),
       border = NA,
       bty = "n"
      )

# Add labels below each point
#text(mtAA_ca$li[, 1] + sample(10,1)/1000, 
#     mtAA_ca$li[, 2] - sample(10,1)/1000,       # Dynamic offset below the dots
 #    labels = rownames(mtAA_ca$li),    # Use row names or a custom vector for labels
 #    cex = 0.6,                        # Adjust size of the text
 #    col = "black")   



# Extract and plot column scores
c1_scores <- mtAA_ca$c1

# Plot the contributions to the first two axes
plot(c1_scores[, 1], c1_scores[, 2], type = "n", xlab = "CA Axis 1", ylab = "CA Axis 2")
text(c1_scores[, 1], c1_scores[, 2], labels = rownames(c1_scores), col = "blue")
abline(h = 0, v = 0, col = "grey", lty = 2)

