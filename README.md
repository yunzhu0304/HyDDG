# HyDDG User‘s Guide  <a href="[https://jbengler.github.io/tidyplots/](https://maze-icicle-277.notion.site/HyDDG-User-s-Guide-175b4db0634380048bd4c28d39f0ec05)"><img src="man/figure/logo.svg" align="right" height="139" alt="HyDDG website" /></a>

# HyDDG: Hypergeometric Distribution for Detecting Differentially Expressed Genes (DEGs)  




## 📚Overview


HyDDG is an R package designed for researchers in bioinformatics and computational biology to detect differentially expressed genes (DEGs) using a hypergeometric distribution-based approach. The package is tailored for analyzing transcriptomics data and provides statistical measures for identifying up- and down-regulated genes with precise p-value adjustments and fold change calculations.

---

## 🏗️Installation

To install and load the HyDDG package, use the following commands:

```r
# Install the HyDDG package from GitHub 

if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("yunzhu0304/HyDDG")

# Load the package

library(HyDDG)
```

---

## 📌Input Requirements

HyDDG requires the following inputs:

1. **Data Matrix (`data`)**: A normalized data frame with rows representing probes/genes and columns representing samples.
2. **Group List (`group.list`)**: A factor variable containing two levels, representing control and treatment groups, corresponding to the columns of the data matrix.

Ensure that the data matrix is normalized prior to analysis and that the `group.list` levels are ordered with the control group as the first level and the treatment group as the second level.

---

## 📏Main Functions

### 📍1. `hyfit`

The `hyfit` function calculates fold changes for each sample in the treatment group relative to the mean expression in the control group.

### Example Usage:

```r
# Example data from the CLL package

library(CLL)
data("CLLbatch")
CLLrma <- rma(CLLbatch)
ourData <- as.data.frame(exprs(CLLrma))[c(1:100),]

# Define group list

group <- factor(c(rep("Control", 12), rep("Treat", 12)), levels = c("Control", "Treat"))

# Calculate fold changes

fit <- hyfit(data = ourData, group.list = group)
```

### Output:

A data frame with rows representing genes and columns representing treatment samples, where each value is the fold change relative to the control group mean.

---

### 📊2. `HyDDG`

The `HyDDG` function performs hypergeometric distribution-based statistical analysis to identify DEGs.

### Parameters:

- `data`: Normalized data matrix (same as `hyfit`).
- `group.list`: Factor variable (same as `hyfit`).
- `adj.p.method`: Method for adjusting p-values (default: “BH”). Supported methods include “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, and “none”.
- `BV`: Boundary value for identifying up- or down-regulated genes (default: 1).

### Example Usage:

```r
# Perform DEG analysis
HyDDGResult <- HyDDG(data = ourData, group.list = group)
# View results
head(HyDDGResult)
```

### Output:

`HyDDGResult` is a data frame containing the following columns:
- `ID`: Gene/probe ID.
- `ave.expr`: Average expression of the gene across all samples.
- `lg.p`: Log-transformed p-value.
- `p.value`: Raw p-value for differential expression.
- `p.adj`: Adjusted p-value using the specified method.
- `FC`: Fold change (mean treatment expression / mean control expression).
- `logFC`: Log2-transformed fold change.

---

## 🧭Recommended Thresholds

To identify significant DEGs, the following thresholds are recommended:
- Absolute log-transformed p-value (`abs(lg.p)`）≥3
- Adjusted p-value (`p.adj`) < 0.05

Example:

```r
significant_DEGs <- subset(HyDDGResult, abs(lg.p) >= 3 & p.adj < 0.05)
head(significant_DEGs)
```

---

## 📝Example Workflow

```r
# Step 1: Load required packages
library(CLL)
data("CLLbatch")
CLLrma <- rma(CLLbatch)
ourData <- as.data.frame(exprs(CLLrma))[c(1:100),]

# Step 2: Define group list
group <- factor(c(rep("Control", 12), rep("Treat", 12)), levels = c("Control", "Treat"))

# Step 3: Perform analysis
HyDDGResult <- HyDDG(data = ourData, group.list = group)

# Step 4: View significant results
significant_DEGs <- subset(HyDDGResult, abs(lg.p) >= 3 & p.adj < 0.05)
print(significant_DEGs)
```

---

## ⚠️Notes

- Ensure the input data is normalized (e.g., using `rma` or similar methods).
- Use appropriate thresholds for filtering DEGs based on the experimental context.

---

## 📖References

- Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
Nucleic Acids Research 43(7), e47.
- Kondrakhin YV, Sharipov RN, Keld AE, Kolpakov FA. Identification of differentially expressed genes by meta-analysis of microarray data on breast cancer. In Silico Biol. 2008;8(5-6):383-411. PMID: 19374127.

---

## ☎️Contact

For questions or issues, please contact the package maintainer at Zhu.Yun@mh-hannover.de；Xu.Yanzhe@mh-hannover.de or visit the GitHub repository: [https://github.com/yunzhu0304/HyDDG](https://github.com/yunzhu0304/HyDDG).
