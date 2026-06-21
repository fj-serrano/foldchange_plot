# Fold Change Plot (Python)

This Python script processes differential expression matrices from multiple studies and generates four graphical representations: a bar plot, concentration distribution plot, MA plot and volcano plot. The resulting figures support expression analysis and quality assessment across datasets.

---

## Input files

- **Directory containing sRNAde results (differential expression matrices)**
  - Selected by the user 
  - The folder for the assignment method must be named: 
	-	de_rcsa 
	-	de_rcadj
  - The matrices must be named: 
	-	`DESeq2_[Condition 1]_vs_[Condition 2].tsv`
	-	`edgeR_[Condition 1]_vs_[Condition 2].tsv`
	-	`limma_[Condition 1]_vs_[Condition 2].tsv`

---

## Requirements

- **Python 3**
- Required libraries:
  - `pandas`
  - `matplotlib`
  - `numpy`
  - `adjustText`
  - `seaborn`

No additional external dependencies are required.

---

## Script usage

Run the script from the terminal or a Python environment:

```bash
python foldchange_plot.py -i <input_file> -o <output_file> [-s <study>] [-a <assignment>]
````

Arguments
| Argument | Alternative_Argument | Required | Description | Default |
| --- | --- | --- | --- | --- |
| -i | --input | Yes | Path to the input file | - |
| -o | --output | Yes | Path to the output file | - |
| -s | --study | No | Select the name of the study | "SRP" |
| -a | --assignment | No | Select the assignment | "de_rcadj" |

---

## Output

The script generates five output files for every study - assignment - differential expression framework:

1. **Tab-separated values file (tsv)**
- Contains the processed expression data and the values required for downstream analysis.

2. **Bar plot (png)**
- Bar chart showing the percentage of significant miRNAs per log2FoldChange interval they fall into.

3. **Concentration distribution plot (png)**
- Distribution plot showing the base-10 logarithm of relative expression in one condition against the base-10 logarithm of relative expression in the other condition.

4. **MA plot (png)**
- Plots the base-2 logarithm of the mean expression of a miRNA across the two studied conditions against its log2FoldChange value.

5. **Volcano plot (png)**
- Plots negative log10 of adjusted p-value (padj) against log2FoldChange to identify upregulated and downregulated miRNAs, highlighting the five most significant ones.

---

## How the script works

The script is divided into three main stages:

### 1️. File selection and loading

- Searches the user-selected folder for all differential expression matrices corresponding to the selected assignment method and stores them in a list.
- Each file is iterated over individually and the corresponding `tsv` file is read into a `DataFrame`.  

### 2. Graphical representations

- From the generated `DataFrame`, all the graphical representations described above are created and saved in the selected output directory.

### 3. Tab-separated values file

- A copy of the expression matrix is saved, keeping only significantly and differentially expressed miRNAs.
