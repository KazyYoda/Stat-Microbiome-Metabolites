HAllA Installation and Analysis

HAllA (Hierarchical All-against-All association) is used to identify associations between two high-dimensional datasets, such as:

* Gut microbiota and functional profiles
* Gut microbiota and metabolites
* Metabolites and functional profiles

HAllA supports multiple association measures, including Spearman correlation, Pearson correlation, distance correlation, and XICOR. It also applies hierarchical false discovery rate correction for high-dimensional association testing. (Huttenhower Lab)

For the official HAllA documentation and installation instructions, please refer to:

Huttenhower Lab: HAllA (https://huttenhower.sph.harvard.edu/halla/)

⸻

1. Requirements

HAllA requires:

* Python ≥ 3.7
* R ≥ 3.6.1
* HAllA Python package
* R package eva
* R package XICOR

These requirements are based on the official HAllA documentation. (Huttenhower Lab)

For macOS, the following workflow can be performed using the Terminal in RStudio.

⸻

2. Activate the Conda Environment

Open the Terminal in RStudio:

Tools → Terminal → New Terminal

Activate the Conda environment used for the analysis:

conda activate myenv

Check that Python is available:

python --version

For example:

Python 3.11.x

⸻

3. Install HAllA

With myenv activated, install HAllA using pip:

pip install halla

This is the recommended installation method provided by the official HAllA documentation. (Huttenhower Lab)

After installation, check that HAllA is available:

halla --help

If the HAllA help message appears, the installation was successful.

⸻

4. Install Required R Packages

HAllA requires two R packages that should be installed manually:

* eva
* XICOR

Open the R Console in RStudio.

Run:

install.packages("eva")
install.packages("XICOR")

Check that both packages can be loaded:

library(eva)
library(XICOR)

If no error is returned, the required R packages are available.

The official HAllA documentation specifically instructs users to install these R packages manually. (Huttenhower Lab)

⸻

5. Return to the RStudio Terminal

After installing the R packages, return to the Terminal in RStudio.

Activate the Conda environment again if necessary:

conda activate myenv

Check HAllA:

halla --help

⸻

6. Set the Working Directory

Navigate to the directory containing the HAllA input files and analysis scripts.

For example:

cd ~/Documents/HillTribe_NGS/6.HAllA
cd ~/Documents/HillTribe_NGS/4.Functional_profile/HAllA_Genus_Function/Geo

Check the current directory:

pwd

Then check the available files:

ls

⸻

7. Run HAllA Analysis

The general HAllA command is:

halla -X DATASET1 -Y DATASET2 --output OUTPUT_DIRECTORY

For example:

halla -X microbiome.txt -Y metabolite.txt -o HAllA_OUTPUT

HAllA uses two tab-delimited input files. The two datasets must contain the same samples, represented as columns, while the features can differ between the two datasets. Features are represented as rows. (Huttenhower Lab)

Therefore, the expected structure is:

              Sample1  Sample2  Sample3  Sample4
Feature_1        ...      ...      ...      ...
Feature_2        ...      ...      ...      ...
Feature_3        ...      ...      ...      ...

For example:

Microbiome dataset
        │
        │ same samples
        ▼
┌───────────────────┐
│ Sample 1          │
│ Sample 2          │
│ Sample 3          │
│ ...               │
└───────────────────┘
        │
        ▼
      HAllA
        ▲
        │
┌───────────────────┐
│ Metabolite data   │
│ Sample 1          │
│ Sample 2          │
│ Sample 3          │
│ ...               │
└───────────────────┘

⸻

8. Example: Microbiome–Metabolite Association

For example, if the input files are:

microbiome.txt
metabolites.txt

run:

halla -X microbiome.txt -Y metabolites.txt -o HAllA_Microbiome_Metabolite

The results will be saved in:

HAllA_Microbiome_Metabolite/

⸻

9. HAllA Output

HAllA generates an associations.txt file containing the significant associations.

For example:

HAllA_Microbiome_Metabolite/
├── associations.txt
├── ...

The association results include information such as:

* Association rank
* Feature cluster from dataset X
* Feature cluster from dataset Y
* Similarity scores
* p-value
* q-value
* Similarity score between clusters

HAllA also generates complementary visualizations, including hallagrams, diagnostic plots, and heatmaps. (Huttenhower Lab)

⸻

10. Complete Installation Workflow

Terminal in RStudio

conda activate myenv
python --version
pip install halla
halla --help

R Console in RStudio

install.packages("eva")
install.packages("XICOR")
library(eva)
library(XICOR)

Return to Terminal

conda activate myenv
halla --help
cd ~/Documents/HillTribe_NGS/6.HAllA
halla -X microbiome.txt -Y metabolites.txt -o HAllA_OUTPUT

⸻

11. Analysis Workflow

                 HAllA
                   │
        ┌──────────┴──────────┐
        │                     │
        ▼                     ▼
   Dataset X             Dataset Y
   Microbiome            Metabolites
        │                     │
        └──────────┬──────────┘
                   │
                   ▼
          Hierarchical clustering
                   │
                   ▼
          Association testing
                   │
                   ▼
       Multiple-testing correction
                   │
                   ▼
          Significant associations
                   │
        ┌──────────┼──────────┐
        ▼          ▼          ▼
 associations   Hallagram   Heatmaps
    .txt

⸻

12. Important Notes

Same samples are required

The two input datasets must contain the same samples and must be aligned consistently across columns. (Huttenhower Lab)

Before running HAllA, verify:

* Sample IDs are consistent.
* Sample order is consistent.
* Features are rows.
* Samples are columns.
* The datasets do not contain unexpected missing or non-numeric values for the selected analysis.

Conda environment

The myenv environment only needs to be created once.

For subsequent analyses, activate it with:

conda activate myenv

You do not need to reinstall HAllA every time.

R packages

The eva and XICOR packages only need to be installed once in the R environment.

⸻

Quick Reference

Task	Command
Activate Conda environment	conda activate myenv
Check Python	python --version
Install HAllA	pip install halla
Check HAllA	halla --help
Install eva	install.packages("eva")
Install XICOR	install.packages("XICOR")
Set directory	cd ~/Documents/HillTribe_NGS/6.HAllA
Run HAllA	halla -X DATASET1 -Y DATASET2 -o OUTPUT

Official documentation: HAllA – Huttenhower Lab
