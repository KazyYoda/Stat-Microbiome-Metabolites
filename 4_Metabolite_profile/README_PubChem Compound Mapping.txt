PubChem Compound Mapping

This workflow maps a list of compound names to PubChem identifiers and chemical information using the PubChem Python API.

The script retrieves:

* CID: PubChem Compound ID
* IUPAC Name: Systematic chemical name
* InChIKey: Standardized chemical structure identifier

This workflow is useful for metabolomics studies where compound names need to be converted into standardized chemical identifiers for downstream annotation, pathway analysis, database integration, or multi-omics analysis.

⸻

1. Requirements

The following files are required:

Compound_Mapping/
├── Pubchem_mapping.py
└── compound_list.txt

Input file

The input should be a plain-text .txt file containing a list of compound names.

Example:

compound_list.txt
Glucose
Citric acid
Lactic acid
Phenylpyruvic acid
Docosane

Each compound should be listed on a separate line.

Python script

The mapping script:

Pubchem_mapping.py

The script reads the compound names, queries PubChem, and generates a mapped output file.

⸻

2. Install Conda

Miniconda is sufficient for this workflow (https://continuumio-docs.readthedocs-hosted.com/miniconda/).

If Conda is already installed on your computer, you can skip this step.

Check whether Conda is available:

conda --version

If Conda is installed, the terminal should return something similar to:

conda 25.x.x

⸻

3. Create a Conda Environment

It is recommended to create a separate Conda environment for this workflow rather than 
installing Python packages into the system Python installation.

Create an environment called myenv:

conda create -n myenv python=3.11 (or latest version)

Activate the environment:

conda activate myenv

Check the Python version:

python --version

You should see something similar to:

Python 3.11.x

Why use a Conda environment?

The environment keeps the Python version and required packages isolated from the rest of your computer. 
This makes the workflow easier to reproduce and reduces the risk of conflicts with other Python projects.

⸻

4. Install Required Python Packages

After activating myenv, install the packages required by the script.

For example:

pip install pubchempy

If the script uses additional packages, install them in the same environment.

You can check installed packages with:

pip list

⸻

5. Set the Working Directory

Open Terminal and navigate to the directory containing the Python script and input file.

For example:

cd ~/Documents/HillTribe_NGS/5.Metabolites/Compound_Mapping

Check the files in the directory:

ls

You should see something similar to:

Pubchem_mapping.py
compound_list.txt

⸻

6. Activate the Conda Environment

Before running the Python script, activate the environment:

conda activate myenv

Confirm that the correct Python installation is being used:

python --version

You can also check which Python executable is being used:

which python

On macOS/Linux, the path should point to your Conda environment, for example:

.../miniconda3/envs/myenv/bin/python

⸻

7. Configure the Input and Output Files

Open:

Pubchem_mapping.py

At the bottom of the script, specify the input and output filenames:

if __name__ == "__main__":
    input_file = "compound_list.txt"
    output_file = "Mapped_PubChem_output.txt"

Input file

Change:

input_file = "compound_list.txt"

to the name of your compound list.

For example:

input_file = "compound_list.txt"

Output file

Change:

output_file = "Mapped_PubChem_out_put.txt"

to your preferred output filename.

For example:

output_file = "Mapped_PubChem_metabolites.txt"

Because the input and output files are specified using filenames only, they should normally be located in the same working directory as the Python script.

⸻

8. Run the Mapping

From the terminal (via R Studio), run:

python Pubchem_mapping.py

The script will read the compound names from the input .txt file and query PubChem for the corresponding compound information.

Example:

(myenv) $ python Pubchem_mapping.py

⸻

9. Output

The mapped output file will be saved in the working directory.

For example:

Compound_Mapping/
├── Pubchem_mapping.py
├── compound_list.txt
└── Mapped_PubChem_output.txt

The output contains the mapped PubChem information, including:

Field	Description
Compound Name	Name provided in the input file
CID	PubChem Compound ID
IUPAC Name	Systematic chemical name
InChIKey	Standardized chemical structure identifier

The resulting file can then be imported into R or other downstream analysis software.

⸻

10. Complete Workflow

For a new user, the complete workflow is:

First-time setup

conda create -n myenv python=3.11
conda activate myenv
pip install pubchempy

Run the analysis

cd ~/Documents/HillTribe_NGS/5.Metabolites/Compound_Mapping
conda activate myenv
python --version
python Pubchem_mapping.py

The mapped output will be generated in the same directory.

⸻

11. Using the Workflow for a New Compound List

To reuse this workflow for another metabolomics dataset:

1. Prepare a .txt file containing compound names.
2. Place the file in the Compound_Mapping directory.
3. Open Pubchem_mapping.py
4. Change input_file.
5. Change output_file.
6. Activate the Conda environment in the terminal via R Studio.
7. Run the Python script.

For example:

if __name__ == "__main__":
    input_file = "new_compounds.txt"
    output_file = "Mapped_PubChem_new_compounds.txt"

Then run:

python Pubchem_mapping.py

⸻

12. Notes

Conda environment

The Conda environment only needs to be created once.

After that, you only need to activate it before running the script:

conda activate myenv

You do not need to recreate the environment every time.

Python installation

Python should be installed inside the Conda environment:

conda create -n myenv python=3.11

This is preferable to installing Python separately on the system.

Package installation

If pubchempy has already been installed in myenv, you do not need to reinstall it each time.

⸻

14. Summary

Compound names (.txt)
        │
        ▼
Pubchem_mapping.py
        │
        ▼
PubChem query
        │
        ├── CID
        ├── IUPAC Name
        └── InChIKey
        │
        ▼
Mapped_PubChem_*.txt
        │
        ▼
Downstream metabolomics analysis

Main command:

python Pubchem_mapping.py
