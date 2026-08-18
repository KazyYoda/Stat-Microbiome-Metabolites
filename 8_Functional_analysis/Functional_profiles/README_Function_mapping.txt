Functional Profile Mapping: KO and EC Mapping

This workflow is used to map functional profile identifiers to standardized functional annotations using Python scripts.

Two mapping workflows are available:

* KO mapping using KO_mapping.py
* EC mapping using EC_mapping.py

The scripts can be executed directly from the Terminal in RStudio using the Conda environment created for the functional mapping workflow.

⸻

1. Open the Terminal in RStudio

In RStudio, open:

Tools → Terminal → New Terminal

The commands below should be executed in the RStudio Terminal.

The mapping itself is performed by Python. RStudio is only used as the working environment for accessing the Terminal and managing the analysis project.

⸻

2. Activate the Conda Environment

Activate the Conda environment containing Python and the required packages:

conda activate myenv

The terminal should show the environment name at the beginning of the command line:

(myenv)

⸻

3. Check the Python Version

Confirm that Python is available in the activated environment:

python --version

For example:

Python 3.11.x

If the Python version is displayed correctly, the environment is ready to run the mapping scripts.

⸻

4. Set the Working Directory

Navigate to the directory containing the functional mapping scripts:

cd ~/Documents/HillTribe_NGS/4.Functional_profile/Function_Mapping

Check the contents of the directory:

ls

You should see the mapping scripts and the corresponding input files.

⸻

5. KO Mapping

To perform KEGG Orthology (KO) mapping, run:

python KO_mapping.py

The script will process the input functional profile and generate the mapped KO output.

The output file will be saved in the current working directory.

⸻

6. EC Mapping

To perform Enzyme Commission (EC) mapping, run:

python EC_mapping.py

The script will process the input functional profile and generate the mapped EC output.

The output file will also be saved in the current working directory.

Note: Make sure the filename exactly matches the Python script in your repository. For example, if the script is named EC_mappping.py with three ps, use that exact filename in the command. For consistency, I recommend naming it EC_mapping.py.

⸻

7. Complete Workflow

KO mapping

conda activate myenv
python --version
cd ~/Documents/HillTribe_NGS/4.Functional_profile/Function_Mapping
python KO_mapping.py

EC mapping

conda activate myenv
python --version
cd ~/Documents/HillTribe_NGS/4.Functional_profile/Function_Mapping
python EC_mapping.py

⸻

8. Output

The generated mapping files will be stored in the same directory:

~/Documents/HillTribe_NGS/4.Functional_profile/Function_Mapping/

The general workflow is:

Functional profile
       │
       ├───────────────┐
       │               │
       ▼               ▼
KO_mapping.py      EC_mapping.py
       │               │
       ▼               ▼
   KO mapping       EC mapping
       │               │
       └───────┬───────┘
               ▼
        Mapped output
               │
               ▼
       Downstream analysis

⸻

9. Troubleshooting

Conda environment not found

If you receive:

conda: command not found

Conda may not be initialized in the current Terminal session.

Restart RStudio or initialize Conda according to your Miniconda installation.

Python not found

If:

python --version

does not work, make sure the environment is activated:

conda activate myenv

Then check again:

python --version

Script not found

If you receive:

python: can't open file 'KO_mapping.py'

check that you are in the correct directory:

pwd

Then:

ls

Make sure KO_mapping.py is present in the directory.

⸻

10. Quick Reference

Task	Command
Activate environment	conda activate myenv
Check Python	python --version
Set working directory	cd ~/Documents/HillTribe_NGS/4.Functional_profile/Function_Mapping
KO mapping	python KO_mapping.py
EC mapping	python EC_mapping.py
Check current directory	pwd
List files	ls

In short:

RStudio Terminal
      │
      ▼
conda activate myenv
      │
      ▼
cd Function_Mapping/
      │
      ├── python KO_mapping.py
      │
      └── python EC_mapping.py
              │
              ▼
        Mapping output
