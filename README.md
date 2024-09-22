# Project: Analysis of antibiotic database
## Description
This project is designed to analyze and process a specific dataset related to the biological activity of an antibiotic. The code is designed to work exclusively with this data and is not intended to be generalized to other datasets.
## Project structure
**/data/** # Project data  
**/features/** # Code for creating and processing features  
**/models/** # Saved models and their parameters  
**/preprocessing/** # Scripts for preprocessing data  
**/statistics/** # Statistical analysis of data
**/test/** # Testing the model  
**/visualization/** # Visualization of data and results  
**README.md** # Current manual
## Data
The dataset goes through several processing stages, including preprocessing and feature engineering. This section describes each stage of the data.
1. Initial data  
The initial data contains the following columns:  

| Column                               | Description                                                                   |
| ------------------------------------ | ----------------------------------------------------------------------------- |
| molecule (Canonical/Isomeric SMILES) | SMILES code representing the molecule structure (canonical or isomeric)       |
| bacteria                             | The bacterial strain tested against the molecule                              |
| MIC (µmol mL−1)                      | Minimum Inhibitory Concentration (µmol/mL) needed to inhibit bacterial growth |
| MIC (µg mL−1)                        | Minimum Inhibitory Concentration (µg/mL) needed to inhibit bacterial growth   |
| doi                                  | Link to the publication where the data is sourced from                        |
| id                                   | Unique identifier of the molecule                                             |
| type                                 | Type of molecule                                                              |
| group                                | A specific fragment of the molecule (requires clarification)                  |
| group_name                           | A part of the molecule (requires clarification)                               |
| Unnamed: 9                           | Undefined column (requires clarification)                                     |
| Molecular weight                     | Molecular weight of the substance                                             |
| MIC                                  | Minimum Inhibitory Concentration (nmol/mL) needed to inhibit bacterial growth |
| log MIC                              | Logarithm of the MIC value                                                    |
| S. aureus                            | Indicator for S. aureus strain (where '-' indicates exclusion from analysis)  |
| E. coli                              | Indicator for E. coli strain (where '+' indicates inclusion in analysis)      |
| Unnamed: 15                          | Undefined column (requires clarification)                                     |

Example of data entry in CSV format:
O=C(OCC)C1=C(NC(=O)NC1C=2OC=CC2)CN3C=NC=4C=CC=CC43,"E. c. 25922, E. coli25922",,64,10.1002/cjoc.202200326 (https://onlinelibrary.wiley.com/doi/epdf/10.1002/cjoc.202200326),5a,hybrid,O=C1NCC=C(CN2C=NC3=C2C=CC=C3)N1,pyrimidinone,,366.37700000000007,174.68345447448937,2.242251771779038,-,+
 

1. Data after preprocessing

After the preprocessing step, the data can be cleaned, standardizated or transformed.  

| Column        | Description                                                                |
| ------------- | -------------------------------------------------------------------------- |
| SMILES        | Canonical/Isomeric SMILES representation of the molecule                   |
| bacteria      | The bacterial strain tested against the molecule                           |
| MIC (µg/mL)   | Minimum Inhibitory Concentration in micrograms per milliliter (µg/mL)      |
| MIC (µmol/mL) | Minimum Inhibitory Concentration in micromoles per milliliter (µmol/mL)    |
| log MIC       | Logarithm of the MIC value                                                 |
| S. aureus     | Indicator for S. aureus strain (where '+' indicates inclusion in analysis) |
| molar_mass    | Molecular weight of the substance (in g/mol)                               |

1. Data after feature engineering
   
After the feature engineering stage, new features may emerge based on the original data.

| Column                     | Description                                                                 |
| -------------------------- | --------------------------------------------------------------------------- |
| molecule_id                | Unique molecule identifier                                                  |
| molecular_weight_scaled    | Scaled molecular weight                                                     |
| logP_scaled                | Scaled logarithm of the distribution coefficient                            |
| H_bond_donors_binarized    | Binary indicator of whether a molecule has at least one hydrogen bond donor |
| H_bond_acceptors_binarized | Binary feature for hydrogen bond acceptors                                  |
| num_rotatable_bonds        | Total number of rotating links                                              |
| smiles_vector_pca          | Reduced feature vector after PCA                                            |
| new_feature_X              | A new feature built on a combination of data                                |
| target_activity            | The target variable remains unchanged                                       |

## Installation  
## Requirements
- Python 3.x
- Necessary dependencies are listed in requirements.txt  

Installation steps:   
1. Clone the repository:
git clone https://github.com/yourusername/yourproject.git
cd yourproject
2. Create and activate the virtual environment:
python -m venv venv
source venv/bin/activate # for Linux/macOS
venv\Scripts\activate # for Windows
3. Install dependencies:
pip install -r requirements.txt  

## Usage

1. Place the data in the /data/ directory.
2. Run the script to preprocess the data:
python preprocessing/preprocess.py
3. Run the model:
python models/train.py
4. View the visualization results:
python visualization/plot_results.py  

## Restrictions

**Important Note**: The code is designed to work only with the current data set. Using other data is not supported and may lead to incorrect results.  
## Authors

| Names           | Part in project |
| --------------- | --------------- |
| Pavel Bogdanov  |                 |
| Mikhail Budanov |                 |
| Anita Bagramian |                 |
