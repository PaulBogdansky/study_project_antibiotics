# Project: Analysis of antibiotic database
## Description
This project is designed to analyze and process a specific dataset related to the biological activity of an antibiotic. The code is designed to work exclusively with this data and is not intended to be generalized to other datasets.
## Project structure
**/data/** # Project data  
**/features/** # Code for creating and processing features  
**/models/** # Saved models and their parameters  
**/preprocessing/** # Scripts for preprocessing data  
**/statistics/** # Statistical analysis of data
**/test/**  # Testing the model  
**/visualization/** # Visualization of data and results  
**README.md** # Current manual
## Data
The dataset goes through several processing stages, including preprocessing and feature engineering. This section describes each stage of the data.
1. Initial data  
The initial data contains the following columns:  

| Column 1         | Description                                                                                                 |
| ---------------- | ----------------------------------------------------------------------------------------------------------- |
| molecule_ID      | Molecule or substance unique identifier                                                                     |
| Smiles           | SMILES string describing molecular structure                                                                |
| molecular_weight | Molecular mass of a substance in daltons                                                                    |
| logP             | Logarithm of the partition coefficient of a substance between the aqueous and lipid phases (hydrophobicity) |
| H_bond_donors    | The number of hydrogen bond donors in a molecule                                                            |
| H_bond_acceptors | Number of hydrogen bond acceptors                                                                           |
| rotatable_bonds  | The number of rotating bonds in a molecule                                                                  |
| target_activity  | A target variable representing the activity of a molecule against a specific biological target              |

Example of data entry in CSV format:
molecule_id,smiles,molecular_weight,logP,H_bond_donors,H_bond_acceptors,rotatable_bonds,target_activity
1,C1CCCCC1,84.16,2.29,0,0,0,1.23
2,C1=CC=CC=C1,78.11,1.70,0,0,0,0.95   

2. Data after preprocessing

After the preprocessing step, the data can be cleaned, normalized or transformed.  

| Column 1         | Description                                                        |
| ---------------- | ------------------------------------------------------------------ |
| molecule_id      | Unique identifier of the molecule (unchanged)                      |
| molecular_weight | Molecular weight (can be normalized or scaled)                     |
| logP             | Logarithm of the distribution coefficient (can also be normalized) |
| H_bond_donors    | Number of hydrogen bond donors (after gap processing)              |
| rotatable_bonds  | Number of rotating bonds (after normalization)                     |
| smiles_vector    | A vector representing a molecule in numerical format               |
| target_activity  | The target variable remains unchanged                              |

3. Data after feature engineering
   
After the feature engineering stage, new features may emerge based on the original data.

| Column 1                   | Description                                                                 |
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
