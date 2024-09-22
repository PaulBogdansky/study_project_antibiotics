import unittest
import pandas as pd
from preprocessing.preprocessing import (
    find_project_root, load_dataset_from_project, clean_dataset, contains_cyrillic,
    replace_cyrillic, calculate_molar_mass, validate_and_correct_MIC, standardize_features
)
import os
from unittest.mock import patch, MagicMock

class TestFunctionOfPreprocessing(unittest.TestCase):

    def test_find_project_root(self):
        with patch('os.listdir', return_value=['study_project_antibiotics']):
            root_dir = find_project_root(os.getcwd(), 'study_project_antibiotics')
            self.assertTrue(root_dir.endswith('study_project_antibiotics'))

    def test_load_dataset_from_project(self):
        df = pd.DataFrame({
            'SMILES': ['CCO', 'CC', 'COC'],
            'S. aureus': ['+', '-', '+'],
            'MIC (µmol/mL)': [10, 20, 30]
        })
        with patch('pandas.read_csv', return_value=df):
            dataset = load_dataset_from_project('study_project_antibiotics', 'data/test.csv')
            self.assertEqual(dataset.shape, (3, 3))  # проверяем, что загружены 3 строки и 3 столбца

    def test_clean_dataset(self):
        df = pd.DataFrame({
            'molecule (Canonical/Isomeric SMILES)': ['CCO', 'CC', 'COC'],
            'S. aureus': ['+', '-', '+'],
            'MIC': [10, 20, 30],
            'MIC (µg mL−1)': [100, 200, 300],
            'id': [1, 2, 3]
        })
        columns_to_drop = ['id']
        cleaned_df = clean_dataset(df, columns_to_drop)
        self.assertNotIn('id', cleaned_df.columns)
        self.assertEqual(cleaned_df.shape[0], 2)  # проверяем, что удалена одна строка (где 'S. aureus' == '-')

    def test_contains_cyrillic(self):
        self.assertTrue(contains_cyrillic("Привет"))
        self.assertFalse(contains_cyrillic("Hello"))

    def test_replace_cyrillic(self):
        self.assertEqual(replace_cyrillic("Привет"), "Privet")
        self.assertEqual(replace_cyrillic("Hello"), "Hello")

    def test_calculate_molar_mass(self):
        self.assertEqual(round(calculate_molar_mass("CCO"), 2), 46.07)  # проверяем молекулярную массу для этанола
        self.assertIsNone(calculate_molar_mass("invalid_smiles"))  # проверяем некорректный SMILES
    
    def test_validate_and_correct_MIC(self):
        df = pd.DataFrame({
            'MIC (µmol/mL)': [1, 2],
            'molar_mass': [46.07, 18.02],
            'MIC (µg/mL)': [46.07, 36.04],
            'Molecular weight': [50.00, 20.00]  # Значения отличаются от molar_mass
        })
        corrected_df = validate_and_correct_MIC(df)
        
        # Проверяем, что колонка была удалена
        self.assertNotIn('Molecular weight', corrected_df.columns)
        
        # Проверяем, что значения MIC (µg/mL) остались корректными
        self.assertEqual(corrected_df['MIC (µg/mL)'].tolist(), [0.05, 0.04])


    def test_standardize_features(self):
        df = pd.DataFrame({
            'MIC (µmol/mL)': [1, 2, 3],
            'log MIC': [0.1, 0.2, 0.3],
            'molar_mass': [46.07, 18.02, 78.11],
            'MIC (µg/mL)': [46.07, 36.04, 78.11]
        })
        numerical_features = ['MIC (µmol/mL)', 'log MIC', 'molar_mass', 'MIC (µg/mL)']
        standardized_df = standardize_features(df, numerical_features)
        self.assertAlmostEqual(standardized_df['MIC (µmol/mL)'].mean(), 0)
        self.assertAlmostEqual(standardized_df['molar_mass'].std(), 1, delta=0.25)

if __name__ == '__main__':
    unittest.main()