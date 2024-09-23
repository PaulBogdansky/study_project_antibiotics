import unittest
import pandas as pd
from features.feature_engineering_descriptors import (
    load_dataset, generate_basic_descriptors, add_descriptors_to_dataset, smiles_to_graph, generate_graph_metrics_from_graph, add_graph_metrics_to_dataset
)

class TestFunctionOfFeatureEngineering(unittest.TestCase):

    def setUp(self):
        self.test_data = {
            'SMILES': [
                'CCO',
                'C1=CC=CC=C1',
                'C(C(=O)O)N' 
            ],
            'bacteria': [
                'Methicillin-resistant S. aureus (ATCC 43300)',
                'methicillin-resistant S. aureus (MRSA)',
                'methicillin-resistant S. aureus (MRSA)'
            ],
            'MIC (µg/mL)': [50.0, 12.5, 6.25],
            'MIC (µmol/mL)': [0.089133, 0.020647, 0.010324],
            'log MIC': [1.950037, 1.314859, 1.013829],
            'S. aureus': ['+', '+', '+'],
            'molar_mass': [46.07, 78.11, 75.07]
        }
        self.df = pd.DataFrame(self.test_data)

    def test_load_dataset(self):
        df_loaded = load_dataset('df_cleaned_antibody.csv')
        self.assertIsInstance(df_loaded, pd.DataFrame)
        self.assertGreater(len(df_loaded), 0)

    def test_generate_basic_descriptors(self):
        smiles = 'CCO'
        descriptors = generate_basic_descriptors(smiles)
        self.assertIn('mol_MR', descriptors)
        self.assertIsInstance(descriptors['mol_MR'], float)

    def test_add_descriptors_to_dataset(self):
        df_with_descriptors = add_descriptors_to_dataset(self.df, 'SMILES')
        self.assertEqual(df_with_descriptors.shape[1], self.df.shape[1] + 9)  # 9 дескрипторов

    def test_smiles_to_graph(self):
        smiles = 'CCO'
        graph = smiles_to_graph(smiles)
        self.assertEqual(len(graph.nodes), 3)  # 3 атома
        self.assertEqual(len(graph.edges), 2)  # 2 связи

    def test_generate_graph_metrics_from_graph(self):
        smiles = 'CCO'
        graph = smiles_to_graph(smiles)
        metrics = generate_graph_metrics_from_graph(graph)
        self.assertIn('mean_pagerank', metrics)
        self.assertIsInstance(metrics['mean_pagerank'], float)

    def test_add_graph_metrics_to_dataset(self):
        df_with_metrics = add_graph_metrics_to_dataset(self.df, 'SMILES')
        self.assertEqual(df_with_metrics.shape[1], self.df.shape[1] + 4)  # 4 метрики

if __name__ == '__main__':
    unittest.main()