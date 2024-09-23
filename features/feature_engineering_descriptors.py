import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import (
    Descriptors, rdMolDescriptors, AllChem 
    )
import networkx as nx



def load_dataset(filename: str) -> pd.DataFrame:
    """
    Загружает датасет из папки 'data' в структуре проекта.

    :param filename: Название файла, который нужно загрузить (например, 'df_cleaned_antibody.csv').
    :return: DataFrame с загруженными данными.
    """
    current_dir = os.getcwd()
    
    # Формируем путь к файлу, который находится в папке 'data'
    data_path = os.path.join(current_dir, 'data', filename)
    
    # Загружаем CSV файл в DataFrame
    try:
        df = pd.read_csv(data_path)
        print(f"Файл успешно загружен: {data_path}")
        return df
    except FileNotFoundError:
        print(f"Файл не найден: {data_path}")
        return None

    
def generate_basic_descriptors(smiles_str: str) -> dict:
    """
    Генерирует основные молекулярные дескрипторы из SMILES строки.
    
    Возвращает словарь с дескрипторами.
    """
    mol = Chem.MolFromSmiles(smiles_str)
    
    descriptors = {
        'mol_MR': Descriptors.MolMR(mol), # Молекулярный рефрактивный индекс
        'logP': Descriptors.MolLogP(mol), # Липофильность (коэффициент распределения LogP)
        'num_h_donors': Descriptors.NumHDonors(mol), # Число доноров водородных связей
        'num_h_acceptors': Descriptors.NumHAcceptors(mol), # Число акцепторов водородных связей
        'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol), # Число вращаемых связей
        'tpsa': Descriptors.TPSA(mol), # Топологическая площадь полярной поверхности (TPSA)
        'fraction_csp3': rdMolDescriptors.CalcFractionCSP3(mol), # Доля атомов углерода в sp3-гибридизации
        'ring_count': mol.GetRingInfo().NumRings(), # Число колец в молекуле
        'aromatic_ring_count': rdMolDescriptors.CalcNumAromaticRings(mol) # Число ароматических колец
    }
    
    return descriptors


def add_descriptors_to_dataset(df: pd.DataFrame, smiles_column: str) -> pd.DataFrame:
    """
    Добавляет дескрипторы для каждого соединения на основе SMILES строки.
    
    Параметры:
    df: pandas DataFrame - исходный датасет.
    smiles_column: str - название столбца с SMILES.

    Возвращает:
    pandas DataFrame - датасет с новыми колонками дескрипторов.
    """
    descriptors_list = []

    for smiles in df[smiles_column]:
        descriptors = generate_basic_descriptors(smiles)
        descriptors_list.append(descriptors)

    descriptors_df = pd.DataFrame(descriptors_list)

    df_with_descriptors = pd.concat([df, descriptors_df], axis=1)

    return df_with_descriptors


def smiles_to_graph(smiles_str: str) -> nx.Graph:
    """
    Преобразует SMILES в граф с использованием RDKit и NetworkX.
    Вершины — это атомы, рёбра — это связи.
    
    Параметры:
    smiles_str: str - строка SMILES, представляющая молекулу.
    
    Возвращает:
    graph: nx.Graph - граф молекулы, где атомы являются вершинами, а связи рёбрами.
    """
    mol = Chem.MolFromSmiles(smiles_str)
    
    graph = nx.Graph()
    
    # Добавляем вершины (атомы)
    for atom in mol.GetAtoms():
        graph.add_node(atom.GetIdx(), element=atom.GetSymbol())
    
    # Добавляем рёбра (связи)
    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondType())
    
    return graph


def generate_graph_metrics_from_graph(graph: nx.Graph) -> dict:
    """
    Генерирует графовые метрики для уже созданного графа: PageRank, Degree Centrality,
    Closeness Centrality и Betweenness Centrality.

    Возвращает словарь с агрегированными метриками для всей молекулы.
    """

    # PageRank: показывает "вес" (важность) атомов в структуре
    pagerank = nx.pagerank(graph)

    # Degree Centrality: показывает, сколько связей имеет каждый атом
    degree_centrality = nx.degree_centrality(graph)

    # Closeness Centrality: указывает, насколько атом "близок" к другим атомам в молекуле
    closeness_centrality = nx.closeness_centrality(graph)

    # Betweenness Centrality: измеряет, насколько часто атом лежит на кратчайших путях между другими атомами
    betweenness_centrality = nx.betweenness_centrality(graph)

    # Агрегация
    metrics = {
        'mean_pagerank': sum(pagerank.values()) / len(pagerank),                   # Среднее значение PageRank по всей молекуле
        'mean_degree_centrality': sum(degree_centrality.values()) / len(degree_centrality),  # Среднее значение Degree Centrality по всей молекуле
        'mean_closeness_centrality': sum(closeness_centrality.values()) / len(closeness_centrality),  # Среднее значение Closeness Centrality
        'mean_betweenness_centrality': sum(betweenness_centrality.values()) / len(betweenness_centrality)  # Среднее значение Betweenness Centrality
    }

    return metrics


def add_graph_metrics_to_dataset(df: pd.DataFrame, smiles_column: str) -> pd.DataFrame:
    """
    Добавляет графовые метрики (PageRank, Degree Centrality, Closeness Centrality, Betweenness Centrality)
    к каждому SMILES в датасете.
    """
    graph_metrics_list = []

    for smiles in df[smiles_column]:
        graph = smiles_to_graph(smiles)
        metrics = generate_graph_metrics_from_graph(graph)
        graph_metrics_list.append(metrics)

    metrics_df = pd.DataFrame(graph_metrics_list)
    
    df_with_metrics = pd.concat([df, metrics_df], axis=1)
    
    return df_with_metrics