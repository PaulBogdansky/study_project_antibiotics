import pandas as pd
import numpy as np
import re
from rdkit import Chem
from rdkit.Chem import Descriptors
import os
from sklearn.preprocessing import StandardScaler
from transliterate import translit


def find_project_root(start_dir: str, project_name: str) -> str:
    """
    Ищет корневую директорию проекта.

    :param start_dir: Стартовая директория, откуда начинается поиск (тип str).
    :param project_name: Имя папки проекта, которую нужно найти (тип str).
    :return: Путь к корневой директории проекта (тип str).
    :raises FileNotFoundError: Если проектная папка не найдена.
    """
    current_dir = start_dir
    while True:
        if project_name in os.listdir(current_dir):
            return os.path.join(current_dir, project_name)
        new_dir = os.path.dirname(current_dir)
        if new_dir == current_dir:
            raise FileNotFoundError(f"Папка проекта '{project_name}' не найдена.")
        current_dir = new_dir


def load_dataset_from_project(
    project_name: str, relative_file_path: str
) -> (
    pd.DataFrame
):  #  на вход строка, на выход df. И так для каждой ф-ции (что принимает и что выдает)
    """
    Загружает CSV файл из директории проекта.

    :param project_name: Имя папки проекта (тип str).
    :param relative_file_path: Относительный путь к файлу внутри проекта (тип str).
    :return: DataFrame с загруженными данными (тип pd.DataFrame).
    :raises FileNotFoundError: Если файл не найден по указанному пути.
    """
    start_directory = os.getcwd()

    # Ищем корневую папку проекта
    project_root = find_project_root(start_directory, project_name)

    # Полный путь к файлу
    file_path = os.path.join(project_root, relative_file_path)

    # Загружаем файл
    try:
        df = pd.read_csv(file_path)
        print(f"Файл успешно загружен: {file_path}")
        return df
    except FileNotFoundError:
        print(f"Файл не найден по указанному пути: {file_path}")
        return None


def clean_dataset(df: pd.DataFrame, columns_to_drop: list) -> pd.DataFrame:
    """
    Очищает датасет: удаляет ненужные колонки, строки с пропусками и переименовывает колонки.

    :param df: Входной DataFrame (тип pd.DataFrame).
    :param columns_to_drop: Список колонок, которые нужно удалить (тип list).
    :return: Очищенный DataFrame (тип pd.DataFrame).
    """
    df_copy = df.copy()

    # Удаляем указанные колонки
    df_cleaned = df_copy.drop(columns_to_drop, axis=1)

    # Удаляем строки с пустыми или отсутствующими значениями
    df_cleaned = df_cleaned[df_cleaned["S. aureus"] != "-"].reset_index(drop=True)
    df_cleaned = df_cleaned.dropna(
        subset=["molecule (Canonical/Isomeric SMILES)", "S. aureus"]
    ).reset_index(drop=True)

    # Переименовываем столбцы
    df_cleaned = df_cleaned.rename(
        columns={
            "molecule (Canonical/Isomeric SMILES)": "SMILES",
            "MIC": "MIC (µmol/mL)",
            "MIC (µg mL−1)": "MIC (µg/mL)",
        }
    )

    return df_cleaned


def contains_cyrillic(text: str) -> bool:
    """
    Проверяет наличие кириллических символов в тексте.

    :param text: Текст, который нужно проверить (тип str).
    :return: True, если найдены кириллические символы, иначе False (тип bool).
    """
    return bool(re.search("[а-яА-Я]", str(text)))


def detect_cyrillic_in_dataframe(df: pd.DataFrame):
    """
    Проверяет весь DataFrame на наличие кириллических символов и выводит их местоположение.

    :param df: Входной DataFrame (тип pd.DataFrame).
    """
    found_cyrillic = False  # Переменная для отслеживания наличия кириллицы

    for column in df.columns:
        for idx, value in df[column].items():
            if contains_cyrillic(value):
                found_cyrillic = True
                print(
                    f"Кириллица найдена в колонке '{column}', строка {idx}, значение: {value}"
                )

    if not found_cyrillic:
        print("Кириллица не найдена.")


def replace_cyrillic(text: str) -> str:
    """
    Заменяет кириллические символы на латинские с использованием транслитерации.

    :param text: Текст для замены (тип str).
    :return: Текст с заменённой кириллицей на латиницу (тип str).
    """
    if contains_cyrillic(text):
        return translit(str(text), "ru", reversed=True)
    return text


def replace_cyrillic_in_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Применяет замену кириллических символов на латиницу ко всему DataFrame.

    :param df: Входной DataFrame (тип pd.DataFrame).
    :return: DataFrame с заменённой кириллицей на латиницу (тип pd.DataFrame).
    """
    return df.apply(lambda col: col.map(replace_cyrillic))


def calculate_molar_mass(smiles: str) -> float:
    """
    Рассчитывает молекулярную массу по SMILES.

    :param smiles: Строка SMILES (тип str).
    :return: Молекулярная масса (тип float) или None, если SMILES некорректен.
    """
    if isinstance(smiles, str):  # Проверяем, является ли SMILES строкой
        molecule = Chem.MolFromSmiles(smiles)  # Преобразуем SMILES в молекулу
        if molecule:  # Если молекула корректная
            return Descriptors.MolWt(molecule)  # Возвращаем молекулярную массу
        else:
            return None  # Если SMILES некорректен
    return None  # Если SMILES не является строкой


def validate_and_correct_MIC(df: pd.DataFrame) -> pd.DataFrame:
    """
    Проверяет корректность MIC (µg/mL) и пересчитывает его при необходимости.

    :param df: Входной DataFrame (тип pd.DataFrame).
    :return: DataFrame с проверенными и откорректированными значениями MIC (µg/mL) (тип pd.DataFrame).
    """
    # Проверяем расхождения в колонке 'Molecular weight' и удаляем её, если они есть
    if (df["Molecular weight"] != df["molar_mass"]).any():
        df = df.drop(columns=["Molecular weight"])

    # Преобразуем MIC (µmol/mL) в mL (делим на 1000)
    df["MIC (µmol/mL)"] = df["MIC (µmol/mL)"] / 1000

    # Рассчитываем MIC (µg/mL) на основе молекулярной массы
    MIC_mcg_ml = (df["molar_mass"] * df["MIC (µmol/mL)"]).round(2)

    # Проверяем расхождения между 'MIC (µg/mL)' и рассчитанными значениями, заменяем при необходимости
    mask = df["MIC (µg/mL)"] != MIC_mcg_ml
    df.loc[mask, "MIC (µg/mL)"] = MIC_mcg_ml[mask]

    return df


def standardize_features(df: pd.DataFrame, features: list) -> pd.DataFrame:
    """
    Стандартизирует указанные числовые признаки.

    :param df: Входной DataFrame (тип pd.DataFrame).
    :param features: Список числовых признаков для стандартизации (тип list).
    :return: DataFrame с стандартизированными признаками (тип pd.DataFrame).
    """
    # Создаём объект стандартизатора
    scaler = StandardScaler()

    # Применяем стандартизацию к числовым признакам
    df[features] = scaler.fit_transform(df[features])

    return df
