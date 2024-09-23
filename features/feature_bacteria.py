import pandas as pd


def generate_bacteria_features(bacteria_str: str) -> dict:
    """
    Генерирует фичи на основе штамма S. aureus.

    Параметры:
    bacteria_str: str - строка с описанием штамма

    Возвращает:
    dict - словарь с бинарными и категориальными фичами. None заменены на подходящие значения.
    """
    features = {}

    # Бинарные признаки: MRSA, MSSA
    features["is_mrsa"] = 1 if "MRSA" in bacteria_str else 0
    features["is_mssa"] = 1 if "MSSA" in bacteria_str else 0

    # Бинарные признаки: ATCC, MTCC, CIP коллекции
    features["is_atcc"] = 1 if "ATCC" in bacteria_str else 0
    features["is_mtcc"] = 1 if "MTCC" in bacteria_str else 0
    features["is_cip"] = 1 if "CIP" in bacteria_str else 0

    return features


def add_bacteria_features_to_dataset(
    df: pd.DataFrame, bacteria_column: str
) -> pd.DataFrame:
    """
    Добавляет фичи на основе штаммов бактерий к каждому элементу в датасете.

    Параметры:
    df: pd.DataFrame - исходный DataFrame.
    bacteria_column: str - название столбца с информацией о штаммах бактерий.

    Возвращает:
    pd.DataFrame - исходный DataFrame с добавленными фичами.
    """
    bacteria_features_list = []

    # Проходим по каждому штамму бактерий в датасете
    for bacteria_info in df[bacteria_column]:
        features = generate_bacteria_features(
            bacteria_info
        )  # Генерация фичей для каждого штамма
        bacteria_features_list.append(features)

    # Преобразуем список фичей в DataFrame
    features_df = pd.DataFrame(bacteria_features_list)

    # Объединяем с исходным DataFrame
    df_with_features = pd.concat([df, features_df], axis=1)

    return df_with_features
