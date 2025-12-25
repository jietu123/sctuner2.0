import pandas as pd
from pathlib import Path

p = Path("data/processed/real_brca_simS0_seed0/stage1_preprocess/exported/sc_expression_normalized.csv")
sc = pd.read_csv(p, sep=None, engine="python", nrows=3)
print(sc.shape)
print(sc.columns[:5])
print(sc.head())

