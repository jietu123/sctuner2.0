import pandas as pd

sc_path = "data/processed/real_brca_simS0_seed0/stage1_preprocess/exported/sc_expression_normalized.csv"
st_path = "data/processed/real_brca_simS0_seed0/stage1_preprocess/exported/st_expression_normalized.csv"

sc = pd.read_csv(sc_path, sep=None, engine="python")
sc_t = sc.set_index(sc.columns[0]).T
sc_t.to_csv(sc_path)
print("sc_expr transposed:", sc_t.shape)

st = pd.read_csv(st_path, sep=None, engine="python")
st_t = st.set_index(st.columns[0]).T
st_t.to_csv(st_path)
print("st_expr transposed:", st_t.shape)

