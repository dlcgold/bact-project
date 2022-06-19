import pickle

import pandas as pd

with open(f"db_patho/extra.txt", "rb") as f:
    patho_extra = pickle.load(f)

data_patho_extra = {}
for elem in patho_extra:
    key = f"{elem[0]}|{elem[2]}"
    if key not in data_patho_extra.keys():
        data_patho_extra[key] = 1
    else:
        data_patho_extra[key] += 1
df_patho_extra = pd.DataFrame(columns=["pathogen", "subcategory", "frequency"])

for key, val in data_patho_extra.items():
    tmp = key.split("|")
    df_patho_extra.loc[len(df_patho_extra.index)] = [tmp[0],
                             tmp[1],
                             val]
df_patho_extra = df_patho_extra.sort_values(by="frequency", ascending=False).reset_index(drop=True)
print(df_patho_extra)
