from numpy import true_divide
from utils import *
import pandas as pd
from drugbank_sql_lite import *
from kegg_helper import *

# type = "db_extra"
# with open(f"{type}/extra.txt", "rb") as f:
#     bacts = pickle.load(f)
bacts = []
type = "db_extra"
for filename in os.listdir(type):
    with open(f"{type}/{filename}", "rb") as f:
        bacts.append(pickle.load(f))
ids = []
for elem in bacts:
    if elem.drugs != []:
        ids.append((elem.id_bact).replace("ds:", ""))
print(ids)
# ids hardcoded
ids = ['H00277', 'H00300', 'H00301', 'H00298', 'H00303', 'H01067', 'H00304',
        'H01317', 'H01069', 'H00307', 'H00313', 'H01441', 'H00309', 'H01074',
        'H01080', 'H00317', 'H01144', 'H01151', 'H01166', 'H00329', 'H00331',
        'H01401', 'H01405', 'H01426', 'H01423', 'H01406', 'H01443', 'H00335',
        'H01444', 'H00341', 'H00388', 'H02029', 'H01462', 'H01051', 'H01446',
        'H00355', 'H02076', 'H01442', 'H01455', 'H01408', 'H01422']
# data_patho_extra = {}
# for elem in patho_extra:
#     print(elem.id_bact)
#     print(elem.name)
#     key = f"{elem[0]}|{elem[2]}"
#     if key not in data_patho_extra.keys():
#         data_patho_extra[key] = 1
#     else:
#         data_patho_extra[key] += 1
# df_patho_extra = pd.DataFrame(columns=["pathogen", "subcategory", "frequency"])

# for key, val in data_patho_extra.items():
#     tmp = key.split("|")
#     df_patho_extra.loc[len(df_patho_extra.index)] = [tmp[0],
#                                                      tmp[1],
#                                                      val]
# df_patho_extra = df_patho_extra.sort_values(by="frequency", ascending=False).reset_index(drop=True)
# print(df_patho_extra)

# ids = ['H01408', 'H01422', 'H00277', 'H01462', 'H00309', 'H00388', 'H01051',
#        'H01423', 'H01446', 'H00331', 'H01401', 'H01426', 'H02076', 'H01067',
#        'H02029', 'H01442', 'H00335', 'H00300', 'H01074', 'H00304', 'H00301',
#        'H01441', 'H01455', 'H00329', 'H00313', 'H01444', 'H01405']



# for filename in os.listdir("kegg_get"):
#     if filename.replace(".txt", "") in ids:
#         with open("kegg_get/" + filename, "r") as f:
#             bacts.append(parse(str(f.read())))
# pathos = []
# for tmp_bact in bacts:
#     print(tmp_bact.id_bact)
#     if len(tmp_bact.pathogens) > 0:
#         tmp_extra_drugs = []
#         drug_patho = {}
#         for patho_tmp in tmp_bact.pathogens:
#             db_drugs = get_drugs_for_all(patho_tmp)
#             tmp_drugs = []
#             for tmp in db_drugs:
#                 ## Avoid 400 bad requests
#                 try:
#                     id_tmp = conv_db_id_to_kegg_id(tmp)
#                 except:
#                     id_tmp = None
#                 ##
#                 print(id_tmp)
#                 if id_tmp != "None":
#                     tmp_drugs.append(id_tmp)
#             drug_patho[patho_tmp] = tmp_drugs
#             tmp_extra_drugs += tmp_drugs
#         if len(tmp_extra_drugs) > 0:
#             tmp_drugs = []
#             for elem in tmp_bact.drugs:
#                 tmp_drugs.append(elem.id_drug)
#             diffs = list_diff(tmp_extra_drugs, tmp_drugs)
#             print("adding extra drugs from drugbank by pathogens")
#             if len(diffs) > 0:
#                 for elem in diffs:
#                     if elem not in tmp_drugs:
#                         if elem is not None:
#                             for k, v in drug_patho.items():
#                                 if elem in v:
#                                     pathos.append((k, "", tmp_bact.sub))
#                             tmp_drug = get_drug_kegg(elem)
#                             tmp_drug.origin = "DrugBank pathogen"
#                             tmp_bact.drugs.append(tmp_drug)
#             print(pathos)
# with open(f"{type}/extra.txt", "wb") as f:
#     pickle.dump(pathos, f)

pathos=[]
ok = False
for tmp_bact in bacts:
    if len(tmp_bact.pathogens) > 0:
        drug_patho = {}
        for patho_tmp in tmp_bact.pathogens:
            for elem in tmp_bact.drugs:
                if elem.id_drug != "":
                    ok = True
            if ok:
                drug_patho[patho_tmp] = tmp_bact.drugs
                ok = False
            for k, v in drug_patho.items():
                pathos.append((k, "", tmp_bact.sub))



# type = "db_patho"
with open("extra_1.txt", "rb") as f:
    patho_extra = pickle.load(f)
data_patho_extra = {}
for elem in patho_extra:
        key = f"{elem[0]}|{elem[2]}"
        if key not in data_patho_extra.keys():
            data_patho_extra[key] = 1
        else:
            data_patho_extra[key] += 1
df_patho_extra = pd.DataFrame(columns=["pathogen", "subcategory", "frequency"])

val_tot = 0
for key, val in data_patho_extra.items():
    tmp = key.split("|")
    df_patho_extra.loc[len(df_patho_extra.index)] = [tmp[0],
                                                     tmp[1],
                                                     val]
    val_tot += val

df_patho_extra = df_patho_extra.sort_values(by="frequency", ascending=False).reset_index(drop=True)
df_patho_extra.loc[len(df_patho_extra.index)] = ["total",
                                                 "",
                                                 val_tot]
print(df_patho_extra)


with open("extra_pathos.txt", "rb") as f:
    patho_extra = pickle.load(f)
data_patho_extra = {}
for elem in patho_extra:
        key = f"{elem[0]}|{elem[2]}"
        if key not in data_patho_extra.keys():
            data_patho_extra[key] = 1
        else:
            data_patho_extra[key] += 1
df_patho_extra = pd.DataFrame(columns=["pathogen", "subcategory", "frequency"])

val_tot = 0
for key, val in data_patho_extra.items():
    tmp = key.split("|")
    df_patho_extra.loc[len(df_patho_extra.index)] = [tmp[0],
                                                     tmp[1],
                                                     val]
    val_tot += val

df_patho_extra = df_patho_extra.sort_values(by="frequency", ascending=False).reset_index(drop=True)
df_patho_extra.loc[len(df_patho_extra.index)] = ["total",
                                                 "",
                                                 val_tot]
print(df_patho_extra)

