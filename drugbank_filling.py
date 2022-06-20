import os
from utils import *
from drugbank_sql_lite import * 
import plotly as plt
import numpy as np
import pandas as pd
from kegg_helper import *

### 26 FILLATI ISOLATI da drugbank
### filling extra1
### extra bruteforce

### filling extra2
### ogni aggiunta salvo un diz pathogeno-drug
### e ogni volta che fillo, controllo in quali pathogeni ho la drug, ed eventualmente aggiungo la tripla,
### faccio il dataframe,

### Controllo dei merge

###
# Lista di batteri con le drug annotate per patogeno

# Retrivo la lista delle drug su db dato il patogeno,
# Tolgo quelle già presenti nel bact, per id
# Aggiungo quelle di drugbank inserendo il nome (l'ID della drug è riferito a kegg)
# Serializzo per tutti i bact come le altre iterazioni


##import bacts
count = 0
count_drug = 0
bacts = []
for filename in os.listdir("genomejp_patho"):
    with open("genomejp_patho/" + filename, "rb") as f:
        bacts.append(pickle.load(f))

type = "db_extra"
print("Extra Filling from Drugbank")
bacts_db = []
pathos = []
pathos_extra = []
bact_count = 0
if not os.path.exists(type) or not os.listdir(type):
    if not os.path.exists(type):
        os.makedirs(type)
    for tmp_bact in bacts:
        bact_count += 1
        print("Current bact", bact_count, "of", len(bacts), "---------------------------------")
        if len(tmp_bact.pathogens) > 0:
            tmp_extra_drugs = []
            tmp_extra_db = []
            drug_patho = {}
            drug_patho_extra = {}

            for patho_tmp in tmp_bact.pathogens:
                print(patho_tmp)
                db_drugs = get_drugs_for_all(patho_tmp)
                print(db_drugs)
                tmp_drugs = []
                extra_drugs = []
                for tmp in db_drugs:
                    ## Avoid 400 bad requests
                    try:
                        id_tmp = conv_db_id_to_kegg_id(tmp)
                    except:
                        id_tmp = None
                    ##
                    print(tmp, "converted_to", id_tmp)
                    if id_tmp != "None":
                        tmp_drugs.append(id_tmp)
                    else:
                        extra_drugs.append(tmp)
                drug_patho[patho_tmp] = tmp_drugs
                drug_patho_extra[patho_tmp] = extra_drugs

                tmp_extra_drugs += tmp_drugs
                tmp_extra_db += extra_drugs
            if len(tmp_extra_drugs) > 0:
                tmp_drugs = []
                for elem in tmp_bact.drugs:
                    tmp_drugs.append(elem.id_drug)
                diffs = list_diff(tmp_extra_drugs, tmp_drugs)
                print("adding extra drugs from drugbank by pathogens")
                if len(diffs) > 0:
                    for elem in diffs:
                        print(elem)
                        if elem not in tmp_drugs:
                            if elem != None:
                                tmp_drug = get_drug_kegg(elem)
                                tmp_drug.origin = "DrugBank pathogen"
                                tmp_bact.drugs.append(tmp_drug)
   
                                for k, v in drug_patho.items():
                                    if elem in v:
                                        pathos.append((k, "", tmp_bact.sub))

            if len(tmp_extra_db) > 0 :
                count +=1
                count_drug += len(tmp_extra_db)
                for elem in tmp_extra_db:
                    print("Drugbank extra FILLING")
                    tmp_name = get_name_drug(elem)
                    print(tmp_name)
                    tmp_bact.drugs.append(Drug(tmp_name, ""))
                    for k, v in drug_patho_extra.items():
                        if elem in v:
                            pathos_extra.append((k, "", tmp_bact.sub))

        with open(f"{type}/{tmp_bact.id_bact}.txt", "wb") as f:
            pickle.dump(tmp_bact, f)
        bacts_db.append(tmp_bact)
    with open(f"{type}/extra_1.txt", "wb") as f:
        pickle.dump(pathos, f)
    with open(f"{type}/extra_pathos.txt", "wb") as f:
        pickle.dump(pathos_extra, f)


else:
    for filename in os.listdir(type):
        with open(f"{type}/{filename}", "rb") as f:
            bacts_db.append(pickle.load(f))
    # with open(f"db_patho/extra.txt", "rb") as f:
    #     patho_extra=pickle.load(f)
bacts_db_without_drugs = []
bacts_db_with_drugs = []
for bact_tmp in bacts_db:
    if len(bact_tmp.drugs) == 0:
        bacts_db_without_drugs.append(bact_tmp)
    else:
        bacts_db_with_drugs.append(bact_tmp)
print(
    f"{len(bacts_db_without_drugs)} extra FIlling from DrugBank"
    f"by pathogen")
print("TOTAL count - filled Drugs", len(bacts_db_with_drugs))

print("filled drugs", count_drug)
print("filled bact entries with at least 1 drug:", count)


# with open(f"extra_db.txt", "wb") as f:
#     pickle.dump(bacts_db_with_drugs, f)












# bacts = []
# bacts_db = []
# for filename in os.listdir("db_patho"):
#     with open("db_patho/" + filename, "rb") as f:
#         bacts.append(pickle.load(f))

# print("Filling from drugbank by pathogen - Extra Fill")
# if not os.path.exists('db_patho_extra'):
#     os.makedirs('db_patho_extra')

# for tmp_bact in bacts:   
#     if len(tmp_bact.pathogens) > 0:
#         tmp_extra_drugs = []
#         for patho_tmp in tmp_bact.pathogens:
#             print(patho_tmp)
#             db_drugs = get_drugs_for_all(patho_tmp)
#             print(db_drugs)
#             tmp_drugs = []
#             for tmp in db_drugs:
#                 try:
#                     tmp_name = get_name_drug(tmp)
#                 except:
#                     tmp_name = ""
#                 #
#                 print(tmp_name)
#                 if tmp_name != "":
#                     tmp_drugs.append(tmp_name)
       
#             tmp_extra_drugs += name_drugs
#             tmp_bact.drugs.append(tmp_extra_drugs)
#             # bacts_db = bacts_db.append(tmp_extra_drugs)
#     with open(f"db_patho_extra/{tmp_bact.id_bact}.txt", "wb") as f:
#         pickle.dump(tmp_bact, f)
#     bacts_db.append(tmp_bact)
# print(bacts)
# bacts_db_without_drugs = []
# bacts_db_with_drugs = []
# for bact_tmp in bacts_db:
#     if len(bact_tmp.drugs) == 0:
#         bacts_db_without_drugs.append(bact_tmp)
#     else:
#         bacts_db_with_drugs.append(bact_tmp)
# print(
#     f"{len(bacts_db_without_drugs)} found w/out drugs after drugbank filling "
#     f"by pathogen")

#             # if len(tmp_extra_drugs) > 0:
#             #     tmp_drugs = []
#             #     for elem in tmp_bact.drugs:
#             #         tmp_drugs.append(elem.id_drug)
#             #     diffs = list_diff(tmp_extra_drugs, tmp_drugs)
#             #     print("adding extra drugs from drugbank by pathogens")
#             #     if len(diffs) > 0:
#             #         for elem in diffs:
#             #             print(elem)
#             #             if elem not in tmp_drugs:
#             #                 if elem != None:
#             #                     tmp_drug = get_drug_kegg(elem)
#             #                     tmp_drug.origin = "DrugBank pathogen"
#             #                     tmp_bact.drugs.append(tmp_drug)
# # with open(f"db_patho/{tmp_bact.id_bact}.txt", "wb") as f:
# #     pickle.dump(tmp_bact, f)
# # bacts_db.append(tmp_bact)

# # bacts_db_without_drugs = []
# # bacts_db_with_drugs = []
# # for bact_tmp in bacts_db:
# #     if len(bact_tmp.drugs) == 0:
# #         bacts_db_without_drugs.append(bact_tmp)
# #     else:
# #         bacts_db_with_drugs.append(bact_tmp)
# # print(
# #     f"{len(bacts_db_without_drugs)}  found w/out drugs after drugbank filling "
# #     f"by pathogen")


    

# # pathos = []
# # for tmp_bact in bacts:
# #     print(tmp_bact.id_bact)
# #     if len(tmp_bact.pathogens) > 0:
# #         tmp_extra_drugs = []
# #         drug_patho = {}
# #         for patho_tmp in tmp_bact.pathogens:
# #             db_drugs = get_drugs_for_all(patho_tmp)
# #             tmp_drugs = []
# #             for tmp in db_drugs:
# #                 ## Avoid 400 bad requests
# #                 try:
# #                     id_tmp = conv_db_id_to_kegg_id(tmp)
# #                 except:
# #                     id_tmp = None
# #                 ##
# #                 print(id_tmp)
# #                 if id_tmp != "None":
# #                     tmp_drugs.append(id_tmp)
# #             drug_patho[patho_tmp] = tmp_drugs
# #             tmp_extra_drugs += tmp_drugs
# #         if len(tmp_extra_drugs) > 0:
# #             tmp_drugs = []
# #             for elem in tmp_bact.drugs:
# #                 tmp_drugs.append(elem.id_drug)
# #             diffs = list_diff(tmp_extra_drugs, tmp_drugs)
# #             print("adding extra drugs from drugbank by pathogens")
# #             if len(diffs) > 0:
# #                 for elem in diffs:
# #                     if elem not in tmp_drugs:
# #                         if elem is not None:tmp_extra_drugs
# #                                     pathos.append((k, "", tmp_bact.sub))
# #                             tmp_drug = get_drug_kegg(elem)
# #                             tmp_drug.origin = "DrugBank pathogen"
# #                             tmp_bact.drugs.append(tmp_drug)
# #             print(pathos)

