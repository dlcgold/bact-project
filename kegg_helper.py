from Bio import Entrez
from Bio.KEGG import REST
from drugbank_sql_lite import *

Entrez.email = "d.cozzi@campus.unimib.it"


def get_pathway_target_from_kegg(drug_id):
    data = REST.kegg_get(drug_id).read()
    pathways = []
    path_bool = False
    for line in data.splitlines():
        spl = line.strip().split()
        if spl[0] == "PATHWAY":
            name = "".join(spl[i] + " " for i in range(2, len(spl)))
            pathways.append((spl[1].split("(")[0], name))
            path_bool = True
        elif spl[0].upper() != spl[0] and path_bool:
            name = "".join(spl[i] + " " for i in range(1, len(spl)))
            pathways.append((spl[0].split("(")[0], name))
            path_bool = True
        elif spl[0].upper() == spl[0] and path_bool:
            break
    return pathways


def conv_name_to_kegg_id(drug_name):
    drug_name = drug_name.replace(" ", "_")
    data = REST.kegg_find("Drug", drug_name).read()
    if data.splitlines()[0] != "":
        return data.splitlines()[0].split()[0].split(":")[1]
    else:
        return "None"


def conv_id_to_kegg_name(drug_id):
    drug_id = drug_id.replace(" ", "_")
    data = REST.kegg_get(drug_id).read()
    names = []
    name_bool = False
    for line in data.splitlines():
        spl = line.strip().split()
        if spl[0] == "NAME":
            names.append((spl[1].split("(")[0]))
            name_bool = True
        elif spl[0].upper() != spl[0] and name_bool:
            names.append((spl[0].split("(")[0]))
            name_bool = True
        elif spl[0].upper() == spl[0] and name_bool:
            break
    return list(set(names))


def conv_db_id_to_kegg_id(db_id):
    name = get_name_drug(db_id)
    return conv_name_to_kegg_id(name)


def conv_kegg_id_to_db_id(kegg_id):
    kegg_names = conv_id_to_kegg_name(kegg_id)
    names = []
    for tmp_name in kegg_names:
        if get_id_drug(tmp_name) != "":
            names.append(get_id_drug(tmp_name))
    return list(set(names))


# tmp = get_pathway_target_from_kegg("D11713")
# name = get_name_drug("DB00135")
# print(name)
# tmp = conv_name_to_kegg_id(name)
# tmp = conv_db_id_to_kegg_id("DB00135")
# tmp = conv_id_to_kegg_name("D11713")
tmp = conv_kegg_id_to_db_id("D11713")
print(tmp)
