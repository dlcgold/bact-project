import sqlite3
import re
import pandas as pd
import requests as req
from bs4 import BeautifulSoup as bfs

import urllib


# TODO unire query per ottenere oggetti Target senza passare da scraping

class Target:
    def __init__(self, db_id, uniprot_id, uniprot_name):
        self.db_id = db_id
        self.uniprot_id = uniprot_id
        self.uniprot_name = uniprot_name

    def __repr__(self):
        return f"{self.uniprot_name}, ({self.db_id, self.uniprot_id})"


def connect_db():
    return sqlite3.connect('drugbank_5.1.9.db')


def db_uniprof_df():
    return pd.read_csv('drug_target_uniprot_links_5.1.2.csv')


def db_pubchem_df():
    return pd.read_csv('db_pubchem.tsv', sep="\t")


def get_target_name(target_id):
    # STA ROBA é LENTA NON SO SCRAPARE
    url = f"https://go.drugbank.com/bio_entities/{target_id}"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0'
    }
    response = req.get(url, headers)
    soup = bfs(response.content, 'html.parser')
    # for tmp in soup.find(class_="card-content px-md-4 px-sm-2 pb-md-4 pb-sm-2").find_all(class_="col-xl-10 col-md-9 col-sm-8"):
    #    print(tmp)
    return soup.title(text=True)[0].split("|")[0].strip()


def db_target_uniprot(drug_id):
    df = db_uniprof_df()
    uniprots = []
    df_tmp = df[df["DrugBank ID"] == drug_id]
    for idx, row in df_tmp[["UniProt ID", "UniProt Name"]].iterrows():
        uniprots.append((row["UniProt ID"], row["UniProt Name"]))
    return uniprots


def db_pubchem_conv(drug_id):
    df = db_pubchem_df()
    pubchem = ""
    df_tmp = df[df["drugbank_id"] == drug_id]
    for idx, row in df_tmp[["pubchem_cid"]].iterrows():
        pubchem = row["pubchem_cid"]
    return pubchem


def get_id_drug(drug_name):
    conn = connect_db()
    sql = "SELECT `drugbank-id` FROM dbdf WHERE name = ?"
    cursor = conn.execute(sql, (drug_name,))
    id = ""
    for row in cursor:
        id = row[0]
    return id


def get_name_drug(drug_id):
    conn = connect_db()
    sql = "SELECT name FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    name = ""
    for row in cursor:
        name = row[0]
    return name


def get_description_drug(drug_id):
    conn = connect_db()
    sql = "SELECT description FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    des = ""
    for row in cursor:
        des = row[0]
    return des


def get_indication_drug(drug_id):
    conn = connect_db()
    sql = "SELECT indication FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    ind = ""
    for row in cursor:
        ind = row[0]
    return ind


def get_pharmacodyn_drug(drug_id):
    conn = connect_db()
    sql = "SELECT pharmacodynamics FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    phd = ""
    for row in cursor:
        phd = row[0]
    return phd


def get_groups_drug(drug_id):
    conn = connect_db()
    sql = "SELECT groups FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    group = ""
    for row in cursor:
        group = row[0]
    return group


def get_targets_drug(drug_id):
    conn = connect_db()
    targets = []
    sql = "SELECT targets FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    for tmp in cursor:
        # print(tmp[0])
        targets_id = re.findall("BE[0-9]{7}", tmp[0])
        targets += targets_id
    return targets


def get_targets_doi_drug(drug_id):
    conn = connect_db()
    targets = []
    sql = "SELECT targets FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    for tmp in cursor:
        # print(tmp[0])
        targets_id = re.findall("10[.][0-9]{4,9}[/][-._;()\/:A-Z0-9a-z]+[.]", tmp[0])
        targets += targets_id
    final_targets = []
    for elem in targets:
        if elem[-1] == '.':
            final_targets.append(elem.strip()[:len(elem) - 1])
        else:
            final_targets.append(elem)
    return final_targets


def get_targets_name_drug(drug_id):
    conn = connect_db()
    targets = []
    sql = "SELECT targets FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    for tmp in cursor:
        # print(tmp[0])
        targets_id = re.findall("BE[0-9]{7}[A-Z][a-z0-9\ ]*", tmp[0])
        targets += targets_id
    final_targets = []
    for elem in targets:
        final_targets.append(elem.strip()[9:len(elem)])
    return final_targets


def get_targets_full_drug(drug_id):
    # TODO vedere se si riesce a prendere l'id da uniprot
    conn = connect_db()
    uniprot = db_target_uniprot(drug_id)
    targets = []
    sql = "SELECT targets FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    count = 0
    for tmp in cursor:
        targets_id = re.findall("BE[0-9]{7}", tmp[0])
        for tmp in targets_id:
            if count < len(uniprot):
                targets.append(Target(tmp, uniprot[count][0], uniprot[count][1]))
            else:
                name = get_target_name(tmp)
                if name is not None:
                    targets.append(Target(tmp, "NA", name))
                else:
                    targets.append(Target(tmp, "NA", "NA"))
            count += 1
    return targets


def get_pathways_drug(drug_id):
    conn = connect_db()
    pathways = []
    sql = "SELECT pathways FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    for tmp in cursor:
        # print(tmp[0])
        targets_id = re.findall("SMP[0-9]{7}", tmp[0])
        pathways += targets_id
    return pathways
    

def get_pathways_name_drug(drug_id):
    words_to_remove = ["diseaseD","Disease", "Type","Action","disease","type","action","Metabolic","metabolic", "Pathway","signaling","Signaling","Physiological","physiological","drug_action","drug_metabolism"]
    
    conn = connect_db()
    pathways = []
    sql = "SELECT pathways FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    for tmp in cursor:
        # print(tmp[0])
        targets_id = re.findall("SMP[0-9]{7}[A-Z][a-z0-9]*[ and ]*[A-Z][a-z0-9]*", tmp[0])
        pathways += targets_id
        
    final_targets = []
    for elem in pathways:
        raw = elem.strip()[10:len(elem)]
        for w in words_to_remove:
            raw = raw.replace(w,"")
        clean = raw.strip()
        final_targets.append(clean)
    
    return set(final_targets)

def get_enzymes_drug(drug_id):
    conn = connect_db()
    enzymes = []
    sql = "SELECT enzymes FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    for tmp in cursor:
        # print(tmp[0])
        enzymes_id = re.findall("BE[0-9]{7}", tmp[0])
        enzymes += enzymes_id
    return enzymes


def get_carriers_transporters_drug(drug_id):
    conn = connect_db()
    carriers = []
    trasporters = []
    sql = "SELECT carriers, transporters FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    for tmp in cursor:
        carriers_id = re.findall("BE[0-9]{7}", tmp[0])
        carriers += carriers_id
        trasporters_id = re.findall("BE[0-9]{7}", tmp[1])
        trasporters += trasporters_id
    return [carriers, trasporters]


def get_drugs_for_target(target):
    conn = connect_db()
    drugs = []
    target_wild = f"%{target}%"
    sql = "SELECT `drugbank-id` FROM dbdf WHERE targets LIKE ?"
    cursor = conn.execute(sql, (target_wild,))
    for tmp in cursor:
        drugs.append(tmp[0])
    return drugs


def get_drugs_for_enzyme(enzyme):
    conn = connect_db()
    drugs = []
    enzyme_wild = f"%{enzyme}%"
    sql = "SELECT `drugbank-id` FROM dbdf WHERE enzymes LIKE ?"
    cursor = conn.execute(sql, (enzyme_wild,))
    for tmp in cursor:
        drugs.append(tmp[0])
    return drugs

def get_drugs_for_pathway(pathway):
    conn = connect_db()
    drugs = []
    enzyme_wild = f"%{pathway}%"
    sql = "SELECT `drugbank-id` FROM dbdf WHERE pathways LIKE ?"
    cursor = conn.execute(sql, (enzyme_wild,))
    for tmp in cursor:
        drugs.append(tmp[0])
    return drugs

def get_drugs_for_carriers_transporters(ct):
    conn = connect_db()
    drugs = []
    ct_wild = f"%{ct}%"
    sql = "SELECT `drugbank-id` FROM dbdf WHERE carriers LIKE ? OR transporters LIKE ?"
    cursor = conn.execute(sql, (ct_wild, ct_wild,))
    for tmp in cursor:
        drugs.append(tmp[0])
    return drugs


def get_drugs_for_all(target):
    conn = connect_db()
    drugs = []
    target_wild = f"%{target}%"
    sql = "SELECT `drugbank-id` FROM dbdf " \
          "WHERE targets LIKE ? OR enzymes LIKE ? " \
          "OR carriers LIKE ? OR transporters LIKE ?"
    cursor = conn.execute(sql, (target_wild, target_wild, target_wild, target_wild))
    for tmp in cursor:
        drugs.append(tmp[0])
    return drugs


def get_drugs_inter_for_drug(drug_id):
    conn = connect_db()
    drugs = []
    sql = "SELECT `drug-interactions` FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    for tmp in cursor:
        # print(tmp[0])
        drugs_id = re.findall("DB[0-9]{5}", tmp[0])
        drugs += drugs_id
    return drugs


def get_patents_drug(drug_id):
    conn = connect_db()
    patents = []
    sql = "SELECT `patents` FROM dbdf WHERE `drugbank-id` = ?"
    cursor = conn.execute(sql, (drug_id,))
    for tmp in cursor:
        # print(tmp[0])
        patents_id = re.findall("[0-9]{7}", tmp[0])
        patents += patents_id
        patents_id = re.findall("RE[0-9]{5}", tmp[0])
        patents += patents_id
    return patents

# tmp = db_target_uniprot("DB00004")
# tmp = db_pubchem_conv("DB00014")
# tmp = get_name_drug("DB14738")
# tmp = get_id_drug("Zinc")
# tmp = get_targets_drug("DB00014")
# tmp = get_targets_full_drug("DB14738")
# tmp = get_description_drug("DB14738")
# tmp = get_pharmacodyn_drug("DB14738")
# tmp = get_groups_drug("DB14738")
# tmp = get_indication_drug("DB01175")
# tmp = get_pathways_drug("DB07718")
# tmp = get_pathways_name_drug("DB01175")
# tmp = get_enzymes_drug("DB09130")
# tmp = get_carriers_transporters_drug("DB09130")
# tmp = get_drugs_for_target("BE0005831")
# tmp = get_drugs_for_enzyme("BE0000267")
# tmp = get_drugs_for_carriers_transporters("BE0000530")
# tmp = get_drugs_for_all("BE0000530")
# tmp = get_target_name("BE0000530")
# tmp = get_drugs_inter_for_drug("DB15865")
# tmp = get_patents_drug("DB01175")
# tmp = get_drugs_for_all("BE0005831")
# tmp = get_targets_doi_drug("DB00002")
# tmp = get_drugs_for_pathway("SMP0000006")
tmp = get_drugs_for_pathway("hsa05144")

print(tmp)
