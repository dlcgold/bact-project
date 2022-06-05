import sqlite3
import re
import pandas as pd
import requests as req
from bs4 import BeautifulSoup as bfs


def connect_db():
    return sqlite3.connect('drugbank_all_full_database.xml/drugbank_5.1.9.db')


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


def db_uniprot_conv(drug_id):
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


def get_drugs_for_target(target_id):
    conn = connect_db()
    drugs = []
    target_wild = f"%{target_id}%"
    sql = "SELECT `drugbank-id` FROM dbdf WHERE targets LIKE ?"
    cursor = conn.execute(sql, (target_wild,))
    for tmp in cursor:
        drugs.append(tmp[0])
    return drugs


def get_drugs_for_enzyme(enzyme_id):
    conn = connect_db()
    drugs = []
    enzyme_wild = f"%{enzyme_id}%"
    sql = "SELECT `drugbank-id` FROM dbdf WHERE enzymes LIKE ?"
    cursor = conn.execute(sql, (enzyme_wild,))
    for tmp in cursor:
        drugs.append(tmp[0])
    return drugs


def get_drugs_for_carriers_transporters(ct_id):
    conn = connect_db()
    drugs = []
    ct_wild = f"%{ct_id}%"
    sql = "SELECT `drugbank-id` FROM dbdf WHERE carriers LIKE ? OR transporters LIKE ?"
    cursor = conn.execute(sql, (ct_wild, ct_wild,))
    for tmp in cursor:
        drugs.append(tmp[0])
    return drugs


def get_drugs_for_all(target_id):
    conn = connect_db()
    drugs = []
    target_wild = f"%{target_id}%"
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


# tmp = db_uniprot_conv("DB00014")
# tmp = db_pubchem_conv("DB00014")
# tmp = get_name_drug("DB14738")
# tmp = get_id_drug("Zinc")
# tmp = get_targets_drug("DB14738")
# tmp = get_description_drug("DB14738")
# tmp = get_pharmacodyn_drug("DB14738")
# tmp = get_indication_drug("DB01175")
# tmp = get_pathways_drug("DB07718")
# tmp = get_enzymes_drug("DB09130")
# tmp = get_carriers_transporters_drug("DB09130")
# tmp = get_drugs_for_target("BE0005831")
# tmp = get_drugs_for_enzyme("BE0000267")
# tmp = get_drugs_for_carriers_transporters("BE0000530")
# tmp = get_drugs_for_all("BE0000530")
# tmp = get_target_name("BE0000530")
# tmp = get_drugs_inter_for_drug("DB15865")
# tmp = get_patents_drug("DB01175")

print(tmp)
