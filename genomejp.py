import pandas as pd
import requests as req
from bs4 import BeautifulSoup as bfs

import os
import pickle
from Bio.KEGG import REST

def get_drugs_for_id(id):
    url = f"https://www.genome.jp/entry/{id}"

    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0'
    }

    response = req.get(url, headers)
    soup = bfs(response.content, 'html.parser')
    links = soup.select('a[href*="/entry/D"]')
    drugs = [l['href'][-6:] for l in links]
    
    return drugs

def get_drugs_for_disease_name(disease):
    query = disease.replace(" ", "+")
    url = f"https://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&dbkey=drug_all&keywords={query}"

    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0'
    }

    response = req.get(url, headers)
    soup = bfs(response.content, 'html.parser')
    links = soup.select('a[href*="/entry/D"]')
    drugs = [l['href'][-6:] for l in links]
    
    return drugs

def get_drugs_for_pathway_kegg(pathway):
    if not os.path.exists('pathways_obj'):
        return get_drugs_for_id(pathway)
    else:
        with open(f"pathways_obj/{pathway}.ser", "rb") as f:
            p = pickle.load(f)
            return p['drugs']

#disease, malaria
#id = "H00361"
#pathway
id = "hsa04340"
tmp = get_drugs_for_pathway_kegg(id)
#tmp = get_drugs_for_disease_name("japanese encephalitis")

print(tmp)