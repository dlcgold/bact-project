import os
import pickle

import requests as req
from bs4 import BeautifulSoup as bfs


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
    if not os.path.exists(f"pathways_obj/{pathway}.ser"):
        print("not serialized")
        return get_drugs_for_id(pathway)
    else:
        with open(f"pathways_obj/{pathway}.ser", "rb") as f:
            p = pickle.load(f)
            return p['drugs']

# kegg_drugs GET by id or name of patogen
def get_genomejp_drugs(query):
    query = query.replace(" ", "+")
    db = "drug"
    url = f"https://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=gn&dbkey={db}&keywords={query}"

    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0'
    }

    response = req.get(url, headers)
    soup = bfs(response.content, 'html.parser')
    links = []
    for link in soup.find_all('a'):
        links.append(link.get('href'))
    drugs = []
    for entry in links:
        if "entry" in entry:
            tmp = entry.split('/')[-1]
            if tmp[0] == "D":
                drugs.append(tmp)
    return list(set(drugs))
    
# disease, malaria
# id = "H00361"
# pathway
# id = "hsa05216"
# tmp = get_drugs_for_pathway_kegg(id)
# tmp = get_drugs_for_disease_name("japanese encephalitis")

# print(tmp)
