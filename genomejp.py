from cmath import e
import os
import pickle

import requests as req
from bs4 import BeautifulSoup as bfs

def get_genome_id_from_disease_id(disease_id):
    query = f"ds:{disease_id}"  
    url = f"https://www.genome.jp/dbget-bin/get_linkdb?-t+genome+{query}"
    response = req.get(url)
    soup = bfs(response.content, 'html.parser')
    try: 
        links = soup.select('a[href*="/entry/gn:"]')
        genomes = [l['href'][-6:] for l in links]
    except:
        genomes = []
    return genomes

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


def get_assembly_biosample(id):
    url = f"https://www.genome.jp/entry/{id}"

    response = req.get(url)
    soup = bfs(response.content, 'html.parser')
    try:
        links = soup.select('a[href*="/assembly"]')
        assembly = [l['href'][38:] for l in links]
    except:
        assembly = ''
    try:
        links = soup.select('a[href*="/bioproject"]')
        bios = [l['href'][40:] for l in links]
    except:
        bios = ''
    try:
        return(assembly[0], bios[0])
    except:
        return("","")

def get_assembly_for_id(id):
    url = f"https://www.genome.jp/entry/{id}"

    response = req.get(url)
    soup = bfs(response.content, 'html.parser')
    try:
        links = soup.select('a[href*="/assembly"]')
        assembly = [l['href'][38:] for l in links]
    except:
        assembly = ''
    return assembly[0]

def get_biosample_for_id(id):
    url = f"https://www.genome.jp/entry/{id}"

    response = req.get(url)
    soup = bfs(response.content, 'html.parser')
    try:
        links = soup.select('a[href*="/bioproject"]')
        bios = [l['href'][40:] for l in links]
    except:
        bios = ''
    return bios[0]

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
    if not os.path.exists(f"data/pathways_obj/{pathway}.ser"):
        print("not serialized")
        return get_drugs_for_id(pathway)
    else:
        with open(f"data/pathways_obj/{pathway}.ser", "rb") as f:
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
