import pickle
from sysconfig import get_path
import pandas as pd
import requests as req
from bs4 import BeautifulSoup as bfs
import os
import re
import urllib
from Bio.KEGG import REST

def all_pathways():
    url = f"https://rest.kegg.jp/link/hsa/pathway"

    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0'
    }

    response = req.get(url, headers)
    soup = bfs(response.content, 'html.parser')

    all = re.split("\t|\n", soup.text)
    del all[-1]
    all = list(set(all))
    print("len all ", len(all))
    p = [a[5:] for a in all if a[0] == 'p']
    h = [a for a in all if a[0] == 'h']

    return p,h,all
    
def get_pathway_page(pathway):
    url = f"https://www.genome.jp/entry/{pathway}"

    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0'
    }

    response = req.get(url, headers)
    soup = bfs(response.content, 'html.parser')
    name = soup.select('td[class*="td31 defd"]')[0].text
    links = soup.select('a[href*="/entry/D"]')
    drugs = [l['href'][-6:] for l in links]
    page = soup.get_text()

    return name,drugs,page

def main():
    p,a,all = all_pathways()
    if not os.path.exists('pathways_txt') and not os.path.exists('pathways_obj'):
        os.makedirs('pathways_txt')
        os.makedirs('pathways_obj')
        for pathway in p:
            name, drugs, page = get_pathway_page(pathway)
            obj = {"id" : pathway, "name" : name, "drugs" : drugs}
            with open(f"pathways_obj/{pathway}.ser", "wb") as f:
                pickle.dump(obj, f)
            with open(f"pathways_txt/{pathway}.txt", "w") as f:
                f.write(page)
    """else:
        """

if __name__ == "__main__":
    main()