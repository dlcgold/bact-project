## lista di disease con drug
## lista di disease senza drug


# entrez


#1. kegg genome con codice disease
#esempio https://www.genome.jp/dbget-bin/get_linkdb?-t+genome+ds:H00113


#2. kegg Genome entry
    #prendo i bioProject (ho un submission e i biosample)

    # prendo assembly, ho biosample e bioproject, e genome
#3. prendo l'assembly con le statistics (non esiste su ENTREZ, devo passare da Genome)
    # BIOsample con Geotagging (per i bacterical)  (posso arrivarci direttamente da Entrez)
    # GENOME
        # gc content
        # protein count
        # protein list 
            # protein sample con geotagging (ridondante con biosample)

        

# biosamples

from Bio import Entrez
from Bio.KEGG import REST
from numpy import result_type
# from drugbank_sql_lite import *
# from bact_classes import * 
import requests as req
from bs4 import BeautifulSoup as bfs
from genomejp import *
import pickle
from utils import *

Entrez.email = "m.sgro2@campus.unimib.it"

def main():


    print("getting bacterial infections list")
    type_infection = "Bacterial infections"
    bactetial_infections = REST.kegg_get("br:br08401").read()
    H_list = []
    check = True
    bact_bool = False

    # get all H ID
    for line in bactetial_infections.splitlines():
        if not check:
            break
        if line == f"A{type_infection}":
            bact_bool = True
            continue
        if bact_bool:
            if line[0] == 'C':
                H_list.append(line.split()[1])
        if bact_bool and line[0] == 'A':
            check = False

    # fetch KEGG data and divide into the two subgroups and save in files
    print("Fetch data from KEGG")
    drug_list = []
    no_drug_list = []
    if not os.path.exists('kegg_get') or not os.listdir("kegg_get"):
        if not os.path.exists('kegg_get'):
            os.makedirs('kegg_get')
        if not os.path.exists('kegg_get_drug'):
            os.makedirs('kegg_get_drug')
        for h_id in H_list:
            bact_tmp = REST.kegg_get(h_id).read()
            if 'DRUG' not in bact_tmp:
                with open(f"kegg_get/{h_id}.txt", "w") as f:
                    f.write(bact_tmp)
                drug_list.append(h_id)
            else:
                with open(f"kegg_get_drug/{h_id}.txt", "w") as f:
                    f.write(bact_tmp)
                no_drug_list.append(h_id)

    # folder where plots will be stored
    if not os.path.exists('plot_print'):
        os.makedirs('plot_print')

    # fill lists of bacterials with all the data parsed from KEGG
    print("parse KEGG data")

    folder = 'tmp_bact'      
    kegg_folder = 'kegg_get'
    bacts = []

    if not os.path.exists(folder):
        os.makedirs(folder)
    if len(os.listdir(folder)) == 0:
        for filename in os.listdir(kegg_folder):
            with open(kegg_folder+"/" + filename, "r") as f:
                tmp_bact = parse(str(f.read()))
                bacts.append(tmp_bact)
                with open(folder + "/" + filename, "wb") as f:
                    pickle.dump(tmp_bact, f)
    else:
        for filename in os.listdir(folder):
            with open(folder +"/" + filename, "rb") as f:
                bacts.append(pickle.load(f))
        diff = list_diff(os.listdir(kegg_folder), os.listdir(folder))
        for filename in diff:
            with open(kegg_folder + "/" + filename, "r") as f:
                tmp_bact = parse(str(f.read()))
                bacts.append(tmp_bact)
                with open(folder+ "/" + filename, "wb") as f:
                    pickle.dump(tmp_bact, f)

    print(f"{len(bacts)} bacterial infections found w/out drugs")

    
    # Filling list with disease with drugs
    print("parse KEGG data with drugs")

    folder = 'tmp_bact_drug'      
    kegg_folder = 'kegg_get_drug'  

    bacts_drug = []
    if not os.path.exists(folder):
        os.makedirs(folder)
    if len(os.listdir(folder)) == 0:
        for filename in os.listdir(kegg_folder):
            with open(kegg_folder+"/" + filename, "r") as f:
                tmp_bact = parse(str(f.read()))
                bacts_drug.append(tmp_bact)
                with open(folder + "/" + filename, "wb") as f:
                    pickle.dump(tmp_bact, f)
    else:
        for filename in os.listdir(folder):
            with open(folder +"/" + filename, "rb") as f:
                bacts_drug.append(pickle.load(f))
        diff = list_diff(os.listdir(kegg_folder), os.listdir(folder))
        for filename in diff:
            with open(kegg_folder + "/" + filename, "r") as f:
                tmp_bact = parse(str(f.read()))
                bacts_drug.append(tmp_bact)
                with open(folder+ "/" + filename, "wb") as f:
                    pickle.dump(tmp_bact, f)

    
    print(f"{len(bacts_drug)} bacterial infections found w/ drugs")



    # genomes = get_genome_id_from_disease_id("H01451")
    # assemblies = []
    # biosamples = []
    # for ids in genomes:
    #     ass = get_assembly_for_id(ids)
    #     bios = get_biosample_for_id(ids)
    #     assemblies.append(ass)
    #     biosamples.append(bios)

    # # KEGG.GENOME from KEGG
    # print(genomes)
    # # ASSEMBLIES form KEGG.GENOME
    # print(assemblies)
    # # BIOSAMPLES form KEGG.GENOME
    # print(biosamples)

    # # assembly_entry = Entrez.efetch(db="all", id=assemblies[0])
    # # print(assembly_entry.read())
    # # assembly_entry = Entrez.esearch(term=assemblies[0], db="pubmed", retmax=1)
    # # print(assembly_entry.read())
    # for assembly in assemblies:
    #     query = f"https://www.ncbi.nlm.nih.gov/assembly/{assembly}"
    #     res = req.get(query)

    #     # print(res.text)
    #     soup = bfs(res.content, 'html.parser')
    #     soup_find = soup.find(class_="assembly_summary_new margin_t0")
    #     titles = []
    #     descriptions = []
    #     for title in soup_find.find_all("dt"):
    #         titles.append(title.text)
    #     for description in soup_find.find_all("dd"):
    #         descriptions.append(description.text)
    #     i = 0
    #     bioid = ""
    #     for _ in titles:
    #         if titles[i] == "BioSample: ":
    #             bioid = descriptions[i]
    #             print(titles[i], descriptions[i])
    #         elif titles[i] == "Submitter: ":
    #             print(titles[i], descriptions[i])
    #         elif titles[i] == "Organism name: ":
    #             print(titles[i], descriptions[i])
    #         elif titles[i] == "Date: ":
    #             print(titles[i], descriptions[i])
    #         i += 1


    #     bio3 = Entrez.efetch(id=bioid, db="biosample", rettmode='text')
    #     bio3_text = bio3.read()
    #     bio3_str = str(bio3_text)
    #     position_geotag = bio3_str.find('geo_loc_name')
    #     position_geotag_begin = bio3_str.find('>', position_geotag)
    #     position_geotag_end = bio3_str.find('</Attribute>', position_geotag)
    #     if (position_geotag != -1):
    #         print("Geotag:  ", bio3_str[position_geotag_begin+1:position_geotag_end]) 
            


        # for i in soup.find_all(class_="assembly_summary_new margin_t0"):
        #     print(i.text)
        # try: 
        #     list = soup.select('a[href*="/entry/gn:"]')
        #     genomes = [l['href'][-6:] for l in links]
        # except:
        #     print("fail")
        
    # bio_proj = Entrez.esearch(term="GCA_000092985.1", db="genome", rettype="gb", rettmode='text')
    # genome_id = Entrez.read(bio_proj)["IdList"][0]
    # print(genome_id)
    # genome_ret = Entrez.efetch(db="genome", id=genome_id, retmode="xml")
    # # print(Entrez.read(genome_ret))

    # bio2 = Entrez.efetch(id=biosamples[0], db="bioproject", rettmode='text')
    # # print(Entrez.read(bio2))
    # # print(bio2.read())
    # bio3 = Entrez.efetch(id="SAMN02603853", db="biosample", rettmode='text')
    # # print(bio3.read())

if __name__ == "__main__":
    main()