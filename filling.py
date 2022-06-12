from ast import dump
from pkg_resources import BINARY_DIST
from kegg_helper import *
from scrape_pathways import *
from drugbank_sql_lite import *
from utils import *
from bact_classes import *
from genomejp import *
import matplotlib.pyplot as plt
import math
import numpy as np


Entrez.email = "d.cozzi@campus.unimib.it"
Entrez.email = 'm.sgro2@campus.unimib.it'
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

    # fill lists of bacterials with all the data parsed from KEGG
    print("parse KEGG data")
    bacts = []
    for filename in os.listdir("kegg_get"):
        with open("kegg_get/" + filename, "r") as f:
            tmp_bact = parse(str(f.read()))
            bacts.append(tmp_bact)
    print(f"{len(bacts)} bacterial infections found w/out drugs")

    print("parse KEGG data with drugs")
    bacts_drug = []
    for filename in os.listdir("kegg_get_drug"):
        with open("kegg_get_drug/" + filename, "r") as f:
            tmp_bact = parse(str(f.read()))
            bacts_drug.append(tmp_bact)
    print(f"{len(bacts_drug)} bacterial infections found w/ drugs")

    # bar plot by quantity of infection by subgroups 
    bar_data = {}
    for bact in bacts:
        if bact.sub not in bar_data.keys():
            bar_data[bact.sub] = 1
        else:
            bar_data[bact.sub] += 1
    
    bar_data_drugs = {}
    for bact in bacts_drug:
        if bact.sub not in bar_data_drugs.keys():
            bar_data_drugs[bact.sub] = 1
        else:
            bar_data_drugs[bact.sub] += 1

    labels = []
    for x in bar_data.keys():
        labels.append(x.replace("Infections caused by", ""))
    labels_drug = []
    for x in bar_data_drugs.keys():
        labels_drug.append(x.replace("Infections caused by", ""))

    delta = np.arange(len(labels))
    width = 0.15
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.bar(delta, bar_data.values(), color='red', width = width, label = "Withoud drug")
    ax.bar(delta + width, bar_data_drugs.values(), color='green', width = width, label = "With drug")
    plt.title("Bacterial infections with/without drugs by subgroup")
    plt.ylabel("Quantity of infections") 
    plt.xticks(delta + width / 2, labels, rotation=45)
    plt.rc('xtick', labelsize=3)
    y_int = range(math.floor(min(bar_data.values())), math.ceil(max(bar_data_drugs.values()))+1)
    plt.yticks(y_int)
    ax.legend()
    fname = './plot_print/comparative_drugs_bar_plot.png'
    # plt.show()
    fig.set_size_inches((16,12), forward=False)
    fig.savefig(fname, dpi=500)
    

    # Filling from genomejp by pathogen query
    print("Filling from genomejp by pathogen query")
    bact_names = "" # for geo parsing
    bacts_genomejp = []
    if not os.path.exists('genomejp') or not os.listdir("genomejp"):
        if not os.path.exists('genomejp'):
            os.makedirs('genomejp')
        for tmp_bact in bacts:
            bact_names += (' ' + tmp_bact.name.lower())
            tmp_extra_drugs = get_drugs_for_disease_name(tmp_bact.name)
            if len(tmp_extra_drugs) > 0:
                tmp_drugs = []
                for elem in tmp_bact.drugs:
                    tmp_drugs.append(elem.id_drug)
                diffs = list_diff(tmp_extra_drugs, tmp_drugs)
                print("adding extra drugs from genomejp")
                if len(diffs) > 0:
                    for elem in diffs:
                        if elem not in tmp_drugs:
                            tmp_bact.drugs.append(get_drug_kegg(elem))
            with open(f"genomejp/{tmp_bact.id_bact}.txt", "wb") as f:
                    pickle.dump(tmp_bact, f)
            bacts_genomejp.append(tmp_bact)
    else:
        for filename in os.listdir("genomejp"):
            with open("genomejp/" + filename, "rb") as f:
                bacts_genomejp.append(pickle.load(f))
    print(len(bacts))
    no_drug_after_genomejp = []
    for bact in bacts_genomejp:
        if len(bact.drugs) == 0:
            no_drug_after_genomejp.append(bact)
    print(f"{len(no_drug_after_genomejp)} bacterial infections found w/out drugs after genomejp")


    
    if len(no_drug_after_genomejp) > 0:
        # bar plot by quantity of infection by subgroups after genomejp
        bar_data_genomejp = {}
        for bact in bacts:
            if bact.sub not in bar_data_genomejp.keys():
                bar_data_genomejp[bact.sub] = 1
            else:
                bar_data_genomejp[bact.sub] += 1
        print(bar_data_genomejp)
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.bar(delta, bar_data.values(), color='red', width = width, label = "Withoud drug")
        ax.bar(delta, bar_data_genomejp.values(), color='orange', width = width, label = "Filled with genomejp")
        ax.bar(delta + width, bar_data_drugs.values(), color='green', width = width, label = "With drug")
        plt.title("Bacterial infections with/without drugs by subgroup")
        plt.ylabel("Quantity of infections") 
        plt.xticks(delta + width / 2, labels, rotation=45)
        plt.rc('xtick', labelsize=3)
        y_int = range(math.floor(min(bar_data.values())), math.ceil(max(bar_data_drugs.values()))+1)
        plt.yticks(y_int)
        ax.legend()
        fname = './plot_print/comparative_drugs_bar_plot_after_genomejp.png'
        plt.show()
        fig.set_size_inches((16,12), forward=False)
        fig.savefig(fname, dpi=500)

    print("Filling from genomejp by pathogen query")
    bact_names = "" # for geo parsing
    bacts_genomejp_pathway = []
    if not os.path.exists('genomejp') or not os.listdir("genomejp"):
        if not os.path.exists('genomejp'):
            os.makedirs('genomejp')
        for tmp_bact in bacts:
            bact_names += (' ' + tmp_bact.name.lower())
            tmp_extra_drugs = get_drugs_for_disease_name(tmp_bact.name)
            if len(tmp_extra_drugs) > 0:
                tmp_drugs = []
                for elem in tmp_bact.drugs:
                    tmp_drugs.append(elem.id_drug)
                diffs = list_diff(tmp_extra_drugs, tmp_drugs)
                print("adding extra drugs from genomejp")
                if len(diffs) > 0:
                    for elem in diffs:
                        if elem not in tmp_drugs:
                            tmp_bact.drugs.append(get_drug_kegg(elem))
            with open(f"genomejp/{tmp_bact.id_bact}.txt", "wb") as f:
                    pickle.dump(tmp_bact, f)
                bacts_genomejp_pathway.append(tmp_bact)
    else:
        for filename in os.listdir("genomejp"):
            with open("genomejp/" + filename, "rb") as f:
                bacts_genomejp_pathway.append(pickle.load(f))
    # print(len(bacts))
    no_drug_after_genomejp = []
    for bact in bacts_genomejp_pathway:
        if len(bact.drugs) == 0:
            no_drug_after_genomejp.append(bact)
    print(f"{len(no_drug_after_genomejp)} {type_infection} found w/out drugs after genomejp")



    print("Filling from drugbank by pathogen")
    bacts_db = []
    if not os.path.exists('db_patho') or not os.listdir("db_patho"):
        if not os.path.exists('db_patho'):
            os.makedirs('db_patho')
        for tmp_bact in bacts:
            if len(tmp_bact.pathogens) > 0:
                tmp_extra_drugs = []
                for patho_tmp in tmp_bact.pathogens:
                    db_drugs = get_drugs_for_all(patho_tmp)
                    print(db_drugs)
                    tmp_drugs = []
                    for tmp in db_drugs:
                        print(tmp)
                        ## Aggiunto per evitare 400 bad requests
                        try: id_tmp = conv_db_id_to_kegg_id(tmp)
                        except: id_tmp = None
                        ##
                        print(id_tmp)
                        if id_tmp != "None":
                            tmp_drugs.append(id_tmp)
                    print(tmp_drugs)
                    tmp_extra_drugs += tmp_drugs
                if len(tmp_extra_drugs) > 0:
                    tmp_drugs = []
                    for elem in tmp_bact.drugs:
                        tmp_drugs.append(elem.id_drug)
                    diffs = list_diff(tmp_extra_drugs, tmp_drugs)
                    print("adding extra drugs from drugbank by pathogens")
                    if len(diffs) > 0:
                        for elem in diffs:
                            if elem not in tmp_drugs:
                                tmp_bact.drugs.append(get_drug_kegg(elem))
            with open(f"db_patho/{tmp_bact.id_bact}.txt", "wb") as f:
                    pickle.dump(tmp_bact, f)
            bacts_db.append(tmp_bact)
    else:
        for filename in os.listdir("db_patho"):
            with open("db_patho/" + filename, "rb") as f:
                bacts_db.append(pickle.load(f))

if __name__ == "__main__":
    main()