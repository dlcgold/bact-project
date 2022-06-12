from kegg_helper import *
from scrape_pathways import *
from drugbank_sql_lite import *
from utils import *
from bact_classes import *
import matplotlib.pyplot as plt
import math
import numpy as np


Entrez.email = "d.cozzi@campus.unimib.it"

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
            if 'DRUG' in bact_tmp:
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
    fig.set_size_inches((12, 16), forward=False)
    fig.savefig(fname, dpi=1000)
    

    # Filling from genomejp by pathogen query
    bact_names = "" # for geo parsing
    for tmp_bact in bacts:
        bact_names += (' ' + tmp_bact.name.lower())
        tmp_extra_drugs = get_genomejp_drugs(tmp_bact.name)
        if len(tmp_extra_drugs) > 0:
            tmp_drugs = []
            for elem in tmp_bact.drugs:
                tmp_drugs.append(elem.id_drug)
            diffs = list_diff(tmp_extra_drugs, tmp_drugs)
            print("add extra drugs")
            if len(diffs) > 0:
                for elem in diffs:
                    if elem not in tmp_drugs:
                        tmp_bact.drugs.append(get_drug_kegg(elem))


if __name__ == "__main__":
    main()