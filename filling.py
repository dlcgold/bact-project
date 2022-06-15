import math

from wordcloud import WordCloud
import matplotlib.pyplot as plt
import numpy as np
from gensim.parsing.preprocessing import remove_stopwords

from genomejp import *
from geo_parsing import *
from kegg_helper import *

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

    # display map of diseases without drugs
    bacts_nodrug_df = geo_parse(bacts, "nodrug")
    bacts_nodrug_map_display = display_map(bacts_nodrug_df,
                                           "Locations associated to bacteria without drugs",
                                           "type",
                                           ["red", "blue", "green"])

    # Filling list with disease with drugs
    print("parse KEGG data with drugs")
    bacts_drug = []
    for filename in os.listdir("kegg_get_drug"):
        with open("kegg_get_drug/" + filename, "r") as f:
            tmp_bact = parse(str(f.read()))
            bacts_drug.append(tmp_bact)
    print(f"{len(bacts_drug)} bacterial infections found w/ drugs")

    # display map of diseases with drugs
    bacts_drug_df = geo_parse(bacts_drug, "drug")
    bacts_drug_map_display = display_map(bacts_drug_df,
                                         "Locations associated to bacteria with drugs", "type",
                                         ["red", "blue", "green"])

    geo_df_total = merge_geo_df(bacts_nodrug_df, bacts_drug_df, "no drug", "drug")
    bacts_all_map_display = display_map(geo_df_total,
                                        "Locations associated to all bacteria",
                                        "type",
                                        ["red", "blue", "green"])

    # bar plot by quantity of infection by subgroups

    # Data for the plots
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

    # charge np array for the bar and plot them
    labels = []
    keys = list(set(list(bar_data.keys()) + list(bar_data_drugs.keys())))
    np_drugs = np.zeros(len(keys))
    np_no_drugs = np.zeros(len(keys))
    i = 0
    for key in keys:
        labels.append(key.replace("Infections caused by", "").strip())
        if key in bar_data.keys():
            np_no_drugs[i] = bar_data[key]
        if key in bar_data_drugs.keys():
            np_drugs[i] = bar_data_drugs[key]
        i += 1
    delta = np.arange(len(labels))
    width = 0.4
    fig, ax = plt.subplots()
    ax.bar(delta, np_no_drugs, color='red', width=width, label="Without drug")
    ax.bar(delta + width, np_drugs, color='green', width=width, label="With drug")
    plt.title("Bacterial infections with/without drugs by subgroup")
    plt.ylabel("Quantity of infections")
    plt.xticks(delta + width / 2, labels, rotation=45, fontsize=7)
    y_int = range(0, math.ceil(max(np_no_drugs) + 1))
    plt.yticks(y_int)
    ax.legend()
    if not os.path.exists('plot_print'):
        os.makedirs('plot_print')
    fname = './plot_print/comparative_drugs_bar_plot.png'
    plt.show()
    fig.set_size_inches((16, 12), forward=False)
    fig.savefig(fname, dpi=500)

    # Filling from genomejp by pathogen query
    print("Filling from genomejp by pathogen query")
    bact_names = ""  # for geo parsing
    bacts_genomejp_patho = []
    if not os.path.exists('genomejp_patho') or not os.listdir("genomejp_patho"):
        if not os.path.exists('genomejp_patho'):
            os.makedirs('genomejp_patho')
        for tmp_bact in bacts:
            bact_names += (' ' + tmp_bact.name.lower())
            tmp_extra_drugs = get_drugs_for_disease_name(tmp_bact.name)
            if len(tmp_extra_drugs) > 0:
                tmp_drugs = []
                for elem in tmp_bact.drugs:
                    tmp_drugs.append(elem.id_drug)
                diffs = list_diff(tmp_extra_drugs, tmp_drugs)
                print("adding extra drugs from genomejp by pathogen")
                if len(diffs) > 0:
                    for elem in diffs:
                        if elem not in tmp_drugs:
                            tmp_drug = get_drug_kegg(elem)
                            tmp_drug.origin = "Genomejp Pathogen"
                            tmp_bact.drugs.append(tmp_drug)
            with open(f"genomejp_patho/{tmp_bact.id_bact}.txt", "wb") as f:
                pickle.dump(tmp_bact, f)
            bacts_genomejp_patho.append(tmp_bact)
    else:
        for filename in os.listdir("genomejp_patho"):
            with open("genomejp_patho/" + filename, "rb") as f:
                bacts_genomejp_patho.append(pickle.load(f))

    no_drug_after_genomejp_patho = []
    for bact in bacts_genomejp_patho:
        if len(bact.drugs) == 0:
            no_drug_after_genomejp_patho.append(bact)
    print(f"{len(no_drug_after_genomejp_patho)} bacterial infections found w/out drugs after "
          f"genomejp by pathogen")

    # proceed to plot iff there was some filling
    if len(no_drug_after_genomejp_patho) != len(bacts):
        # bar plot by quantity of infection by subgroups after genomejp
        bar_data_genome_patho = {}

        for bact in bacts_genomejp_patho:
            if bact.sub not in bar_data_genome_patho.keys():
                bar_data_genome_patho[bact.sub] = 0
            if len(bact.drugs) > 0:
                bar_data_genome_patho[bact.sub] += 1
        bar_data_no_genome = {}
        for key, value in bar_data_genome_patho.items():
            bar_data_no_genome[key] = bar_data[key] - value

        fig, ax = plt.subplots()
        np_no_genome_values = np.zeros(len(keys))
        np_drugs_values = np.zeros(len(keys))
        np_genome_values = np.zeros(len(keys))
        labels_genome = []
        i = 0
        for key in keys:
            labels_genome.append(key.replace("Infections caused by", "").strip())
            if key in bar_data_no_genome.keys():
                np_no_genome_values[i] = bar_data_no_genome[key]
            if key in bar_data_drugs.keys():
                np_drugs_values[i] = bar_data_drugs[key]
            if key in bar_data_genome_patho.keys():
                np_genome_values[i] = bar_data_genome_patho[key]
            i += 1

        ax.bar(delta, np_no_genome_values, color='red', width=width, label="Without drug")
        ax.bar(delta + width, np_drugs_values, color='green', width=width,
               label="With drug")
        ax.bar(delta + width, np_genome_values, bottom=np_drugs_values,
               color='orange', width=width, label="Filled with genomejp")
        plt.title("Bacterial infections with/without drugs by subgroup after\n"
                  "filling with genomejp by pathogen")
        plt.ylabel("Quantity of infections")
        plt.xticks(delta + width / 2, labels_genome, rotation=45, fontsize=7)
        y_int = range(0,
                      math.ceil(max(np_drugs_values + np_genome_values)) + 1)
        plt.yticks(y_int)
        ax.legend()
        fname = './plot_print/comparative_drugs_bar_plot_after_genomejp_patho.png'
        plt.show()
        fig.set_size_inches((16, 12), forward=False)
        fig.savefig(fname, dpi=500)

        # update with/without drug data
        for key in keys:
            if key in bar_data_genome_patho.keys():
                if key not in bar_data_drugs.keys():
                    bar_data_drugs[key] = bar_data_genome_patho[key]
                else:
                    bar_data_drugs[key] = bar_data_drugs[key] + bar_data_genome_patho[key]
        for key in keys:
            if key in bar_data.keys() and key in bar_data_genome_patho.keys():
                bar_data[key] = bar_data[key] - bar_data_genome_patho[key]
            if key in bar_data.keys() and bar_data[key] == 0:
                bar_data.pop(key)

    # update bacts list
    bacts = no_drug_after_genomejp_patho

    print("Filling from genomejp by pathway query")
    bact_names_path = ""  # for geo parsing
    bacts_genomejp_path = []
    if not os.path.exists('genomejp_path') or not os.listdir("genomejp_path"):
        if not os.path.exists('genomejp_path'):
            os.makedirs('genomejp_path')
        for tmp_bact in bacts:
            bact_names_path += (' ' + tmp_bact.name.lower())
            tmp_extra_drugs = []
            for path in tmp_bact.pathways:
                tmp_extra_drugs += get_drugs_for_pathway_kegg(path.id_path)
            if len(tmp_extra_drugs) > 0:
                tmp_drugs = []
                for elem in tmp_bact.drugs:
                    tmp_drugs.append(elem.id_drug)
                diffs = list_diff(tmp_extra_drugs, tmp_drugs)
                print("adding extra drugs from genomejp by pathway")
                if len(diffs) > 0:
                    for elem in diffs:
                        if elem not in tmp_drugs:
                            tmp_drug = get_drug_kegg(elem)
                            tmp_drug.origin = "Origin Pathway"
                            tmp_bact.drugs.append(tmp_drug)
            with open(f"genomejp_path/{tmp_bact.id_bact}.txt", "wb") as f:
                pickle.dump(tmp_bact, f)
            bacts_genomejp_path.append(tmp_bact)
    else:
        for filename in os.listdir("genomejp_path"):
            with open("genomejp_path/" + filename, "rb") as f:
                bacts_genomejp_path.append(pickle.load(f))
    # print(len(bacts))
    no_drug_after_genomejp_path = []
    for bact in bacts_genomejp_path:
        if len(bact.drugs) == 0:
            no_drug_after_genomejp_path.append(bact)
    print(f"{len(no_drug_after_genomejp_path)} {type_infection} found w/out drugs after genomejp "
          f"pathway")

    # proceed to plot iff there was some filling
    if len(no_drug_after_genomejp_path) != len(bacts):
        # bar plot by quantity of infection by subgroups after genomejp
        bar_data_genome_path = {}
        for bact in bacts_genomejp_path:
            if bact.sub not in bar_data_genome_path.keys():
                bar_data_genome_path[bact.sub] = 0
            if len(bact.drugs) > 0:
                bar_data_genome_path[bact.sub] += 1
        bar_data_no_genome = {}
        for key, value in bar_data_genome_path.items():
            bar_data_no_genome[key] = bar_data[key] - value

        fig, ax = plt.subplots()
        np_no_genome_values = np.zeros(len(keys))
        np_drugs_values = np.zeros(len(keys))
        np_genome_values = np.zeros(len(keys))
        labels_genome = []
        i = 0
        for key in keys:
            labels_genome.append(key.replace("Infections caused by", "").strip())
            if key in bar_data_no_genome.keys():
                np_no_genome_values[i] = bar_data_no_genome[key]
            if key in bar_data_drugs.keys():
                np_drugs_values[i] = bar_data_drugs[key]
            if key in bar_data_genome_path.keys():
                np_genome_values[i] = bar_data_genome_path[key]
            i += 1

        ax.bar(delta, np_no_genome_values, color='red', width=width, label="Without drug")
        ax.bar(delta + width, np_drugs_values, color='green', width=width,
               label="With drug")
        ax.bar(delta + width, np_genome_values, bottom=np_drugs_values,
               color='orange', width=width, label="Filled with genomejp")
        plt.title("Bacterial infections with/without drugs by subgroup after\n"
                  "filling with genomejp by pathway")
        plt.ylabel("Quantity of infections")
        plt.xticks(delta + width / 2, labels_genome, rotation=45, fontsize=7)
        y_int = range(0,
                      math.ceil(max(np_drugs_values + np_genome_values)) + 1)
        plt.yticks(y_int)
        ax.legend()
        fname = './plot_print/comparative_drugs_bar_plot_after_genomejp_path.png'
        plt.show()
        fig.set_size_inches((16, 12), forward=False)
        fig.savefig(fname, dpi=500)

        for key in keys:
            if key in bar_data_genome_path.keys():
                if key not in bar_data_drugs.keys():
                    bar_data_drugs[key] = bar_data_genome_path[key]
                else:
                    bar_data_drugs[key] = bar_data_drugs[key] + bar_data_genome_path[key]
        for key in keys:
            if key in bar_data.keys() and key in bar_data_genome_path.keys():
                bar_data[key] = bar_data[key] - bar_data_genome_path[key]
            if key in bar_data.keys() and bar_data[key] == 0:
                bar_data.pop(key)
    bacts = no_drug_after_genomejp_path

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
                        ## Avoid 400 bad requests
                        try:
                            id_tmp = conv_db_id_to_kegg_id(tmp)
                        except:
                            id_tmp = None
                        ##
                        print(id_tmp)
                        if id_tmp != "None":
                            tmp_drugs.append(id_tmp)
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
                                tmp_drug = get_drug_kegg(elem)
                                tmp_drug.origin = "DrugBank pathogen"
                                tmp_bact.drugs.append(tmp_drug)
            with open(f"db_patho/{tmp_bact.id_bact}.txt", "wb") as f:
                pickle.dump(tmp_bact, f)
            bacts_db.append(tmp_bact)
    else:
        for filename in os.listdir("db_patho"):
            with open("db_patho/" + filename, "rb") as f:
                bacts_db.append(pickle.load(f))
    bacts_db_without_drugs = []

    for bact_tmp in bacts_db:
        if len(bact_tmp.drugs) == 0:
            bacts_db_without_drugs.append(bact_tmp)
    print(
        f"{len(bacts_db_without_drugs)} {type_infection} found w/out drugs after drugbank filling "
        f"by pathogen")

    # proceed to plot iff there was some filling
    if len(bacts_db_without_drugs) != len(bacts):
        bar_data_db = {}
        for bact in bacts_db:
            if bact.sub not in bar_data_db.keys():
                bar_data_db[bact.sub] = 0
            if len(bact.drugs) > 0:
                bar_data_db[bact.sub] += 1
        bar_data_no_db = {}
        for key, value in bar_data_db.items():
            bar_data_no_db[key] = bar_data[key] - value

        fig, ax = plt.subplots()
        np_no_db_values = np.zeros(len(keys))
        np_drugs_values = np.zeros(len(keys))
        np_db_values = np.zeros(len(keys))
        labels_db = []
        i = 0
        for key in keys:
            labels_db.append(key.replace("Infections caused by", "").strip())
            if key in bar_data_no_db.keys():
                np_no_db_values[i] = bar_data_no_db[key]
            if key in bar_data_drugs.keys():
                np_drugs_values[i] = bar_data_drugs[key]
            if key in bar_data_db.keys():
                np_db_values[i] = bar_data_db[key]
            i += 1
        ax.bar(delta, np_no_db_values, color='red', width=width, label="Without drug")
        ax.bar(delta + width, np_drugs_values, color='green', width=width,
               label="With drug")
        ax.bar(delta + width, np_db_values, bottom=np_drugs_values,
               color='orange', width=width, label="Filled with db")
        plt.title("Bacterial infections with/without drugs by subgroup after\n"
                  "filling with drugbank by pathogen")
        plt.ylabel("Quantity of infections")
        plt.xticks(delta + width / 2, labels_db, rotation=45, fontsize=7)
        y_int = range(0,
                      math.ceil(max(np_drugs_values + np_db_values)) + 1)
        plt.yticks(y_int)
        ax.legend()
        fname = './plot_print/comparative_drugs_bar_plot_after_db_patho.png'
        plt.show()
        fig.set_size_inches((16, 12), forward=False)
        fig.savefig(fname, dpi=500)
        for key in keys:
            if key in bar_data_db.keys():
                if key not in bar_data_drugs.keys():
                    bar_data_drugs[key] = bar_data_db[key]
                else:
                    bar_data_drugs[key] = bar_data_drugs[key] + bar_data_db[key]

        for key in keys:
            if key in bar_data.keys() and key in bar_data_db.keys():
                bar_data[key] = bar_data[key] - bar_data_db[key]
            if key in bar_data.keys() and bar_data[key] == 0:
                bar_data.pop(key)
    bacts = bacts_db_without_drugs

    print("Filling from drugbank by pathways")
    bacts_db_path = []
    if not os.path.exists('db_path') or not os.listdir("db_path"):
        if not os.path.exists('db_path'):
            os.makedirs('db_path')
        for tmp_bact in bacts:
            if len(tmp_bact.pathways) > 0:
                tmp_extra_drugs = []
                for path_tmp in tmp_bact.pathways:
                    db_drugs = get_drugs_for_pathway(path_tmp)
                    # print(db_drugs)
                    tmp_drugs = []
                    for tmp in db_drugs:
                        print(tmp)
                        ## Aggiunto per evitare 400 bad requests
                        try:
                            id_tmp = conv_db_id_to_kegg_id(tmp)
                        except:
                            id_tmp = None
                        ##
                        print(id_tmp)
                        if id_tmp != "None":
                            tmp_drugs.append(id_tmp)
                    tmp_extra_drugs += tmp_drugs
                if len(tmp_extra_drugs) > 0:
                    tmp_drugs = []
                    for elem in tmp_bact.drugs:
                        tmp_drugs.append(elem.id_drug)
                    diffs = list_diff(tmp_extra_drugs, tmp_drugs)
                    print("adding extra drugs from drugbank by pathways")
                    if len(diffs) > 0:
                        for elem in diffs:
                            if elem not in tmp_drugs:
                                tmp_drug = get_drug_kegg(elem)
                                tmp_drug.origin = "DrugBank Pathway"
                                tmp_bact.drugs.append(tmp_drug)
            with open(f"db_path/{tmp_bact.id_bact}.txt", "wb") as f:
                pickle.dump(tmp_bact, f)
            bacts_db_path.append(tmp_bact)
    else:
        for filename in os.listdir("db_path"):
            with open("db_path/" + filename, "rb") as f:
                bacts_db_path.append(pickle.load(f))
    bacts_db_path_without_drugs = []

    for bact_tmp in bacts_db_path:
        if len(bact_tmp.drugs) == 0:
            bacts_db_path_without_drugs.append(bact_tmp)
    print(f"{len(bacts_db_path_without_drugs)} {type_infection} found w/out drugs after drugbank "
          f"filling by pathways")

    # proceed to plot iff there was some filling
    if len(bacts_db_path_without_drugs) != len(bacts):
        bar_data_db_path = {}
        for bact in bacts_db_path:
            if bact.sub not in bar_data_db_path.keys():
                bar_data_db_path[bact.sub] = 0
            if len(bact.drugs) > 0:
                bar_data_db_path[bact.sub] += 1
        bar_data_no_db_path = {}
        for key, value in bar_data_db_path.items():
            bar_data_no_db_path[key] = bar_data[key] - value

        fig, ax = plt.subplots()
        np_no_db_values = np.zeros(len(keys))
        np_drugs_values = np.zeros(len(keys))
        np_db_values = np.zeros(len(keys))
        labels_db = []
        i = 0
        for key in keys:
            labels_db.append(key.replace("Infections caused by", "").strip())
            if key in bar_data_no_db_path.keys():
                np_no_db_values[i] = bar_data_no_db_path[key]
            if key in bar_data_drugs.keys():
                np_drugs_values[i] = bar_data_drugs[key]
            if key in bar_data_db_path.keys():
                np_db_values[i] = bar_data_db_path[key]
            i += 1
        ax.bar(delta, np_no_db_values, color='red', width=width, label="Without drug")
        ax.bar(delta + width, np_drugs_values, color='green', width=width,
               label="With drug")
        ax.bar(delta + width, np_db_values, bottom=np_drugs_values,
               color='orange', width=width, label="Filled with db")
        plt.title("Bacterial infections with/without drugs by subgroup after\n"
                  "filling with drugbank by pathways")
        plt.ylabel("Quantity of infections")
        plt.xticks(delta + width / 2, labels_db, rotation=45, fontsize=7)
        y_int = range(0,
                      math.ceil(max(np_drugs_values + np_db_values)) + 1)
        plt.yticks(y_int)
        ax.legend()
        fname = './plot_print/comparative_drugs_bar_plot_after_db.png'
        plt.show()
        fig.set_size_inches((16, 12), forward=False)
        fig.savefig(fname, dpi=500)
        for key in keys:
            if key in bar_data_db_path.keys():
                if key not in bar_data_drugs.keys():
                    bar_data_drugs[key] = bar_data_db_path[key]
                else:
                    bar_data_drugs[key] = bar_data_drugs[key] + bar_data_db_path[key]
        # TODO CONTROLLARE TUTTI QUESTI FOR
        for key in keys:
            if key in bar_data_db_path.keys() and key in bar_data.keys():
                bar_data[key] = bar_data[key] - bar_data_db_path[key]
            if key in bar_data.keys() and bar_data[key] == 0:
                bar_data.pop(key)
    bacts = bacts_db_path_without_drugs

    # some wordclouds
    word_cloud_name = ""
    word_cloud_des = ""
    word_cloud_pap = ""
    word_cloud_abs_no_drug = get_full_abstract("nodrug")
    word_cloud_abs_drug = get_full_abstract("drug")
    bad_words = ["infection", "Infection", "human", "Human", "pathogen", "Pathogen",
                 "Bacterium", "bacterium", "Disease", "disease", "Clinical", "clinical"]
    for bact_tmp in bacts:
        name = bact_tmp.name
        des = bact_tmp.description
        pap = ""
        for paper in bact_tmp.papers:
            pap += paper.title
        for bw in bad_words:
            name = name.replace(bw, "")
            des = des.replace(bw, "")
            pap = pap.replace(bw, "")
        word_cloud_name += (name + " ")
        word_cloud_des += (des + " ")
        word_cloud_pap += (pap + " ")

    for bw in bad_words:
        word_cloud_abs_drug.replace(bw, "")
        word_cloud_abs_no_drug.replace(bw, "")

    wordcloud = WordCloud(max_font_size=40, background_color="white", contour_color='#5e81ac',
                          contour_width=0.1, scale=2).generate(remove_stopwords(word_cloud_name))
    plt.figure()
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis("off")
    plt.savefig("plot_print/name_wc.png", dpi=500)
    plt.close()

    wordcloud = WordCloud(max_font_size=40, background_color="white", contour_color='#5e81ac',
                          contour_width=0.1, scale=2).generate(remove_stopwords(word_cloud_des))
    plt.figure()
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis("off")
    plt.savefig("plot_print/des_wc.png", dpi=500)

    wordcloud = WordCloud(max_font_size=40, background_color="white", contour_color='#5e81ac',
                          contour_width=0.1, scale=2).generate(remove_stopwords(word_cloud_pap))
    plt.figure()
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis("off")
    plt.savefig("plot_print/pap_wc.png", dpi=500)

    wordcloud = WordCloud(max_font_size=40, background_color="white", contour_color='#5e81ac',
                          contour_width=0.1, scale=2).generate(
        remove_stopwords(word_cloud_abs_no_drug))
    plt.figure()
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis("off")
    plt.savefig("plot_print/no_drug_wc.png", dpi=500)

    wordcloud = WordCloud(max_font_size=40, background_color="white", contour_color='#5e81ac',
                          contour_width=0.1, scale=2).generate(
        remove_stopwords(word_cloud_abs_drug))
    plt.figure()
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis("off")
    plt.savefig("plot_print/drug_wc.png", dpi=500)


if __name__ == "__main__":
    main()
