"""
progetto-ds.ipynb
"""

import os
import pickle

import requests as requests
from Bio import Entrez
from Bio.KEGG import REST
from mordecai import Geoparser
import pandas as pd
import plotly.express as px

Entrez.email = "d.cozzi@campus.unimib.it"


class Paper:
    def __init__(self, pmid, authors, title, journal, doi):
        self.pmid = pmid
        self.authors = authors
        self.title = title
        self.journal = journal
        self.doi = doi

    def __repr__(self):
        return f"{self.title} ({self.pmid})\t{self.authors}\t{self.journal} ({self.doi})"


class Bact:
    def __init__(self, name, id_bact, description, category, sub, drugs, papers):
        self.name = name
        self.id_bact = id_bact
        self.drugs = drugs
        self.papers = papers
        self.description = description
        self.category = category
        self.sub = sub

    def __repr__(self):
        return f"{self.name} ({self.id_bact})\n{self.drugs}\n ({self.papers})"


class Drug:
    def __init__(self, name, id_drug):
        self.name = name
        self.id_drug = id_drug

    def __repr__(self):
        return f"{self.name}, ({self.id_drug})"


class GeoData:
    def __init__(self, id_geo, text, geo_dict):
        self.id_geo = id_geo
        self.text = text
        self.geo_dict = geo_dict


def parse(data):
    """
    Function to parse Kegg data
    :param data: Kegg result as string
    :return: a Bact instance
    """
    sub_bool = False
    sub_count = 0
    bact_sub = ""
    drug_bool = False
    paper_bool = False
    id_tmp = ""
    title_tmp = ""
    authors_tmp = ""
    journal_tmp = ""
    doi_bool = False
    bact_name = ""
    bact_id = ""
    bact_cat = ""
    bact_des = ""
    drugs = []
    papers = []
    for line in data.splitlines():
        spl = line.split()
        if spl[0] == "NAME":
            bact_name = "".join(spl[i] + " " for i in range(1, len(spl)))
        if spl[0] == "ENTRY":
            bact_id = "".join(spl[i] + " " for i in range(1, len(spl))).split()[0]
        if spl[0] == "CATEGORY":
            bact_cat = "".join(spl[i] + " " for i in range(1, len(spl)))
        if spl[0] == "DESCRIPTION":
            bact_des = "".join(spl[i] + " " for i in range(1, len(spl)))

        if not sub_bool and spl[0] == "BRITE":
            sub_bool = True
            sub_count += 1
        elif sub_bool and sub_count == 1:
            sub_count += 1
        elif sub_bool and sub_count == 2:
            sub_bool = False
            bact_sub = "".join(spl[i] + " " for i in range(0, len(spl)))

        if not drug_bool and spl[0] == "DRUG":
            name = "".join(spl[i] + " " for i in range(1, len(spl) - 1))
            drugs.append(Drug(name, spl[-1].replace("[", "").replace("]", "")))
            drug_bool = True
        elif drug_bool and spl[0].upper() != spl[0]:
            name = "".join(spl[i] + " " for i in range(0, len(spl)))
            drugs.append(Drug(name, spl[-1].replace("[", "").replace("]", "")))
        elif drug_bool and spl[0].upper() == spl[0]:
            drug_bool = False

        if not paper_bool and spl[0] == "REFERENCE":
            if len(spl) > 1 and spl[1].split(":")[0] == "PMID":
                id_tmp = spl[1].split(":")[1].replace('\n', ' ')
            paper_bool = True
        elif paper_bool and spl[0] == "AUTHORS":
            authors_tmp = "".join(spl[i] + " " for i in range(1, len(spl)))
        elif paper_bool and spl[0] == "TITLE":
            title_tmp = "".join(spl[i] + " " for i in range(1, len(spl)))
        elif paper_bool and spl[0] == "JOURNAL":
            journal_tmp = "".join(spl[i] + " " for i in range(1, len(spl)))
            doi_bool = True
        elif paper_bool and doi_bool and spl[0].upper() != spl[0]:
            doi_tmp = spl[0].split(":")[1].replace('\n', ' ')
            doi_bool = False
            paper_bool = False
            papers.append(Paper(id_tmp, authors_tmp, title_tmp, journal_tmp, doi_tmp))

    # print(drugs)
    # print(papers)

    bact_tmp = Bact(bact_name, "ds:" + bact_id, bact_des, bact_cat, bact_sub, drugs, papers)
    # print(bact_tmp)
    return bact_tmp

def is_valid(item, bact_names):
    word = item['word']
    only_letters = word.replace(" ", "").isalpha()
    upper_case = word.isupper()
    no_cap_letters = all(w[0].islower() for w in word.split(' '))
    no_country = item['country_predicted'] == ''

    # Abstracts contain terms related to bacterias and species
    # so we exclude geographic location names containing these terms

    is_bacteria = bact_names.find(word.lower()) != -1
    is_species = False

    # Wikipedia pages related to species show an infobox containing
    # their scientific classification

    response = requests.get(url="https://en.wikipedia.org/wiki/" + word.replace(" ", "_"))
    if response.status_code == 200:
        is_species = response.text.find("Scientific classification") != -1

    valid = only_letters and not upper_case and not no_cap_letters and not no_country and not is_bacteria and not is_species

    return valid


def main():
    # Estrazione batteri
    print("getting bacterial infections list")
    bactetial_infections = REST.kegg_get("br:br08401").read()
    H_list = []
    check = True
    bact_bool = False
    for line in bactetial_infections.splitlines():
        if not check:
            break
        if line == "ABacterial infections":
            bact_bool = True
            continue
        if bact_bool:
            if line[0] == 'C':
                H_list.append(line.split()[1])
        if bact_bool and line[0] == 'A':
            check = False

    # sngOrg = REST.kegg_list("ds").read()
    # sngOrg = REST.kegg_get(["H02511", "H02512"]).read()
    # sngOrg = REST.kegg_get("br:br08401").read()
    # print(sngOrg)
    # json_data = json.loads(sngOrg)
    # H_list = []
    # print(json_data['children'][0]['children'][0]['children'][0]['name'])
    # print(len(json_data['children']))
    # for line in sngOrg.splitlines():
    #  if line[0] == 'C':
    #    H_list.append(line.split()[1])
    # H_list = []
    # for line in sngOrg.splitlines():
    #  code = line.split("\t")[0]
    #  code_split = code.split(":")
    #  if code_split[0] == "ds":
    #    H_list.append(code_split[1])
    # print(H_list)
    # print(len(H_list))
    print("Fetch data from KEGG")
    no_paper = 0
    drug_list = []
    no_drug_list = []
    all_get = ""
    if not os.path.exists('kegg_get') or not os.listdir("kegg_get"):
        if not os.path.exists('kegg_get'):
            os.makedirs('kegg_get')
        if not os.path.exists('kegg_get_drug'):
            os.makedirs('kegg_get_drug')
        for h_id in H_list:
            bact_tmp = REST.kegg_get(h_id).read()

            all_get = all_get + bact_tmp
            all_get = all_get + "\n"
            if 'DRUG' in bact_tmp:
                with open(f"kegg_get/{h_id}.txt", "w") as f:
                    f.write(bact_tmp)
                drug_list.append(h_id)
            else:
                with open(f"kegg_get_drug/{h_id}.txt", "w") as f:
                    f.write(bact_tmp)
                no_drug_list.append(h_id)
            if 'PMID' not in bact_tmp:
                no_paper = no_paper + 1

    # print(drug_list)
    # print(no_drug_list)
    # sngOrg = REST.kegg_get(["H00349"]).read()
    # print('DRUG' in sngOrg)

    # print(len(drug_list))
    # print(len(no_drug_list))
    # print(no_paper)

    # sngOrg = REST.kegg_get("gn:ype").read()

    # sngOrg = REST.kegg_get("ds:H00111").read()
    print("parse KEGG data")
    bacts = []
    bact_names = ""
    for filename in os.listdir("kegg_get"):
        with open("kegg_get/" + filename, "r") as f:
            tmp_bact = parse(str(f.read()))
            bact_names += (' ' + tmp_bact.name.lower())
            bacts.append(tmp_bact)

    print(f"Total: {len(bacts)} batteries")

    print("Extract geo_data from abstract")
    if not os.path.exists('abstract_get') or not os.listdir("abstract_get"):
        if not os.path.exists('abstract_get'):
            os.makedirs('abstract_get')
        for bact_tmp in bacts:
            for paper_tmp in bact_tmp.papers:
                if paper_tmp.pmid:
                    with open(f"abstract_get/{paper_tmp.pmid}.txt", "w") as f:
                        pubmed_entry = Entrez.efetch(db="pubmed",
                                                     id=paper_tmp.pmid,
                                                     retmode="xml")
                        result = Entrez.read(pubmed_entry)
                        article = result['PubmedArticle'][0]['MedlineCitation']['Article']
                        if 'Abstract' in article:
                            f.write(article['Abstract']['AbstractText'][0])
    geo_list_abs = []

    if not os.path.exists('ser_abs'):
        os.makedirs('ser_abs')
    if len(os.listdir("abstract_get")) != len(os.listdir("ser_abs")):
        geo = Geoparser()
        for filename in os.listdir("abstract_get"):
            id_file = filename.split(".")[0]
            if not os.path.isfile(f"ser_abs/{id_file}.ser"):
                with open("abstract_get/" + filename, "r") as f:
                    text = str(f.read())
                    # print(filename)
                    tmp = GeoData(id_file, text, geo.geoparse(text))
                    geo_list_abs.append(tmp)
                    with open(f"ser_abs/{tmp.id_geo}.ser", "wb") as fw:
                        pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_abs"):
            with open(f"ser_abs/{filename}", "rb") as f:
                geo_list_abs.append(pickle.load(f))
    # print(len(geo_list_abs))

    geo_df = pd.DataFrame(columns=["id", "lat", "lon", "type"])
    for geo_elem in geo_list_abs:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    geo_df.loc[len(geo_df.index)] = [id_tmp,
                                                     float(single_elem['geo']['lat']),
                                                     float(single_elem['geo']['lon']),
                                                     "Abstract"]

    print("Extract geo_data from titles")
    geo_list_title = []
    if not os.path.exists('ser_title'):
        os.makedirs('ser_title')
    if len(os.listdir("ser_title")) == 0:
        geo = Geoparser()
        for bact_tmp in bacts:
            count = 0
            for paper_tmp in bact_tmp.papers:
                tmp = GeoData(paper_tmp.pmid, paper_tmp.title, geo.geoparse(paper_tmp.title))
                geo_list_title.append(tmp)
                with open(f"ser_title/{paper_tmp.pmid}_{count}.ser", "wb") as fw:
                    pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_title"):
            with open(f"ser_title/{filename}", "rb") as f:
                geo_list_title.append(pickle.load(f))

    for geo_elem in geo_list_title:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    geo_df.loc[len(geo_df.index)] = [id_tmp,
                                                     float(single_elem['geo']['lat']),
                                                     float(single_elem['geo']['lon']),
                                                     "Title"]

    print("Extract geo data from descriptions")
    geo_list_bact = []
    if not os.path.exists('ser_bacts'):
        os.makedirs('ser_bacts')
    if len(bacts) != len(os.listdir("ser_bacts")):
        geo = Geoparser()
        for bact_tmp in bacts:
            tmp = GeoData(bact_tmp.id_bact, bact_tmp.description,
                          geo.geoparse(bact_tmp.description))
            geo_list_bact.append(tmp)
            with open(f"ser_bacts/{bact_tmp.id_bact}.ser", "wb") as fw:
                pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_bacts"):
            with open(f"ser_bacts/{filename}", "rb") as f:
                geo_list_bact.append(pickle.load(f))
    for geo_elem in geo_list_bact:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    geo_df.loc[len(geo_df.index)] = [id_tmp,
                                                     float(single_elem['geo']['lat']),
                                                     float(single_elem['geo']['lon']),
                                                     "Description"]

    print("parse KEGG data with drugs")
    bacts_drug = []
    for filename in os.listdir("kegg_get_drug"):
        with open("kegg_get_drug/" + filename, "r") as f:
            bacts_drug.append(parse(str(f.read())))
    print(f"Total: {len(bacts_drug)} batteries with drugs")

    print("Extract geo_data from abstract with drugs")
    if not os.path.exists('abstract_get_drug') or not os.listdir("abstract_get_drug"):
        if not os.path.exists('abstract_get_drug'):
            os.makedirs('abstract_get_drug')
        for bact_tmp in bacts_drug:
            for paper_tmp in bact_tmp.papers:
                if paper_tmp.pmid:
                    with open(f"abstract_get_drug/{paper_tmp.pmid}.txt", "w") as f:
                        pubmed_entry = Entrez.efetch(db="pubmed",
                                                     id=paper_tmp.pmid,
                                                     retmode="xml")
                        result = Entrez.read(pubmed_entry)
                        article = result['PubmedArticle'][0]['MedlineCitation']['Article']
                        if 'Abstract' in article:
                            f.write(article['Abstract']['AbstractText'][0])

    geo_list_abs_drug = []
    if not os.path.exists('ser_abs_drug'):
        os.makedirs('ser_abs_drug')
    if len(os.listdir("abstract_get_drug")) != len(os.listdir("ser_abs_drug")):
        geo = Geoparser()
        for filename in os.listdir("abstract_get_drug"):
            id_file = filename.split(".")[0]
            if not os.path.isfile(f"ser_abs_drug/{id_file}.ser"):
                with open("abstract_get_drug/" + filename, "r") as f:
                    text = str(f.read())
                    tmp = GeoData(id_file, text, geo.geoparse(text))
                    geo_list_abs_drug.append(tmp)
                    with open(f"ser_abs_drug/{tmp.id_geo}.ser", "wb") as fw:
                        pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_abs_drug"):
            with open(f"ser_abs_drug/{filename}", "rb") as f:
                geo_list_abs_drug.append(pickle.load(f))
    # print(len(geo_list_abs))

    geo_df_drug = pd.DataFrame(columns=["id", "lat", "lon", "type"])
    for geo_elem in geo_list_abs_drug:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    geo_df_drug.loc[len(geo_df_drug.index)] = [id_tmp,
                                                               float(single_elem['geo']['lat']),
                                                               float(single_elem['geo']['lon']),
                                                               "Abstract"]

    print("Extract geo_data from titles  with drugs")
    geo_list_title_drug = []
    if not os.path.exists('ser_title_drug'):
        os.makedirs('ser_title_drug')
    if len(os.listdir("ser_title_drug")) == 0:
        geo = Geoparser()
        for bact_tmp in bacts_drug:
            count = 0
            for paper_tmp in bact_tmp.papers:
                tmp = GeoData(paper_tmp.pmid, paper_tmp.title, geo.geoparse(paper_tmp.title))
                geo_list_title_drug.append(tmp)
                with open(f"ser_title_drug/{paper_tmp.pmid}_{count}.ser", "wb") as fw:
                    pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_title_drug"):
            with open(f"ser_title_drug/{filename}", "rb") as f:
                geo_list_title_drug.append(pickle.load(f))

    for geo_elem in geo_list_title_drug:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    geo_df_drug.loc[len(geo_df_drug.index)] = [id_tmp,
                                                               float(single_elem['geo']['lat']),
                                                               float(single_elem['geo']['lon']),
                                                               "Title"]

    print("Extract geo data from descriptions with drugs")
    geo_list_bact_drug = []
    if not os.path.exists('ser_bacts_drug'):
        os.makedirs('ser_bacts_drug')
    if len(bacts_drug) != len(os.listdir("ser_bacts_drug")):
        geo = Geoparser()
        for bact_tmp in bacts_drug:
            tmp = GeoData(bact_tmp.id_bact, bact_tmp.description,
                          geo.geoparse(bact_tmp.description))
            geo_list_bact_drug.append(tmp)
            with open(f"ser_bacts_drug/{bact_tmp.id_bact}.ser", "wb") as fw:
                pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_bacts_drug"):
            with open(f"ser_bacts_drug/{filename}", "rb") as f:
                geo_list_bact_drug.append(pickle.load(f))
    for geo_elem in geo_list_bact_drug:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    geo_df_drug.loc[len(geo_df_drug.index)] = [id_tmp,
                                                               float(single_elem['geo']['lat']),
                                                               float(single_elem['geo']['lon']),
                                                               "Description"]

    geo_df_total = pd.DataFrame(columns=["id", "lat", "lon", "type", "Drug"])
    for index, row in geo_df.iterrows():
        geo_df_total.loc[len(geo_df_total.index)] = [row.id,
                                                     row.lat,
                                                     row.lon,
                                                     row.type,
                                                     "No Drug"]
    for index, row in geo_df_drug.iterrows():
        geo_df_total.loc[len(geo_df_total.index)] = [row.id,
                                                     row.lat,
                                                     row.lon,
                                                     row.type,
                                                     "Drug"]
    print("Produce Maps")
    fig = px.scatter_mapbox(geo_df,
                            hover_name="id",
                            lat="lat",
                            lon="lon",
                            color="type",
                            title="Without drugs",
                            # size_max=15,
                            zoom=1,
                            width=1080,
                            height=720,
                            color_discrete_sequence=["red", "blue", "green"],
                            mapbox_style="open-street-map"
                            )
    fig.show()

    fig = px.scatter_mapbox(geo_df_drug,
                            hover_name="id",
                            lat="lat",
                            lon="lon",
                            title="With drugs",
                            color="type",
                            # size_max=15,
                            zoom=1,
                            width=1080,
                            height=720,
                            color_discrete_sequence=["red", "blue", "green"],
                            mapbox_style="open-street-map"
                            )
    fig.show()

    fig = px.scatter_mapbox(geo_df_total,
                            hover_name="id",
                            lat="lat",
                            lon="lon",
                            title="Total",
                            color="Drug",
                            # size_max=15,
                            zoom=1,
                            width=1080,
                            height=720,
                            color_discrete_sequence=["red", "green"],
                            mapbox_style="open-street-map"
                            )
    fig.show()


if __name__ == "__main__":
    main()
