"""
progetto-ds.ipynb
"""
from Bio import Entrez
from Bio.KEGG import REST
import pandas as pd
import plotly.express as px
from mordecai import Geoparser
import os
import pickle

Entrez.email = "d.cozzi@campus.unimib.it"


class paper:
    def __init__(self, pmid, authors, title, journal, doi):
        self.pmid = pmid
        self.authors = authors
        self.title = title
        self.journal = journal
        self.doi = doi

    def __repr__(self):
        return f"{self.title} ({self.pmid})\t{self.authors}\t{self.journal} ({self.doi})"


class bact:
    def __init__(self, name, id, description, category, sub, drugs, papers):
        self.name = name
        self.id = id
        self.drugs = drugs
        self.papers = papers
        self.description = description
        self.category = category
        self.sub = sub

    def __repr__(self):
        return f"{self.name} ({self.id})\n{self.drugs}\n ({self.papers})"


class drug:
    def __init__(self, name, id):
        self.name = name
        self.id = id

    def __repr__(self):
        return f"{self.name}, ({self.id})"


class geo_data:
    def __init__(self, id, text, geo_dict):
        self.id = id
        self.text = text
        self.geo_dict = geo_dict


def parse(data):
    sub_bool = False
    sub_count = 0
    bact_sub = ""
    drug_bool = False
    paper_bool = False
    id_tmp = ""
    title_tmp = ""
    authors_tmp = ""
    journal_tmp = ""
    doi_tmp = ""
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

        if sub_bool == False and spl[0] == "BRITE":
            sub_bool = True
            sub_count += 1
        elif sub_bool == True and sub_count == 1:
            sub_count += 1
        elif sub_bool == True and sub_count == 2:
            sub_bool = False
            bact_sub = "".join(spl[i] + " " for i in range(0, len(spl)))

        if drug_bool == False and spl[0] == "DRUG":
            name = "".join(spl[i] + " " for i in range(1, len(spl) - 1))
            drugs.append(drug(name, spl[-1].replace("[", "").replace("]", "")))
            drug_bool = True
        elif drug_bool == True and spl[0].upper() != spl[0]:
            name = "".join(spl[i] + " " for i in range(0, len(spl)))
            drugs.append(drug(name, spl[-1].replace("[", "").replace("]", "")))
        elif drug_bool == True and spl[0].upper() == spl[0]:
            drug_bool = False

        if paper_bool == False and spl[0] == "REFERENCE":
            if len(spl) > 1 and spl[1].split(":")[0] == "PMID":
                id_tmp = spl[1].split(":")[1].replace('\n', ' ')
            paper_bool = True
        elif paper_bool == True and spl[0] == "AUTHORS":
            authors_tmp = "".join(spl[i] + " " for i in range(1, len(spl)))
        elif paper_bool == True and spl[0] == "TITLE":
            title_tmp = "".join(spl[i] + " " for i in range(1, len(spl)))
        elif paper_bool == True and spl[0] == "JOURNAL":
            journal_tmp = "".join(spl[i] + " " for i in range(1, len(spl)))
            doi_bool = True
        elif paper_bool == True and doi_bool == True and spl[0].upper() != spl[0]:
            doi_tmp = spl[0].split(":")[1].replace('\n', ' ')
            doi_bool = False
            paper_bool = False
            papers.append(paper(id_tmp, authors_tmp, title_tmp, journal_tmp, doi_tmp))

    # print(drugs)
    # print(papers)

    bact_tmp = bact(bact_name, "ds:" + bact_id, bact_des, bact_cat, bact_sub, drugs, papers)
    # print(bact_tmp)
    return bact_tmp


def main():
    # Estrazione batteri
    sngOrg = REST.kegg_get("br:br08401").read()
    H_list = []
    check = True
    bact_bool = False
    for line in sngOrg.splitlines():
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
    no_paper = 0
    drug_list = []
    no_drug_list = []
    all_get = ""
    if not os.path.exists('kegg_get') or not os.listdir("kegg_get"):
        if not os.path.exists('kegg_get'):
            os.makedirs('kegg_get')
        for h_id in H_list:
            bact_tmp = REST.kegg_get(h_id).read()
            with open(f"kegg_get/{h_id}.txt", "w") as f:
                f.write(bact_tmp)
            all_get = all_get + bact_tmp
            all_get = all_get + "\n"
            if 'DRUG' in bact_tmp:
                drug_list.append(h_id)
            else:
                no_drug_list.append(h_id)
            if not 'PMID' in bact_tmp:
                no_paper = no_paper + 1

    # print(drug_list)
    # print(no_drug_list)
    # sngOrg = REST.kegg_get(["H00349"]).read()
    # print('DRUG' in sngOrg)

    # print(len(drug_list))
    # print(len(no_drug_list))
    # print(no_paper)

    # sngOrg = REST.kegg_get("gn:ype").read()

    sngOrg = REST.kegg_get("ds:H00111").read()

    # print(sngOrg)
    bact_tmp = parse(sngOrg)

    # print(bact)

    bacts = []
    for filename in os.listdir("kegg_get"):
        with open("kegg_get/" + filename, "r") as f:
            bacts.append(parse(str(f.read())))

    print(len(bacts))
    if not os.path.exists('abstract_get') or not os.listdir("abstract_get"):
        if not os.path.exists('abstract_get'):
            os.makedirs('abstract_get')
        for bact_tmp in bacts:
            for paper_tmp in bact_tmp.papers:
                if paper_tmp.pmid:
                    print(paper_tmp.pmid)
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
                    tmp = geo_data(id_file, text, geo.geoparse(text))
                    geo_list_abs.append(tmp)
                    with open(f"ser_abs/{tmp.id}.ser", "wb") as fw:
                        pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_abs"):
            with open(f"ser_abs/{filename}", "rb") as f:
                geo_list_abs.append(pickle.load(f))
    print(len(geo_list_abs))

    geo_df = pd.DataFrame(columns=["id", "lat", "lon", "type"])
    for geo_elem in geo_list_abs:
        if geo_elem.geo_dict:
            id = geo_elem.id
            for single_elem in geo_elem.geo_dict:
                print(single_elem)
                if 'geo' in single_elem:
                    geo_df.loc[len(geo_df.index)] = [id,
                                                     float(single_elem['geo']['lat']),
                                                     float(single_elem['geo']['lon']),
                                                     "Abstract"]
    print(geo_df)

    geo_list_title = []
    if not os.path.exists('ser_title'):
        os.makedirs('ser_title')
    if len(os.listdir("ser_title")) == 0:
        geo = Geoparser()
        for bact_tmp in bacts:
            count = 0
            for paper_tmp in bact_tmp.papers:
                tmp = geo_data(paper_tmp.pmid, paper_tmp.title, geo.geoparse(paper_tmp.title))
                geo_list_title.append(tmp)
                with open(f"ser_title/{paper_tmp.pmid}_{count}.ser", "wb") as fw:
                    pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_title"):
            with open(f"ser_title/{filename}", "rb") as f:
                geo_list_title.append(pickle.load(f))

    for geo_elem in geo_list_title:
        if geo_elem.geo_dict:
            id = geo_elem.id
            for single_elem in geo_elem.geo_dict:
                print(single_elem)
                if 'geo' in single_elem:
                    geo_df.loc[len(geo_df.index)] = [id,
                                                     float(single_elem['geo']['lat']),
                                                     float(single_elem['geo']['lon']),
                                                     "Title"]

    geo_list_bact = []
    if not os.path.exists('ser_bacts'):
        os.makedirs('ser_bacts')
    if len(bacts) != len(os.listdir("ser_bacts")):
        geo = Geoparser()
        for bact_tmp in bacts:
            tmp = geo_data(bact_tmp.id, bact_tmp.description, geo.geoparse(bact_tmp.description))
            geo_list_bact.append(tmp)
            with open(f"ser_bacts/{bact_tmp.id}.ser", "wb") as fw:
                pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_bacts"):
            with open(f"ser_bacts/{filename}", "rb") as f:
                geo_list_bact.append(pickle.load(f))
    for geo_elem in geo_list_bact:
        if geo_elem.geo_dict:
            id = geo_elem.id
            for single_elem in geo_elem.geo_dict:
                print(single_elem)
                if 'geo' in single_elem:
                    geo_df.loc[len(geo_df.index)] = [id,
                                                     float(single_elem['geo']['lat']),
                                                     float(single_elem['geo']['lon']),
                                                     "Description"]
    # print(geo_df)

    fig = px.scatter_mapbox(geo_df,
                            hover_name="id",
                            lat="lat",
                            lon="lon",
                            color="type",
                            size_max=15,
                            zoom=1,
                            width=1080,
                            height=720,
                            color_discrete_sequence=["red", "blue", "yellow"],
                            mapbox_style="open-street-map"
                            )
    fig.show()


if __name__ == "__main__":
    main()
