import pickle
import os
import requests as requests
from Bio import Entrez
from Bio.KEGG import REST
from mordecai import Geoparser
import pandas as pd
from bact_classes import *
from utils import *
import plotly.express as px


"""    filtered_bact_drug = []
    filtered_abs_drug = []
    filtered_title_drug = []"""

def geo_parse(bacts):
    # geo locations filtered
    filtered_bact = []
    filtered_abs = []
    filtered_title = []

    bact_names = ""
    for tmp_bact in bacts:
        bact_names += (' ' + tmp_bact.name.lower())

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
                    tmp = GeoData(id_file, text, geo.geoparse(text))
                    geo_list_abs.append(tmp)
                    with open(f"ser_abs/{tmp.id_geo}.ser", "wb") as fw:
                        pickle.dump(tmp, fw)
    else:
        for filename in os.listdir("ser_abs"):
            with open(f"ser_abs/{filename}", "rb") as f:
                geo_list_abs.append(pickle.load(f))

    geo_df = pd.DataFrame(columns=["id", "lat", "lon", "type"])
    for geo_elem in geo_list_abs:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    if not is_valid(single_elem, bact_names):
                        filtered_abs.append(single_elem)
                    else:
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
                    if not is_valid(single_elem, bact_names):
                        filtered_title.append(single_elem)
                    else:
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
                    if not is_valid(single_elem, bact_names):
                        filtered_bact.append(single_elem)
                    else:
                        geo_df.loc[len(geo_df.index)] = [id_tmp,
                                                         float(single_elem['geo']['lat']),
                                                         float(single_elem['geo']['lon']),
                                                         "Description"]

    return geo_df

def merge_geo_df(no_drug, drug):
    geo_df_total = pd.DataFrame(columns=["id", "lat", "lon", "type", "drug"])
    for index, row in no_drug.iterrows():
        geo_df_total.loc[len(geo_df_total.index)] = [row.id,
                                                     row.lat,
                                                     row.lon,
                                                     row.type,
                                                     "no drug"]
    for index, row in drug.iterrows():
        geo_df_total.loc[len(geo_df_total.index)] = [row.id,
                                                     row.lat,
                                                     row.lon,
                                                     row.type,
                                                     "drug"]
    return geo_df_total

def display_map(geo_df, title, color, color_discrete_sequence):
    print(f"Producing map: {title}")
    fig = px.scatter_mapbox(geo_df,
                            hover_name="id",
                            lat="lat",
                            lon="lon",
                            color=color,
                            title=title,
                            # size_max=15,
                            zoom=1,
                            width=1080,
                            height=720,
                            color_discrete_sequence=color_discrete_sequence,
                            mapbox_style="open-street-map"
                            )
    fig.show()

    return fig
