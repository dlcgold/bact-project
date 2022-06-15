import os
import pickle
from matplotlib import pyplot as plt

import pandas as pd
import plotly.express as px
from Bio import Entrez
from mordecai import Geoparser

from utils import *

"""    filtered_bact_drug = []
    filtered_abs_drug = []
    filtered_title_drug = []"""


def geo_parse(bacts, type_print=""):
    geo = Geoparser()
    # geo locations filtered
    filtered_bact = []
    filtered_abs = []
    filtered_title = []

    bact_names = ""
    for tmp_bact in bacts:
        bact_names += (' ' + tmp_bact.name.lower())

    path_abs = f"abstract_get_{type_print}"
    path_ser_abs = f"ser_abs_{type_print}"
    path_ser_title = f"ser_title_{type_print}"
    path_ser_bacts = f"ser_bacts_{type_print}"
    print("Extract geo_data from abstract")

    if not os.path.exists(path_abs) or not os.listdir(path_abs):
        if not os.path.exists(path_abs):
            os.makedirs(path_abs)
        for bact_tmp in bacts:
            bact_name = bact_tmp.name.strip()
            if bact_name[-1] == ';':
                bact_name = bact_name[:-1]
            for paper_tmp in bact_tmp.papers:
                if paper_tmp.pmid:
                    with open(f"abstract_get_{type_print}/{paper_tmp.pmid}.txt", "w") as f:
                        pubmed_entry = Entrez.efetch(db="pubmed",
                                                     id=paper_tmp.pmid,
                                                     retmode="xml")
                        result = Entrez.read(pubmed_entry)
                        article = result['PubmedArticle'][0]['MedlineCitation']['Article']
                        if 'Abstract' in article:
                            f.write(f"{bact_name}\n{article['Abstract']['AbstractText'][0]}")
    geo_list_abs = []
    if not os.path.exists(path_ser_abs):
        os.makedirs(path_ser_abs)
    if len(os.listdir(path_abs)) == 0:
        for filename in os.listdir(path_abs):
            id_file = filename.split(".")[0]
            if not os.path.isfile(f"ser_abs_{type_print}/{id_file}.ser"):
                with open(f"abstract_get_{type_print}/" + filename, "r") as f:
                    lines = f.readlines()
                    if len(lines) >= 2:
                        text = " ".join(lines[1:])
                        tmp = GeoData(id_file, text, geo.geoparse(text), lines[0])
                        geo_list_abs.append(tmp)
                        with open(f"ser_abs_{type_print}/{tmp.id_geo}.ser", "wb") as fw:
                            pickle.dump(tmp, fw)

    else:
        for filename in os.listdir(path_ser_abs):
            with open(f"ser_abs_{type_print}/{filename}", "rb") as f:
                geo_list_abs.append(pickle.load(f))
    geo_df = pd.DataFrame(columns=["id", "name", "lat", "lon", "type"])
    for geo_elem in geo_list_abs:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            name_tmp = geo_elem.name.strip()
            if name_tmp[-1] == ';':
                name_tmp = name_tmp[:-1]
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    if not is_valid(single_elem, bact_names):
                        filtered_abs.append(single_elem)
                    else:
                        geo_df.loc[len(geo_df.index)] = [id_tmp,
                                                         name_tmp,
                                                         float(single_elem['geo']['lat']),
                                                         float(single_elem['geo']['lon']),
                                                         "Abstract"]
    print("Extract geo_data from titles")
    geo_list_title = []
    if not os.path.exists(path_ser_title):
        os.makedirs(path_ser_title)
    if len(os.listdir(path_ser_title)) == 0:
        for bact_tmp in bacts:
            count = 0
            name_tmp = bact_tmp.name.strip()
            if name_tmp[-1] == ';':
                name_tmp = name_tmp[:-1]
            for paper_tmp in bact_tmp.papers:
                tmp = GeoData(paper_tmp.pmid, paper_tmp.title,
                              geo.geoparse(paper_tmp.title), name_tmp)
                geo_list_title.append(tmp)
                with open(f"ser_title_{type_print}/{paper_tmp.pmid}_{count}.ser", "wb") as fw:
                    pickle.dump(tmp, fw)
    else:
        for filename in os.listdir(path_ser_title):
            with open(f"ser_title_{type_print}/{filename}", "rb") as f:
                geo_list_title.append(pickle.load(f))

    for geo_elem in geo_list_title:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            name_tmp = geo_elem.name.strip()
            if name_tmp[-1] == ';':
                name_tmp = name_tmp[:-1]
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    if not is_valid(single_elem, bact_names):
                        filtered_title.append(single_elem)
                    else:
                        geo_df.loc[len(geo_df.index)] = [id_tmp,
                                                         name_tmp,
                                                         float(single_elem['geo']['lat']),
                                                         float(single_elem['geo']['lon']),
                                                         "Title"]

    print("Extract geo data from descriptions")
    geo_list_bact = []
    if not os.path.exists(path_ser_bacts):
        os.makedirs(path_ser_bacts)
    if len(bacts) != len(os.listdir(path_ser_bacts)):
        for bact_tmp in bacts:
            tmp = GeoData(bact_tmp.id_bact, bact_tmp.description,
                          geo.geoparse(bact_tmp.description), bact_tmp.name)
            geo_list_bact.append(tmp)
            with open(f"ser_bacts_{type_print}/{bact_tmp.id_bact}.ser", "wb") as fw:
                pickle.dump(tmp, fw)
    else:
        for filename in os.listdir(path_ser_bacts):
            with open(f"ser_bacts_{type_print}/{filename}", "rb") as f:
                geo_list_bact.append(pickle.load(f))
    for geo_elem in geo_list_bact:
        if geo_elem.geo_dict:
            id_tmp = geo_elem.id_geo
            name_tmp = geo_elem.name.strip()
            if name_tmp[-1] == ';':
                name_tmp = name_tmp[:-1]
            for single_elem in geo_elem.geo_dict:
                if 'geo' in single_elem:
                    if not is_valid(single_elem, bact_names):
                        filtered_bact.append(single_elem)
                    else:
                        geo_df.loc[len(geo_df.index)] = [id_tmp,
                                                         name_tmp,
                                                         float(single_elem['geo']['lat']),
                                                         float(single_elem['geo']['lon']),
                                                         "Description"]

    return fix_lon_lat(geo_df)


def merge_geo_df(df1, df2, label1, label2):
    geo_df_total = pd.DataFrame(columns=["id", "name", "lat", "lon", "type", "label"])
    for index, row in df1.iterrows():
        geo_df_total.loc[len(geo_df_total.index)] = [row.id,
                                                     row.name,
                                                     row.lat,
                                                     row.lon,
                                                     row.type,
                                                     label1]
    for index, row in df2.iterrows():
        geo_df_total.loc[len(geo_df_total.index)] = [row.id,
                                                     row.name,
                                                     row.lat,
                                                     row.lon,
                                                     row.type,
                                                     label2]
    return geo_df_total


def display_map(geo_df, title, color, color_discrete_sequence):
    print(f"Producing map: {title}")
    fig = px.scatter_mapbox(geo_df,
                            hover_name="name",
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


def fix_lon_lat(df):
    df = df.drop_duplicates(subset=["id", "lat", "lon", "type"], keep="first")
    df = df.reset_index(drop=True)
    df = df.drop_duplicates(subset=["name", "lat", "lon", "type"], keep="first")
    df = df.reset_index(drop=True)
    lat_df = df.sort_values(by='lat')
    curr_lat = float(lat_df.loc[[0], "lat"])
    curr_incr = 0.5
    for idx, row in lat_df.iterrows():
        tmp_lat = float(row["lat"])
        if tmp_lat != curr_lat:
            curr_lat = float(row["lat"])
        else:
            lat_df.at[idx, "lat"] = float(row["lat"]) + curr_incr
            if idx != len(lat_df) - 1:
                if float(lat_df.at[idx, "lat"]) == float(lat_df.at[idx + 1, "lat"]):
                    lat_df.at[idx + 1, "lat"] = float(lat_df.at[idx + 1, "lat"]) + 0.5
                    curr_lat = float(lat_df.at[idx + 1, "lat"])
            curr_incr += 0.5
    lat_df = lat_df.sort_index()
    lon_df = lat_df.sort_values(by='lon')
    curr_lon = float(lon_df.loc[[0], "lon"])
    curr_incr = 0.5
    for idx, row in lon_df.iterrows():
        tmp_lon = float(row["lon"])
        if tmp_lon != curr_lon:
            curr_lon = float(row["lon"])
        else:
            lon_df.at[idx, "lon"] = float(row["lon"]) + curr_incr
            curr_incr += 0.5
            if idx != len(lon_df) - 1:
                if float(lon_df.at[idx, "lon"]) == float(lon_df.at[idx + 1, "lon"]):
                    lon_df.at[idx + 1, "lon"] = float(lon_df.at[idx + 1, "lon"]) + 0.5
                    curr_lon = float(lon_df.at[idx + 1, "lon"])
    final_df = lon_df.sort_index()
    return final_df

def geo_bar_chart(bacts_drug, bacts_nodrug, bacts_total):
    count_drug = len(bacts_drug)
    count_nodrug = len(bacts_nodrug)

    count_geo_drug = len(bacts_drug.index)
    count_geo_nodrug = len(bacts_nodrug.index)
    data = {
    'label' : ['With drug','Without drug'],
    'recorded':[count_drug,count_nodrug],
    'geolocalized' : [count_geo_drug,count_geo_nodrug]
    }

    width = 0.4
    fig, ax = plt.subplots()
    ax.bar(data['label'], data['recorded'], color='red', width=width, label='Recorded')
    ax.bar(data['label'], data['geolocalized'], color='green', width=width, label='Geolocalized')
    plt.title("Bacterial infections geolocalized using papers")
    plt.ylabel("Number of infections")
    ax.legend()
    fname = './plot_print/geolocalized_bar_plot.png'
    fig.set_size_inches((16, 12), forward=False)
    fig.savefig(fname, dpi=500)
