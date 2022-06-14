"""
progetto-ds.ipynb
"""

from copy import deepcopy

from Bio.KEGG import REST

from drugbank_sql_lite import *
from geo_parsing import *

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
    def __init__(self, name, id_bact, description, category, sub, drugs, papers, pathways):
        self.name = name
        self.id_bact = id_bact
        self.drugs = drugs
        self.pathways = pathways
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


class Pathway:
    def __init__(self, name, id_path):
        self.name = name
        self.id_path = id_path

    def __repr__(self):
        return f"{self.name}, ({self.id_path})"


class GeoData:
    def __init__(self, id_geo, text, geo_dict):
        self.id_geo = id_geo
        self.text = text
        self.geo_dict = geo_dict


class DrugBankObj():
    def __init__(self, drugname):
        self.drugname = drugname
        self.query_url = 'https://go.drugbank.com/unearth/q?searcher=drugs&query=' + self.drugname
        try:
            self.url = get_drugbank_link(retrive_page(self.query_url))
        except:
            self.url = 'NoURL'

        try:
            self.index = str(get_db_index(retrive_page(self.url)))
        except:
            self.index = 'NoID'

    def get_page(self):
        try:
            return retrive_page(self.url)
        except:
            return 'NoPage'


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
    pathways = []
    pathways_bool = False
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
            name = name.strip()
            drugs.append(Drug(name, spl[-1].replace("[", "").replace("]", "").split(':')[-1]))
            drug_bool = True
        elif drug_bool and spl[0].upper() != spl[0]:
            name = "".join(spl[i] + " " for i in range(0, len(spl) - 1))
            name = name.strip()
            drugs.append(Drug(name, spl[-1].replace("[", "").replace("]", "").split(':')[-1]))
        elif drug_bool and spl[0].upper() == spl[0]:
            drug_bool = False

        if not drug_bool and spl[0] == "PATHWAY":
            name = "".join(spl[i] + " " for i in range(2, len(spl)))
            name = name.strip()
            pathways.append(Pathway(name, spl[1]))
            pathways_bool = True
        elif pathways_bool and spl[0].upper() != spl[0]:
            name = "".join(spl[i] + " " for i in range(1, len(spl)))
            name = name.strip()
            pathways.append(Pathway(name, spl[0]))
        elif pathways_bool and spl[0].upper() == spl[0]:
            pathways_bool = False

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

    bact_tmp = Bact(bact_name, "ds:" + bact_id, bact_des, bact_cat, bact_sub, drugs, papers,
                    pathways)
    return bact_tmp


def list_diff(li1, li2):
    return list(set(li1) - set(li2)) + list(set(li2) - set(li1))


def get_drug_kegg(drug_id):
    sngOrg = REST.kegg_get([drug_id]).read()
    drug_name = ""
    for line in sngOrg.splitlines():
        spl = line.split()
        if spl[0] == "NAME":
            drug_name = "".join(spl[i] + " " for i in range(1, len(spl))).strip()
            if drug_name[-1] == ";":
                drug_name = drug_name[:len(drug_name) - 1]
            break
    return Drug(drug_name, drug_id)


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


def retrive_page(url):
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0'
    }
    response = req.get(url, headers)
    soup = bfs(response.content, 'html.parser')
    return soup


def get_drugbank_link(soup):
    links = []
    for link in soup.find_all('link'):
        if 'drugbank' in link.get('href'):
            links.append(link.get('href'))
    if len(links) == 0:
        return ""
    return links[0]


def get_db_index(soup):
    for link in soup.find_all('link'):
        if 'drugbank' in link.get('href'):
            return link.get('href').split('/')[-1]


# kegg_drugs GET
def get_kegg_drugs(query):
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


def main():
    # Estrazione batteri
    print("getting bacterial infections list")
    type_infection = "Bacterial infections"
    bacterial_infections = REST.kegg_get("br:br08401").read()
    H_list = []
    check = True
    bact_bool = False

    # geo locations filtered
    filtered_bact = []
    filtered_abs = []
    filtered_title = []
    filtered_bact_drug = []
    filtered_abs_drug = []
    filtered_title_drug = []

    for line in bacterial_infections.splitlines():
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
    print(f"Total number of bacteria: {len(H_list)}")
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

    print(drug_list)
    print(no_drug_list)

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
            """tmp_extra_drugs = get_kegg_drugs(tmp_bact.name)
            if len(tmp_extra_drugs) > 0:
                tmp_drugs = []
                for elem in tmp_bact.drugs:
                    tmp_drugs.append(elem.id_drug)
                # print(tmp_drugs)
                diffs = list_diff(tmp_extra_drugs, tmp_drugs)
                print("add extra drugs")
                if len(diffs) > 0:
                    for elem in diffs:
                        if elem not in tmp_drugs:
                            tmp_bact.drugs.append(get_drug_kegg(elem))"""
            bacts.append(tmp_bact)

    # print(f"Total: {len(bacts)} bacteria")

    bacts_nodrug = deepcopy(bacts)

    print("parse KEGG data with drugs")
    bacts_drug = []
    for filename in os.listdir("kegg_get_drug"):
        with open("kegg_get_drug/" + filename, "r") as f:
            tmp_bact = parse(str(f.read()))
            bacts_drug.append(tmp_bact)
            tmp_extra_drugs = get_kegg_drugs(tmp_bact.name)
            if len(tmp_extra_drugs) > 0:
                tmp_drugs = []
                for elem in tmp_bact.drugs:
                    tmp_drugs.append(elem.id_drug)
                # print(tmp_drugs)
                diffs = list_diff(tmp_extra_drugs, tmp_drugs)
                print("add extra drugs")
                if len(diffs) > 0:
                    for elem in diffs:
                        if elem not in tmp_drugs:
                            tmp_bact.drugs.append(get_drug_kegg(elem))
            bacts.append(tmp_bact)

    print(f"Total: {len(bacts_drug)} bacteria with drugs")
    print(f"Total: {len(bacts_nodrug)} bacteria without drugs")

    # Geoparse of bacteria without associated drugs
    geo_df_nodrug = geo_parse(bacts_nodrug)

    # Geoparse of bacteria without associated drugs
    geo_df_drug = geo_parse(bacts_drug)

    geo_df_total = merge_geo_df(geo_df_nodrug, geo_df_drug)

    """
    # Dump data
    if not os.path.exists('dumps'):
        os.makedirs("dumps")

    with open(f"dumps/bacts.ser", "wb") as fw:
        pickle.dump(bacts, fw)
    with open(f"dumps/bacts_drug.ser", "wb") as fw:
        pickle.dump(bacts_drug, fw)
    with open(f"dumps/bacts_nodrug.ser", "wb") as fw:
        pickle.dump(bacts_nodrug, fw)
    with open(f"dumps/geo_df_nodrug.ser", "wb") as fw:
        pickle.dump(geo_df_nodrug, fw)
    with open(f"dumps/geo_df_drug.ser", "wb") as fw:
        pickle.dump(geo_df_drug, fw)
    with open(f"dumps/geo_df_total.ser", "wb") as fw:
        pickle.dump(geo_df_total, fw)
    with open(f"dumps/geo_list_nodrug.ser", "wb") as fw:
        pickle.dump(geo_list_nodrug, fw)
    with open(f"dumps/geo_list_drug.ser", "wb") as fw:
        pickle.dump(geo_list_drug, fw)
    """

    drugs_map = display_map(geo_df_drug, "Locations associated to bacteria with drugs", "type",
                            ["red", "blue", "green"])
    nodrugs_map = display_map(geo_df_nodrug, "Locations associated to bacteria without drugs",
                              "type", ["red", "blue", "green"])
    all_map = display_map(geo_df_total, "Locations associated to all bacteria", "drug",
                          ["red", "green"])

    #
    # if not os.path.exists('ser_drug_bank') or not os.listdir("ser_drug_bank"):
    #     db_count = 0
    #     drug_count = 0
    #     if not os.path.exists('ser_drug_bank'):
    #         os.makedirs('ser_drug_bank')
    #     no_id_count = 0
    #     for bact_tmp in bacts:
    #         for drug_tmp in bact_tmp.drugs:
    #             drug_count += 1
    #             tmp = DrugBankObj(drug_tmp.name)
    #             if tmp.index != "NoID":
    #                 db_count += 1
    #                 with open(f"ser_drug_bank/{tmp.index}.ser", "wb") as fw:
    #                     pickle.dump(tmp, fw)
    #             else:
    #                 with open(f"ser_drug_bank/{tmp.index}_{no_id_count}.ser", "wb") as fw:
    #                     pickle.dump(tmp, fw)
    #                     no_id_count += 1
    #
    #     print(f"{db_count} vs {drug_count} | lost {drug_count - db_count} drugs")

    count_target = 0
    count_all = 0
    for bact_tmp in bacts:
        name = bact_tmp.name.replace("infection", "").strip()
        if name[-1] == ";":
            name = name[:len(name) - 1]

        drugs_db = get_drugs_for_target(name)
        if len(drugs_db) != 0:
            count_target += 1
            print("YES TARGET")
        else:
            print(name)
        drugs_db = get_drugs_for_all(name)
        if len(drugs_db) != 0:
            count_all += 1
            print("YES ALL")
        else:
            print(name)
    print(f"for {len(bacts)} bacts: {count_target} target and {count_all} all")

    print(f"bacts: {len(bacts)}")
    print(f"bacts_drug: {len(bacts_drug)}")


if __name__ == "__main__":
    main()
