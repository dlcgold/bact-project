from bact_classes import *
import requests


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
    pathogens = []
    pathogens_bool = False
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

        if not pathogens_bool and spl[0] == "PATHOGEN":
            name = "".join(spl[i] + " " for i in range(1, len(spl)))
            name = name.strip().split('[')[0].strip()
            pathogens.append(name)
            pathogens_bool = True
        elif pathogens_bool and spl[0].upper() != spl[0]:
            name = "".join(spl[i] + " " for i in range(0, len(spl)))
            name = name.strip().split('[')[0].strip()
            pathogens.append(name)
        elif pathogens_bool and spl[0].upper() == spl[0]:
            pathogens_bool = False

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
                    pathways, pathogens)
    return bact_tmp


def list_diff(li1, li2):
    return list(set(li1) - set(li2)) + list(set(li2) - set(li1))


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
