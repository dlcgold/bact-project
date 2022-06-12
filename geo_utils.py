from mordecai import Geoparser
import requests as req
from bs4 import BeautifulSoup as bfs


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

    response = req.get(url="https://en.wikipedia.org/wiki/" + word.replace(" ", "_"))
    if response.status_code == 200:
        is_species = response.text.find("Scientific classification") != -1

    valid = only_letters and not upper_case and not no_cap_letters and not no_country and not is_bacteria and not is_species

    return valid