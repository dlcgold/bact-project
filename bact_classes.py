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
    def __init__(self, name, id_bact, description, category, sub, drugs, papers, pathways, pathogens):
        self.name = name
        self.id_bact = id_bact
        self.drugs = drugs
        self.pathways = pathways
        self.papers = papers
        self.description = description
        self.category = category
        self.sub = sub
        self.pathogens = pathogens

    def __repr__(self):
        return f"{self.name} ({self.id_bact})\n{self.drugs}\n ({self.papers})"


class Drug:
    def __init__(self, name, id_drug, origin=""):
        self.name = name
        self.id_drug = id_drug
        self.origin = origin

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