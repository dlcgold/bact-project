## lista di disease con drug
## lista di disease senza drug


# entrez


#1. kegg genome con codice disease
#esempio https://www.genome.jp/dbget-bin/get_linkdb?-t+genome+ds:H00113


#2. kegg Genome entry
    #prendo i bioProject (ho un submission e i biosample)

    # prendo assembly, ho biosample e bioproject, e genome
#3. prendo l'assembly con le statistics (non esiste su ENTREZ, devo passare da Genome)
    # BIOsample con Geotagging (per i bacterical)  (posso arrivarci direttamente da Entrez)
    # GENOME
        # gc content
        # protein count
        # protein list 
            # protein sample con geotagging (ridondante con biosample)

        

# biosamples

from Bio import Entrez
from Bio.KEGG import REST
from numpy import result_type
# from drugbank_sql_lite import *
# from bact_classes import * 
import requests as req
from bs4 import BeautifulSoup as bfs
from genomejp import *


def main():
    genomes = get_genome_id_from_disease_id("H01451")
    assemblies = []
    biosamples = []
    for ids in genomes:
        ass = get_assembly_for_id(ids)
        bios = get_biosample_for_id(ids)
        assemblies.append(ass)
        biosamples.append(bios)
    print(biosamples)
    print(assemblies)
    print(genomes)

    Entrez.email = "m.sgro2@campus.unimib.it"
    # assembly_entry = Entrez.efetch(db="all", id=assemblies[0])
    # print(assembly_entry.read())
    # assembly_entry = Entrez.esearch(term=assemblies[0], db="pubmed", retmax=1)
    # print(assembly_entry.read())
    bio_proj = Entrez.esearch(term="GCA_000092985.1", db="genome", rettype="gb", rettmode='text')
    genome_id = Entrez.read(bio_proj)["IdList"][0]
    genome_ret = Entrez.efetch(db="genome", id=genome_id, retmode="xml")
    print(Entrez.read(genome_ret))

    bio2 = Entrez.efetch(id=biosamples[0], db="bioproject", rettmode='text')
    # print(Entrez.read(bio2))
    print(bio2.read())
    bio3 = Entrez.efetch(id="SAMN02603853", db="biosample", rettmode='text')
    print(bio3.read())
if __name__ == "__main__":
    main()