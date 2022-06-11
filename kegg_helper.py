from Bio import Entrez
from Bio.KEGG import REST

Entrez.email = "d.cozzi@campus.unimib.it"


def get_pathway_target_from_kegg(drug_id):
    data = REST.kegg_get(drug_id).read()
    pathways = []
    path_bool = False
    for line in data.splitlines():
        spl = line.strip().split()
        if spl[0] == "PATHWAY":
            name = "".join(spl[i] + " " for i in range(2, len(spl)))
            pathways.append((spl[1].split("(")[0], name))
            path_bool = True
        elif spl[0].upper() != spl[0] and path_bool:
            name = "".join(spl[i] + " " for i in range(1, len(spl)))
            pathways.append((spl[0].split("(")[0], name))
            path_bool = True
        elif spl[0].upper() == spl[0] and path_bool:
            break
    return pathways


tmp = get_pathway_target_from_kegg("D11713")
print(tmp)
