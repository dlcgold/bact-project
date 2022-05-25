# -*- coding: utf-8 -*-
"""progetto-ds.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1Vv2ixKWg8PMJHFMpj0APeFi33Z1z2AiH
"""

from Bio.KEGG import REST
import json
import os

# Estrazione batteri
sngOrg = REST.kegg_get("br:br08401").read()
H_list = []
check = True
bact = False
for line in sngOrg.splitlines():
  if not check:
    break
  if line == "ABacterial infections":
    bact = True
    continue
  if bact:
    if line[0] == 'C':
      #print(line)
      H_list.append(line.split()[1])
  if bact and line[0] == 'A':
    check = False

#sngOrg = REST.kegg_list("ds").read()
#sngOrg = REST.kegg_get(["H02511", "H02512"]).read()
#sngOrg = REST.kegg_get("br:br08401").read()
#print(sngOrg)
#json_data = json.loads(sngOrg)
#H_list = []
#print(json_data['children'][0]['children'][0]['children'][0]['name'])
#print(len(json_data['children']))
#for line in sngOrg.splitlines():
#  if line[0] == 'C':
#    H_list.append(line.split()[1])
# H_list = []
#for line in sngOrg.splitlines():
#  code = line.split("\t")[0]
#  code_split = code.split(":")
#  if code_split[0] == "ds":
#    H_list.append(code_split[1])
print(H_list)
print(len(H_list))
no_paper = 0
drug_list = []
no_drug_list = []
all_get = ""
for h_id in H_list:
  tmp = REST.kegg_get(h_id).read()
  with open(f"kegg_get/{h_id}.txt", "w") as f:
    f.write(tmp)
  all_get = all_get + tmp
  all_get = all_get + "\n"
  if 'DRUG' in tmp:
    drug_list.append(h_id)
  else:
    no_drug_list.append(h_id)
  if not 'PMID' in tmp:
    #print(tmp)
    no_paper = no_paper + 1
print(drug_list)
print(no_drug_list)
#sngOrg = REST.kegg_get(["H00349"]).read()
#print('DRUG' in sngOrg)

print(len(drug_list))
print(len(no_drug_list))
print(no_paper)

#sngOrg = REST.kegg_get("gn:ype").read()

sngOrg = REST.kegg_get("ds:H00111").read()

print(sngOrg)

class paper:
   def __init__(self, pmid, authors, title, journal, doi):
    self.pmid = pmid
    self.authors = authors
    self.title = title
    self.journal = journal
    self.doi = doi

   def __repr__(self):
    return f"{self.title} ({self.pmid})\t{self.authors}\t{self.journal} ({self.doi})"
      
class drug:
  def __init__(self, name, id):
    self.name = name
    self.id = id

  def __repr__(self):
    return f"{self.name}, ({self.id})"

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
      bact_name = "".join(spl[i]+ " "  for i in range(1, len(spl)))
    if spl[0] == "ENTRY":
      bact_id = "".join(spl[i]+ " "  for i in range(1, len(spl))).split()[0]
    if spl[0] == "CATEGORY":
      bact_cat = "".join(spl[i]+ " "  for i in range(1, len(spl)))
    if spl[0] == "DESCRIPTION":
      bact_des = "".join(spl[i]+ " "  for i in range(1, len(spl)))

    if sub_bool == False and spl[0] == "BRITE":
      sub_bool = True
      sub_count += 1
    elif sub_bool == True and sub_count == 1:
      sub_count += 1
    elif sub_bool == True and sub_count == 2:
      sub_bool = False
      bact_sub = "".join(spl[i]+ " "  for i in range(0, len(spl)))

    if drug_bool == False and spl[0] == "DRUG":
      name = "".join(spl[i]+ " "  for i in range(1, len(spl)-1))
      drugs.append(drug(name, spl[-1].replace("[", "").replace("]", "")))
      drug_bool = True
    elif drug_bool == True and spl[0].upper() != spl[0]:
      name = "".join(spl[i]+ " "  for i in range(0, len(spl)))
      drugs.append(drug(name, spl[-1].replace("[", "").replace("]", "")))
    elif drug_bool == True and spl[0].upper() == spl[0]:
      drug_bool = False

    if paper_bool == False and spl[0] == "REFERENCE":
      if len(spl) > 1 and spl[1].split(":")[0] == "PMID":
        id_tmp = spl[1].split(":")[1].replace('\n',' ')
      paper_bool = True
    elif paper_bool == True and spl[0] == "AUTHORS":
      authors_tmp = "".join(spl[i]+ " "  for i in range(1, len(spl)))
    elif paper_bool == True and spl[0] == "TITLE":
      title_tmp = "".join(spl[i]+ " "  for i in range(1, len(spl)))
    elif paper_bool == True and spl[0] == "JOURNAL":
      journal_tmp = "".join(spl[i]+ " "  for i in range(1, len(spl)))
      doi_bool = True
    elif paper_bool == True and doi_bool == True and spl[0].upper() != spl[0]:
      doi_tmp = spl[0].split(":")[1].replace('\n',' ')
      doi_bool = False
      paper_bool = False
      papers.append(paper(id_tmp, authors_tmp, title_tmp, journal_tmp, doi_tmp))

  #print(drugs)
  #print(papers)

  bact_tmp = bact(bact_name, "ds:" + bact_id, bact_des, bact_cat, bact_sub, drugs, papers)
  #print(bact_tmp)
  return bact_tmp

bact_tmp = parse(sngOrg)

#print(bact)

bacts = []
for filename in os.listdir("kegg_get"):
  print(filename)
  with open("kegg_get/"+filename, "r") as f:
    bacts.append(parse(str(f.read())))

print(bacts[0])