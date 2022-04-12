import hashlib
import os
import requests
from typing import List, Dict

from biothings.utils.dataload import tabfile_feeder, dict_sweep, unlist

"""
  | Column Names             | Key Names
--+--------------------------+---------------------------
0 | gene_name                | object.SYMBOL
1 | gene_claim_name          |
2 | entrez_id                | object.NCBIGene
3 | interaction_claim_source | association.provided_by
4 | interaction_types        | association.relation_name
5 | drug_claim_name          |
6 | drug_claim_primary_name  |
7 | drug_name                | subject.name
8 | drug_concept_id          | subject.CHEMBL_COMPOUND
9 | interaction_group_score  | association.interaction_group_score
10| PMIDs                    | association.pubmed
"""


def query_entrez_id(gene_name) -> str:
    """
    Find entrez id associated with given gene name through MyGene API.
    """
    query = "http://mygene.info/v3/query?q=symbol:{}&fields=entrezgene&species=human".format(gene_name)
    response = requests.get(query)
    if response.status_code == 200:
        data = response.json()
        if data['hits']:  # data['hits'] could be an empty list
            entrez_id = data['hits'][0]['entrezgene']
            if entrez_id != "":
                return entrez_id
    return None


def query_chembl_id(drug_name) -> str:
    """
    Find chembl id associated with given drug name through MyChem API.
    """
    query = "http://mychem.info/v1/query?q=chembl.pref_name:{}&fields=chembl.molecule_chembl_id".format(drug_name)
    response = requests.get(query)
    if response.status_code == 200:
        data = response.json()
        if data['hits']:  # data['hits'] could be an empty list
            chembl_id = data['hits'][0]['chembl']['molecule_chembl_id']
            if chembl_id != "":
                return chembl_id
    return None


def create_doc_id(record: List[str]) -> str:
    """
    For each record (received as a list of strings), create a hash string as its ID
    """
    bytestr = bytearray("-".join(record), 'utf-8')
    hashstr = hashlib.blake2b(bytestr, digest_size=8).hexdigest()
    return hashstr


def create_column_index(header: List[str]) -> Dict[str, int]:
    """
    Convert a list of column names into a column-name-to-index dict.

    E.g. with header=['gene_name', 'gene_claim_name', 'entrez_id'], the output index is
    {'gene_name': 0, 'gene_claim_name': 1, 'entrez_id': 2}.
    """
    return {header[i]: i for i in range(len(header))}


def is_empty(val: str) -> bool:
    """
    Check if a string value is empty.
    """
    return (val is None) or (val == "")


def verify_entrez_id(entrez_id: str, gene_name: str) -> str:
    """
    Check if the given entrez id is emtpy. If yes, find the entrez id by querying with the gene name;
    otherwise, return the entrez id as-is.
    """
    if is_empty(entrez_id):
        if is_empty(gene_name):
            return None
        else:
            entrez_id = query_entrez_id(gene_name)

    # return as-is if entrez_id is not empty
    return entrez_id


def create_object_id(entrez_id: str, gene_name: str) -> str:
    """
    Try creating an object id first from the entrez id.
    If the entrez id is empty, create with the gene name.
    If both the entrez id and the gene name are empty, return None
    """
    if is_empty(entrez_id):
        if is_empty(gene_name):
            return None  # a None id indicates discarding the whole record
        else:
            return 'name:' + gene_name

    return 'NCBIGene:' + entrez_id


def verify_drug_concept_id(drug_concept_id: str, drug_name: str) -> str:
    """
    Check if the given drug concept id contains a valid chembl id.

    A given drug concept id could be
        
        1. emtpy,
        2. a CURIE-like wikidata id, e.g. "wikidata:Q419808", or
        3. a CURIE-like chembl id, e.g. chembl:CHEMBL942.

    In case 1 and 2, the drug name is used to query its associated chemble id, without the "chembl:" prefix.
    If the drug name is empty, return None.
    
    In case 3, a chembl id is returned from substringing the given drug concept id.
    """
    if is_empty(drug_concept_id) or drug_concept_id.startswith("wikidata:"):
        if is_empty(drug_name):
            return None
        drug_chembl_id = query_chembl_id(drug_name)
    elif drug_concept_id.startswith("chembl:"):
        drug_chembl_id = drug_concept_id.split(':')[-1]
    else:
        raise ValueError(f"Cannot parse drug concept id {drug_concept_id} (associated drug name is {drug_name})")

    return drug_chembl_id


def create_subject_id(drug_chembl_id: str, drug_name: str) -> str:
    """
    Try creating a subject id first from the chembl id.
    If the chembl id is empty, create with the drug name.
    If both the chembl id and the drug name are empty, return None
    """
    if is_empty(drug_chembl_id):
        if is_empty(drug_name):
            return None  # a None id indicates discarding the whole record
        else:
            return "name:" + drug_name

    return 'CHEMBL.COMPOUND:' + drug_chembl_id


def load_annotations(data_folder):
    data_file = os.path.join(data_folder, "interactions.tsv")
    data = tabfile_feeder(data_file, header=0)
    header = next(data)
    col_index = create_column_index(header)

    for rec in data:
        # Document framework
        doc = {
            "_id": create_doc_id(rec),
            "subject": {},
            "object": {},
            "association": {}
        }

        # Object
        entrez_id = rec[col_index["entrez_id"]]
        gene_name = rec[col_index["gene_name"]]
        entrez_id = verify_entrez_id(entrez_id=entrez_id, gene_name=gene_name)
        object_id = create_object_id(entrez_id=entrez_id, gene_name=gene_name)
        if object_id is None:
            continue
        else:
            doc['object']['NCBIGene'] = entrez_id
            doc['object']['SYMBOL'] = gene_name
            doc['object']['id'] = object_id

        # Subject
        drug_name = rec[col_index["drug_name"]]
        drug_concept_id = rec[col_index["drug_concept_id"]]
        drug_chembl_id = verify_drug_concept_id(drug_concept_id=drug_concept_id, drug_name=drug_name)
        subject_id = create_subject_id(drug_chembl_id=drug_chembl_id, drug_name=drug_name)
        if subject_id is None:
            continue
        else:
            doc['subject']['name'] = drug_name
            doc['subject']['CHEMBL_COMPOUND'] = drug_chembl_id
            doc['subject']['id'] = subject_id

        # Association
        interaction_types = rec[col_index["interaction_types"]].replace(" ", "_").split(",")
        interaction_claim_source = rec[col_index["interaction_claim_source"]]
        interaction_group_score = rec[col_index["interaction_group_score"]]
        pmids = rec[col_index["PMIDs"]].split(",")

        doc['association']['relation_name'] = interaction_types
        doc['association']['provided_by'] = interaction_claim_source
        doc['association']['interaction_group_score'] = float(interaction_group_score)
        doc['association']['pubmed'] = pmids
        
        # Cleanup
        doc = dict_sweep(doc)
        doc = unlist(doc)
        yield doc
