"""
Definitions for the structure of the interactions.tsv file

Data Structure 2022-Feb
   | Column Names             | Key Names
---+--------------------------+---------------------------
0  | gene_name                | object.SYMBOL
1  | gene_claim_name          |
2  | entrez_id                | object.NCBIGene
3  | interaction_claim_source | association.interaction_claim_source
4  | interaction_types        | association.interaction_types
5  | drug_claim_name          |
6  | drug_claim_primary_name  |
7  | drug_name                | subject.drug_name
8  | drug_concept_id          | subject.CHEMBL_COMPOUND
9  | interaction_group_score  | association.interaction_group_score
10 | PMIDs                    | association.pmids


Data Structure 2023-Dec
   | Column Names                  | Key Names
---+------------------------------+---------------------------
0  | gene_claim_name               |
1  | gene_concept_id               | object.id
2  | gene_name                     | object.name
3  | interaction_source_db_name    | interaction.source_db.name
4  | interaction_source_db_version | interaction.source_db.version
5  | interaction_type              | interaction.type
6  | interaction_score             | interaction.score
7  | drug_claim_name               |
8  | drug_concept_id               | subject.id
9  | drug_name                     | subject.name
10 | approved                      | subject.approved
11 | immunotherapy                 | subject.immunotherapy
12 | anti_neoplastic               | subject.anti_neoplastic
"""

import dataclasses


@dataclasses.dataclass(frozen=True)
class InteractionsFileEntry:
    """
    ==> interactions.tsv <==
    [00] gene_claim_name	          | NRAS
    [01] gene_concept_id	          | hgnc:7989
    [02] gene_name                    | NRAS
    [03] interaction_source_db_name	  | JAX-CKB
    [04] interaction_source_db_version| 5-Jun-23
    [05] interaction_type             | NULL
    [06] interaction_score            | 0.065654977
    [07] drug_claim_name	          | Palbociclib
    [08] drug_concept_id	          | rxcui:1601374
    [09] drug_name	                  | Palbociclib
    [10] approved                     | TRUE
    [11] immunotherapy                | FALSE
    [12] anti_neoplastic              | TRUE
    """

    gene_claim_name: str
    gene_concept_id: str
    gene_name: str
    interaction_source_db_name: str
    interaction_source_db_version: str
    interaction_type: str
    interaction_score: float
    drug_claim_name: str
    drug_concept_id: str
    drug_name: str
    approved: bool
    immunotherapy: bool
    anti_neoplastic: bool
