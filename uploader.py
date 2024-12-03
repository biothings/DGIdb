"""
Uploader for the FDA Drugs data plugin

Handles the parsing of the document generated post_dump
from the dumper class instance
"""

from pathlib import Path
from typing import Dict, List, Set, Tuple, Union
import csv
import logging

import httpx

import biothings.hub
from biothings import config
from biothings.hub.dataload.uploader import BaseSourceUploader
from biothings.utils.dataload import dict_sweep, unlist
from biothings.utils.storage import BasicStorage

# from biothings import config
from biothings.utils.common import iter_n
from .file_definitions import InteractionsFileEntry

logger = config.logger


class DGIdbUploader(BaseSourceUploader):

    __database__ = config.DATA_SRC_DATABASE

    storage_class = BasicStorage

    # name of the resource and collection name used to store data
    # (see regex_name though for exceptions)
    name = "DGIdb"
    # if several resources, this one if the main name,
    # it's also the _id of the resource in src_dump collection
    # if set to None, it will be set to the value of variable "name"
    main_source = None
    # in case resource used split collections (so data is spread accross
    # different colleciton, regex_name should be specified so all those split
    # collections can be found using it (used when selecting mappers for instance)
    regex_name = None

    # number of archived collection to keep. Oldest get dropped first.
    keep_archive = 10

    def __init__(
        self,
        db_conn_info: Tuple[str, str],
        collection_name: str = None,
        log_folder: Union[str, Path] = None,
        *args,
        **kwargs,
    ):
        super().__init__(
            db_conn_info,
            collection_name,
            log_folder,
            *args,
            **kwargs,
        )

    def batch_query_entrez_id(
        self, http_client: httpx.Client, gene_concept_ids: List[str]
    ):
        """
        Return the batch query result.
        This is a wrapper for POST query of "/query" service.

        :return: a list of matching objects
        """
        query_scopes = set()
        query_terms = []
        for concept_id in gene_concept_ids:
            if concept_id is not None and concept_id != "NULL":
                scope, _, identifier = concept_id.partition(":")
                query_terms.append(identifier)
                query_scopes.add(scope)

        query_scopes = ",".join([f"{scope}" for scope in query_scopes])

        for query_batch_terms in iter_n(query_terms, 1000, with_cnt=False):
            formatted_batch_terms = ",".join([f"{term}" for term in query_batch_terms])
            query_data = {
                "q": formatted_batch_terms,
                "scopes": query_scopes,
                "fields": ["_id", "entrezgene"],
            }
            endpoint = "https://mygene.info/v3/query"
            response = http_client.post(url=endpoint, data=query_data)
            yield from response.json()

    def batch_query_chembl_id(self, http_client: httpx.Client, drug_names: List[str]):
        """
        Return the batch query result.
        This is a wrapper for POST query of "/query" service.

        :return: a list of matching objects
        """
        for batch_drug_names in iter_n(drug_names, 1000, with_cnt=False):
            query_data = {
                "q": batch_drug_names,
                "scopes": "chembl.pref_name",
                "fields": ["chembl.molecule_chembl_id"],
            }
            endpoint = "https://mychem.info/v1/query"
            response = http_client.post(url=endpoint, json=query_data, timeout=15)
            yield from response.json()

    def _build_gene_lookup_table(
        self, http_client: httpx.Client, gene_concept_ids: List[str]
    ) -> Dict[str, str]:
        """
        Builds the entrezgene lookup table.

        Takes the gene concept ids from the interactions table
        and creates a mapping. The mapping relationship is
        gene_concept_id -> entrezgene

        If no valid entrezgene is found for the concept id we exclude it from the lookup table
        """
        gene_lookup_table = {}
        for response in self.batch_query_entrez_id(http_client, gene_concept_ids):
            entrez_id = response.get("entrezgene", None)
            if entrez_id != "" and entrez_id is not None:
                try:
                    gene_concept_id = gene_concept_ids.pop()
                except IndexError:
                    break
                gene_lookup_table[gene_concept_id] = entrez_id
        return gene_lookup_table

    def _build_drug_lookup_table(
        self, http_client: httpx.Client, drug_names: List[str]
    ) -> Dict[str, str]:
        """
        Builds the chembl id lookup table.

        Takes the drug names from the interactions table
        and creates a mapping. The mapping relationship is
        drug_name -> chembl molecular id

        If no valid chembl molecular id is found for the concept id
        we exclude it from the lookup table
        """
        drug_lookup_table = {}
        for response in self.batch_query_chembl_id(http_client, drug_names):
            chembl_structure = response.get("chembl", None)
            if chembl_structure is not None:
                molecule_chembl_id = chembl_structure.get("molecule_chembl_id", None)
                if molecule_chembl_id is not None:
                    try:
                        drug_name = drug_names.pop()
                    except IndexError:
                        break
                    drug_lookup_table[drug_name] = molecule_chembl_id
        return drug_lookup_table

    def _build_lookup_tables(
        self, interactions_file: Union[str, Path]
    ) -> Tuple[Dict, Dict]:
        """
        Builds the lookup tables for both the chembl id and entrez id

        Batch performs all queries at once to minimize HTTP requests
        """
        gene_concept_ids = set()
        drug_names = set()
        with open(str(interactions_file), "r", encoding="utf-8") as file_handle:
            reader = csv.DictReader(file_handle, delimiter="\t")
            for entry in reader:
                gene_concept_ids.add(entry["gene_concept_id"])
                drug_names.add(entry["drug_name"])

        gene_concept_ids = list(gene_concept_ids)
        drug_names = list(drug_names)

        with httpx.Client() as http_client:
            entrez_id_lookup_table = self._build_gene_lookup_table(
                http_client, gene_concept_ids
            )
            chembl_id_lookup_table = self._build_drug_lookup_table(
                http_client, drug_names
            )
            return (chembl_id_lookup_table, entrez_id_lookup_table)

    def load_data(self, data_path: Union[str, Path]):
        data_path = Path(data_path).resolve().absolute()
        interactions_file = data_path.joinpath("interactions.tsv")
        chembl_lookup_table, entrez_lookup_table = self._build_lookup_tables(
            interactions_file
        )

        with open(str(interactions_file), "r", encoding="utf-8") as file_handle:
            reader = csv.reader(file_handle, delimiter="\t")
            header = next(reader)
            logger.debug("Interactions headers: [%s]", header)

            for row_index, entry in enumerate(reader):
                interaction = InteractionsFileEntry(*entry)

                if interaction.drug_name is None or interaction.drug_name == "NULL":
                    logger.warning(
                        "Invalid drug name %s from entry. Skipping entry #%s",
                        interaction.drug_name,
                        row_index + 1,
                    )
                    continue

                if interaction.gene_name is None or interaction.gene_name == "NULL":
                    logger.warning(
                        "Invalid gene name %s from entry. Skipping entry #%s",
                        interaction.gene_name,
                        row_index + 1,
                    )
                    continue

                chembl_id = chembl_lookup_table.get(interaction.drug_name, None)
                entrezgene = entrez_lookup_table.get(interaction.gene_concept_id, None)

                if entrezgene is None:
                    logger.warning(
                        "Unable to lookup the entrez id from entry. Skipping entry #%s",
                        row_index + 1,
                    )
                    continue

                if chembl_id is None:
                    logger.warning(
                        "Unable to lookup the chembl id from entry. Skipping entry #%s",
                        row_index + 1,
                    )
                    continue

                document_identifier = f"{interaction.drug_name}_{interaction.interaction_type}_{interaction.gene_name}"
                subject_structure = {
                    "id": f"CHEMBL.COMPOUND:{chembl_id}",
                    "name": interaction.drug_name,
                    "type": "Drug",
                    "chembl": chembl_id,
                    "approved": interaction.approved,
                    "immunotherapy": interaction.immunotherapy,
                    "anti_neoplastic": interaction.anti_neoplastic,
                }

                object_structure = (
                    {
                        "id": f"NCBIGene:{entrezgene}",
                        "name": interaction.gene_name,
                        "type": "Gene",
                        "hgnc": entrezgene,
                    },
                )

                document = {
                    "_id": document_identifier,
                    "subject": subject_structure,
                    "object": object_structure,
                    "interaction": {
                        "source_db": {
                            "name": interaction.interaction_source_db_name,
                            "version": interaction.interaction_source_db_version,
                        },
                        "type": interaction.interaction_type,
                        "score": interaction.interaction_score,
                    },
                }
                document = unlist(dict_sweep(document))
                yield document
