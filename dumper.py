"""
Dumper class for handling the DGIdb data download
"""

from pathlib import Path

from biothings import config
import biothings.hub


class DGIdbDumper(biothings.hub.dataload.dumper.LastModifiedHTTPDumper):
    DGIDB_URL = "https://www.dgidb.org/"
    SCHEDULE = "30 0 * * 0"
    SRC_NAME = "DGIdb"
    SRC_ROOT_FOLDER = Path(config.DATA_ARCHIVE_ROOT) / SRC_NAME
    SRC_URLS = ["https://www.dgidb.org/data/2023-Dec/interactions.tsv"]
    UNCOMPRESS = False

    def get_release(self) -> str:
        """
        DGIdb (v5.0.7 - sha 35436d7)
        https://github.com/dgidb/dgidb-v5/releases/tag/v.5.0.7
        Last updated 2024-06-07
        """
        return "5.0.7"
