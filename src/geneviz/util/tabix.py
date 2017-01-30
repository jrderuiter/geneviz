"""Functionality for creating and dealing with tabix-indexed files."""

import contextlib
import itertools
from pathlib import Path
import subprocess
from typing import Any, Callable, Iterable, List

import pysam

GTF_PROXY = pysam.ctabixproxies.GTFProxy


class TabixIterator(object):
    """Iterator that iterates over records in a tabix file using pysam."""

    def __init__(self, file_path: Path) -> None:
        self._file_path = file_path

    @property
    def _parser(self):
        """Returns parser to use for parsing tabix records."""
        raise NotImplementedError()

    @contextlib.contextmanager
    def _open_file(self) -> pysam.TabixFile:
        # Open gtf file.
        tabix_file = pysam.TabixFile(self._file_path, parser=self._parser)

        # Yield file object and ensure it is closed.
        try:
            yield tabix_file
        finally:
            tabix_file.close()

    @property
    def contigs(self) -> List[str]:
        """Contigs present in the tabix file."""

        with self._open_file() as tabix_file:
            return tabix_file.contigs

    def fetch(self,
              reference: str=None,
              start: int=None,
              end: int=None,
              filters: Iterable[Callable[[Any], bool]]=None) -> Iterable[Any]:
        """Fetches tabix records from the tabix file."""

        with self._open_file() as tabix_file:
            # For some reason pysam does not fetch all records if reference
            # is None under Python 2.7. To fix this, here we simply chain all
            # the contig records into one iterable.
            if reference is None:
                contigs = tabix_file.contigs
                records = itertools.chain.from_iterable((tabix_file.fetch(
                    reference=ref, start=start, end=end) for ref in contigs))
            else:
                records = tabix_file.fetch(
                    reference=reference, start=start, end=end)

            # Filter records on additional filters.
            if filters is not None:
                for filter_func in filters:
                    records = (rec for rec in records if filter_func(rec))

            yield from records


class GtfIterator(TabixIterator):
    """Iterator that iterates over records in a GTF file using pysam."""

    @property
    def _parser(self):
        """Returns parser to use for parsing tabix records."""
        return pysam.asGTF()

    def fetch_genes(self, gene_ids: List[str]=None, **kwargs) -> Iterable[Any]:
        """Fetches gene records from the GTF file."""

        # Fetch gene records.
        gene_filter = lambda rec: rec.feature == 'gene'
        records = self.fetch(filters=[gene_filter], **kwargs)

        if gene_ids is not None:
            # Run through records.
            sought_ids = set(gene_ids)

            hits = {}
            for rec in records:
                if rec['gene_id'] in sought_ids:
                    # Record hits.
                    hits[rec['gene_id']] = rec

                    # Stop if we found everything.
                    if len(hits) == len(gene_ids):
                        break

            # Yield results, returning None for genes that weren't found.
            for gene_id in gene_ids:
                yield hits.get(gene_id, None)

        else:
            yield from records


class BedIterator(TabixIterator):
    """Iterator that iterates over records in a BED file using pysam."""

    @property
    def _parser(self):
        return pysam.asBed()
