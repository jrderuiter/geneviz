# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import numpy as np
import pytest

from geneviz.util import tabix

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


@pytest.fixture
def gtf_path():
    return pytest.helpers.data_path('mm10.test.gtf.gz')


@pytest.fixture
def tabix_file(gtf_path):
    return tabix.GtfFile(gtf_path)


class TestTabixFile(object):
    def test_fetch(self, tabix_file):
        records = list(tabix_file.fetch())

        # Assert we fetch all records and check if
        # we have records from both contigs.
        assert len(records) > 0
        assert set(r.contig for r in records) == {'1', '11'}
        assert set(r.strand for r in records) == {'+', '-'}

    def test_fetch_reference(self, tabix_file):
        records = list(tabix_file.fetch('1'))

        # Test if we only fetch records from chrom 1.
        assert len(records) > 0
        assert set(r.contig for r in records) == {'1'}

    def test_fetch_range(self, tabix_file):
        records = list(tabix_file.fetch('1', 182407167, 182464436))

        # Test if all records are within selected range.
        assert len(records) > 0
        assert set(r.contig for r in records) == {'1'}
        assert min(r.start for r in records) >= 182407167
        assert max(r.end for r in records) <= 182464436

    def test_fetch_filter(self, tabix_file):
        records = list(tabix_file.fetch(filters={'feature': 'transcript'}))

        # Test if all records are transcripts.
        assert len(records) > 0
        assert all(r.feature == 'transcript' for r in records)

    def test_fetch_filter2(self, tabix_file):
        records = list(tabix_file.fetch(filters={'strand': '+'}))

        # Test if all records are from the forward strand.
        assert len(records) > 0
        assert all(r.strand == '+' for r in records)

    def test_fetch_incl_left_true(self, tabix_file):
        records = list(
            tabix_file.fetch(
                '1', 182409431, 182464436, incl_left=True))

        # Check if record overlapping left is still included.
        assert len(records) > 0
        assert set(r.contig for r in records) == {'1'}
        assert 182409171 <= min(r.start for r in records) <= 182409431
        assert max(r.end for r in records) <= 182464436

    def test_fetch_incl_left_false(self, tabix_file):
        records = list(
            tabix_file.fetch(
                '1', 182409431, 182464436, incl_left=False))

        # Check if no record falls outside left bound.
        assert len(records) > 0
        assert set(r.contig for r in records) == {'1'}
        assert min(r.start for r in records) >= 182409431
        assert max(r.end for r in records) <= 182464436

    def test_fetch_incl_right_true(self, tabix_file):
        records = list(
            tabix_file.fetch(
                '1', 182409431, 182461830, incl_right=True))

        # Check if record overlapping right is still included.
        assert len(records) > 0
        assert set(r.contig for r in records) == {'1'}
        assert min(r.start for r in records) >= 182409171
        assert 182461830 <= max(r.end for r in records) <= 182464436

    def test_fetch_incl_right_false(self, tabix_file):
        records = list(
            tabix_file.fetch(
                '1', 182409431, 182461830, incl_right=False))

        # Check if no record falls outside right bound.
        assert len(records) > 0
        assert set(r.contig for r in records) == {'1'}
        assert min(r.start for r in records) >= 182409171
        assert max(r.end for r in records) <= 182461830


class TestGtfFile(object):
    def test_get_gene(self, gtf_path):
        """Tests retrieval of a gene."""

        gtf = tabix.GtfFile(gtf_path)
        gene = gtf.get_gene('ENSMUSG00000026510')

        # Check gtf fields.
        assert gene['contig'] == '1'
        assert gene['source'] == 'protein_coding'
        assert gene['feature'] == 'gene'
        assert gene['start'] == 182409172
        assert gene['end'] == 182462432
        assert gene['strand'] == '+'
        assert np.isnan(gene['score'])
        assert gene['frame'] == '.'

        # Check extra fields.
        assert gene['gene_id'] == 'ENSMUSG00000026510'
        assert gene['gene_name'] == 'Trp53bp2'
        assert gene['gene_biotype'] == 'protein_coding'
        assert gene['gene_source'] == 'ensembl_havana'

    def test_get_gene_by_name(self, gtf_path):
        """Tests retrieval of a gene by its name."""

        gtf = tabix.GtfFile(gtf_path)
        gene = gtf.get_gene('Trp53bp2', field_name='gene_name')

        assert gene['gene_name'] == 'Trp53bp2'

    def test_get_gene_missing(self, gtf_path):
        """Tests if retrieval of a missing gene raises a ValueError."""

        gtf = tabix.GtfFile(gtf_path)

        with pytest.raises(ValueError):
            gtf.get_gene('ENSMUSG00000000000')
