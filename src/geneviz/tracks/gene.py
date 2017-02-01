"""This module provides basic tracks and supporting actor classes
for plotting gene/transcript features. The module provides the
following tracks:

    - GeneTrack: Basic track for plotting genes/transcripts from
        a GtfFile or GtfFrame instance on a genomic axis.
    - BiomartGeneTrack: Gene track that downloads its annotation
        directly from Biomart, thus avoiding having to provide
        an annotation file beforehand.

"""

import pandas as pd
import pybiomart
import toolz

from geneviz.tracks.base import Track
from geneviz.util.genomic import merge_intervals
from geneviz.util.tabix import GtfIterator

from .feature import FeatureTrack


class _BaseGeneTrack(Track):
    def __init__(self,
                 gene_id='gene_id',
                 transcript_id='transcript_id',
                 filter=None,
                 collapse=None,
                 hue=None,
                 hue_order=None,
                 palette=None,
                 height=0.9,
                 spacing=0.05,
                 label_kws=None):
        super().__init__()

        self._gene_id = gene_id
        self._transcript_id = transcript_id
        self._collapse = collapse
        self._filter = filter

        # Determine which grouping feature to use.
        collapse_group_map = {
            'transcript': gene_id,
            'gene': gene_id,
            None: transcript_id
        }
        group = collapse_group_map[collapse]

        # Setup kws for feature track.
        self._track_kws = {
            'hue': hue,
            'hue_order': hue_order,
            'palette': palette,
            'height': height,
            'spacing': spacing,
            'label': group,
            'group': group,
            'strand_junctions': True,
            'label_kws': label_kws
        }

    def get_height(self, region, ax):
        track = self._build_track(region)
        return track.get_height(region, ax)

    def _build_track(self, region):
        data = self._fetch_exons(region)
        return FeatureTrack(data=data, **self._track_kws)

    def _fetch_data(self, region):
        raise NotImplementedError()

    def _fetch_exons(self, region):
        exons = self._fetch_data(region)

        if self._filter is not None:
            exons = exons.query(self._filter)

        if self._collapse == 'gene':
            exons = self._collapse_gene(exons)
        elif self._collapse == 'transcript':
            exons = self._collapse_transcripts(exons)
        elif self._collapse is not None:
            raise ValueError('Unexpected value for collapse')

        return exons

    def _collapse_gene(self, exons):
        """Collapses exons into a single gene body."""

        return exons.groupby(self._gene_id).agg({
            'chromosome': 'first',
            'start': 'min',
            'end': 'max',
            'strand': 'first'
        }).reset_index()

    def _collapse_transcripts(self, exons):
        """Collapses multiple transcripts into a single transcript."""

        def _collapse_grp(grp):
            first = grp.iloc[0]
            starts, ends = zip(*merge_intervals(zip(grp.start, grp.end)))
            return pd.DataFrame({
                'chromosome': first.chromosome,
                'start': starts,
                'end': ends,
                'strand': first.strand,
                self._gene_id: first[self._gene_id]
            })

        return pd.concat(
            (_collapse_grp(grp) for _, grp in exons.groupby(self._gene_id)),
            axis=0,
            ignore_index=True)

    def draw(self, region, ax):
        track = self._build_track(region)
        return track.draw(region, ax)


class GeneTrack(_BaseGeneTrack):
    def __init__(self,
                 data,
                 gene_id='gene_id',
                 transcript_id='transcript_id',
                 collapse=None,
                 filter=None,
                 hue=None,
                 hue_order=None,
                 palette=None,
                 height=0.9,
                 spacing=0.05,
                 label_kws=None):

        super().__init__(
            gene_id=gene_id,
            transcript_id=transcript_id,
            collapse=collapse,
            filter=filter,
            hue=hue,
            hue_order=hue_order,
            palette=palette,
            height=height,
            spacing=spacing,
            label_kws=label_kws)
        self._data = data

    def _fetch_data(self, region):
        raise NotImplementedError()


class GtfTrack(_BaseGeneTrack):
    def __init__(self,
                 gtf_path,
                 gene_id='gene_id',
                 transcript_id='transcript_id',
                 collapse=None,
                 filter=None,
                 hue=None,
                 hue_order=None,
                 palette=None,
                 height=0.9,
                 spacing=0.05,
                 label_kws=None):
        super().__init__(
            gene_id=gene_id,
            transcript_id=transcript_id,
            collapse=collapse,
            filter=filter,
            hue=hue,
            hue_order=hue_order,
            palette=palette,
            height=height,
            spacing=spacing,
            label_kws=label_kws)
        self._gtf_path = gtf_path

    def _fetch_data(self, region):
        # TODO: Ensure we fetch the full transcript.

        # Fetch exon records from gtf.
        gtf_iter = GtfIterator(self._gtf_path)
        records = gtf_iter.fetch(*region)

        # Convert to DataFrame.
        exons = pd.DataFrame.from_records(self._records_to_exons(records))

        return exons

    @classmethod
    def _records_to_exons(cls, records):
        for record in records:
            if record.feature == 'exon':
                yield cls._parse_exon_record(record)

    @staticmethod
    def _parse_exon_record(record):
        basic_info = {
            'chromosome': record.contig,
            'start': record.start,
            'end': record.end,
            'strand': numeric_strand(record.strand)
        }
        return toolz.merge(basic_info, dict(record))


def numeric_strand(strand):
    """Converts the +/- strand notation to a numeric 1/-1 format."""
    assert strand in {'+', '-'}
    return 1 if strand == '+' else -1


class BiomartTrack(_BaseGeneTrack):
    """Track for plotting gene/transcript annotations from Biomart.

    Gene draws genes/transcripts on a genomic axis, in a stacked
    fashion if transcripts/genes are overlapping. Gene annotations
    are queried directly from Biomart.

    Args:
        host (str): Biomart host address.
        mart (str): Biomart mart name.
        dataset (str): Biomart dataset name.
        **kwargs: Keywords are passed to the GeneTrack constructor,
            which is used internally to draw transcripts/genes.

    """

    def __init__(self,
                 host='http://www.ensembl.org',
                 mart='ENSEMBL_MART_ENSEMBL',
                 dataset='hsapiens_gene_ensembl',
                 gene_id='gene_id',
                 transcript_id='transcript_id',
                 collapse=None,
                 filter=None,
                 hue=None,
                 hue_order=None,
                 palette=None,
                 height=0.9,
                 spacing=0.05,
                 label_kws=None):
        super().__init__(
            gene_id=gene_id,
            transcript_id=transcript_id,
            collapse=collapse,
            filter=filter,
            hue=hue,
            hue_order=hue_order,
            palette=palette,
            height=height,
            spacing=spacing,
            label_kws=label_kws)

        if pybiomart is None:
            raise ValueError('Pybiomart must be installed to use '
                             'the BiomartGeneTrack ')

        # Fetch dataset from server.
        server = pybiomart.Server(host=host)
        self._dataset = (server.marts[mart].datasets[dataset])

    def _fetch_data(self, region):
        # TODO: Fetch transcript name instead of id (needs extra query).

        chromosome, start, end = region

        # Determine which transcripts are within region.
        transcripts = self._dataset.query(
            attributes=['ensembl_transcript_id'],
            filters={
                'chromosome_name': [chromosome],
                'start': [start],
                'end': [end]
            })
        transcript_ids = list(transcripts[transcripts.columns[0]])

        # Retrieve exons for these transcripts.
        data = self._dataset.query(
            attributes=[
                'external_gene_name', 'ensembl_gene_id',
                'ensembl_transcript_id', 'exon_chrom_start', 'exon_chrom_end',
                'strand'
            ],
            filters={'link_ensembl_transcript_stable_id': transcript_ids},
            use_attr_names=True)

        # Reshape data to conform to expected format.
        data = pd.DataFrame({
            'chromosome': chromosome,
            'start': data['exon_chrom_start'],
            'end': data['exon_chrom_end'],
            'strand': data['strand'],
            'gene_name': data['external_gene_name'],
            'gene_id': data['ensembl_gene_id'],
            'transcript_id': data['ensembl_transcript_id']
        })

        return data
