"""This module provides basic tracks and supporting actor classes
for plotting gene/transcript features. The module provides the
following tracks:

    - GeneTrack: Basic track for plotting genes/transcripts from
        a GtfFile or GtfFrame instance on a genomic axis.
    - BiomartGeneTrack: Gene track that downloads its annotation
        directly from Biomart, thus avoiding having to provide
        an annotation file beforehand.

"""

from __future__ import absolute_import, division, print_function

#pylint: disable=W0622,W0611
from builtins import (ascii, bytes, chr, dict, filter, hex, input, int, map,
                      next, oct, open, pow, range, round, str, super, zip)

import itertools

import pandas as pd
import numpy as np
import toolz

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import matplotlib.collections as mpl_coll

try:
    import pybiomart
except ImportError:
    pybiomart = None

from geneviz.tracks import Track
from geneviz.util.tabix import GtfFile, GtfFrame
from geneviz.util.misc import merge_intervals, get_region

from .base import Actor, StackedTrack, BoundingBox

DEFAULT_COLOR = 'grey'


class Gene(Actor):
    """
    Actor that represents a gene.

    Arguments
    ---------
        seqname : str
            Chromosome of the gene.
        start : int
            Start position of the gene.
        end : int
            End position of the gene.
        strand : int
            Strand of the gene, should either be 1 if feature is on
            forward strand, or -1 for the reverse strand.
        name : str
            Name of the gene.
        height : int
            Height of the drawn gene.

    """

    _default_plot_kws = {'edgecolor': 'black', 'lw': 1}

    #pylint: disable=too-many-arguments
    def __init__(self,
                 seqname,
                 start,
                 end,
                 strand,
                 name=None,
                 height=1,
                 **kwargs):
        super().__init__(seqname=seqname, start=start, end=end)

        self.strand = strand
        self.name = name

        self._height = height
        self._plot_kws = toolz.merge(self._default_plot_kws, kwargs)

    def get_height(self):
        return self._height

    def get_bounding_box(self, ax):
        if self.name is not None:
            txt = self._draw_label(ax)
            plt.draw()

            bbox = txt.get_window_extent()
            start = ax.transAxes.inverted().transform(bbox)[0, 0]

            txt.remove()
        else:
            start = self.start

        return BoundingBox(start, self.end, self._height, self)

    def draw(self, ax, y=0):
        patch = Rectangle(
            xy=(self.start, y),
            width=self.end - self.start,
            height=self._height,
            **self._plot_kws)
        ax.add_patch(patch)

        self._draw_label(ax, y=y)

    def _draw_label(self, ax, y=1):
        # Choose anchor depending on orientation of axes.
        reverse = ax.get_xlim()[0] > ax.get_xlim()[1]

        if reverse:
            anchor = self.end
        else:
            anchor = self.start

        # Draw and return label.
        return ax.annotate(
            xy=(anchor, y + (0.5 * self._height)),
            xycoords='data',
            xytext=(-5, 0),
            textcoords='offset points',
            s=self.name,
            fontsize=12,
            horizontalalignment='right',
            verticalalignment='center',
            clip_on=True)


class Transcript(Actor):
    """Actor that represents a gene transcript.

    Args:
        features : (pandas.DataFrame)
            Transcript features
        name : str
        height : int
        arrow_spacing : int
        arrow_size : int
        **kwargs
            Extra kwargs are used to define visual properties and are passed
            to PatchCollection, LineCollection and RegularPolyCollection to
            style the exons, junctions and junction arrows respectively.
            Examples of possible arguments are 'facecolor', 'edgecolor' and
            'linewidth'. For more arguments see the matplotlib documentation.

    """

    _default_plot_kws = {'edgecolor': 'black', 'lw': 1}

    def __init__(self,
                 features,
                 name=None,
                 height=1,
                 arrow_spacing=30,
                 arrow_size=150,
                 **kwargs):

        seqname = features['seqname'].iloc[0]
        start = features['start'].min()
        end = features['end'].max()

        super().__init__(seqname=seqname, start=start, end=end)

        self.name = name

        self._features = features
        self._height = height

        self._arrow_size = arrow_size
        self._arrow_spacing = arrow_spacing

        self._plot_kws = toolz.merge(self._default_plot_kws, kwargs)

    def get_height(self):
        return self._height

    def get_bounding_box(self, ax):
        if self.name is not None:
            txt = self._draw_label(ax)
            plt.draw()

            bbox = txt.get_window_extent()
            start = ax.transAxes.inverted().transform(bbox)[0, 0]

            txt.remove()
        else:
            start = self.start

        return BoundingBox(start, self.end, self._height, self)

    def draw(self, ax, y=1):
        # Draw gene itself.
        self._draw_exons(ax, y)
        self._draw_junctions(ax, y)

        # Draw gene label.
        if self.name is not None:
            self._draw_label(ax, y)

    def _draw_exons(self, ax, y):
        exons = self._features.query('feature == "exon"')

        patches = (Rectangle(
            xy=(f.start, y), width=f.end - f.start, height=self._height)
                   for f in exons.itertuples())
        collection = PatchCollection(patches, **self._plot_kws)

        ax.add_collection(collection)

    def _draw_junctions(self, ax, y):
        # Get junction ranges.
        exons = (self._features.query('feature == "exon"')
                 .sort_values(by='start'))

        ranges = list(zip(exons['end'].iloc[:-1], exons['start'].iloc[1:]))

        # Draw lines.
        y_line = y + (0.5 * self._height)
        segments = [((start, y_line), (end, y_line)) for start, end in ranges]
        lines = mpl_coll.LineCollection(segments, **self._plot_kws)
        #segments, colors=('grey', ), linewidths=3)
        ax.add_collection(lines)

        # Draw arrows.
        self._draw_junction_arrows(ax, segments)

    def _draw_junction_arrows(self, ax, segments):
        arrow_spacing = self._get_arrow_spacing(ax)

        offsets = []
        for ((start, height), (end, _)) in segments:
            num_arrows = (end - start) // arrow_spacing
            positions = np.linspace(start, end, num_arrows + 2)[1:-1]
            offsets += list(zip(positions, itertools.cycle([height])))

        arrows = mpl_coll.RegularPolyCollection(
            numsides=3,
            rotation=np.pi / 2,
            sizes=(self._arrow_size, ),
            offsets=offsets,
            transOffset=ax.transData,
            zorder=5,
            **self._plot_kws)
        #            facecolors=('grey',),
        #            edgecolors=('grey',), zorder=5)
        ax.add_collection(arrows)

    def _get_arrow_spacing(self, ax):
        pixel_pos = np.vstack([[0, self._arrow_spacing], [0, 0]]).T
        data_pos = ax.transData.inverted().transform(pixel_pos)
        x_pos, _ = data_pos.T
        return max(x_pos) - min(x_pos)

    def _draw_label(self, ax, y=1):
        # Choose anchor depending on orientation of axes.
        reverse = ax.get_xlim()[0] > ax.get_xlim()[1]

        if reverse:
            anchor = self.end
        else:
            anchor = self.start

        # Draw and return label.
        return ax.annotate(
            xy=(anchor, y + (0.5 * self._height)),
            xycoords='data',
            xytext=(-5, 0),
            textcoords='offset points',
            s=self.name,
            fontsize=12,
            horizontalalignment='right',
            verticalalignment='center',
            clip_on=True)


class BaseGeneTrack(Track):
    """Track for plotting gene/transcript annotations.

    Gene draws genes/transcripts on a genomic axis, in a stacked
    fashion if transcripts/genes are overlapping.

    Args:
        gtf (GtfFile or GtfFrame): Gtf file/frame containing the
            gene and/or transcript annotations.
        filters (dict): Dict of filters that are passed to the
            fetch_frame method of the GtfFile/Frame object when
            querying for genes/transripts in a given region.
        collapse (bool): Whether to collapse transcripts per gene.
        transcript_kws: Keywords that define the aesthetics
            of the transcripts, are passed to the Transcript
            actor constructor.
        transcript_id: The gtf column that contains the transcript name.
        gene_id: The gtf column that contains the gene name, is used
            when collapsing transcripts per gene.

    """

    # TODO: Fix case where plot_kws contains facecolor (currently
    # conflicts with the hue implementation).

    def __init__(self,
                 hue=None,
                 palette=None,
                 collapse=None,
                 filter=None,
                 stack_kws=None,
                 gene_id='gene_name',
                 transcript_id='transcript_name',
                 **kwargs):
        super().__init__()

        self._hue = hue
        self._palette = palette

        self._collapse = collapse
        self._filter = filter

        self._stack_kws = stack_kws or {}
        self._plot_kws = kwargs

        self._gene_id = gene_id
        self._transcript_id = transcript_id

    def _fetch_data(self, seqname, start, end):
        """
        Abstract method that fetches data for a given region.

        Should return a Dataframe containing the following columns: seqname,
        start, end, strand, feature (containing the type of feature that the
        row represents, e.g. exon, gene transcript).
        """
        raise NotImplementedError()

    def draw(self, ax, seqname, start, end):
        track = self._build_track(seqname, start, end)
        track.draw(ax, seqname, start, end)
        ax.set_xlabel('Chromosome {}'.format(seqname))

    def get_height(self, ax, seqname, start, end):
        track = self._build_track(seqname, start, end)
        return track.get_height(ax, seqname, start, end)

    def _build_track(self, seqname, start, end):
        data = self._fetch_data(seqname, start, end)
        data = self._preprocess_data(data)
        features = self._convert_to_features(data)
        return StackedTrack(features, **self._stack_kws)

    def _preprocess_data(self, data):
        # Filter data if needed.
        if self._filter is not None:
            data = data.query(self._filter)

        # Ensure we are modifying a copy.
        data = data.copy()

        # Add color.
        if self._hue is not None:
            data['color'] = data[self._hue].map(self._palette)
        else:
            data['color'] = DEFAULT_COLOR

        # Copy id fields to standard locations.
        data['transcript_id'] = data[self._transcript_id]
        data['gene_id'] = data[self._gene_id]

        return data

    def _convert_to_features(self, data):
        # Collapse transcripts if needed.
        if self._collapse == 'gene':
            groups = (data.dropna(subset=['gene_id']).groupby('gene_id'))

            features = [self._data_to_gene(grp) for _, grp in groups]
        else:
            if self._collapse == 'transcript':
                data = self._collapse_transcripts(data)

            # Create list of transcripts.
            if len(data) > 0:
                groups = (data.dropna(subset=['transcript_id'])
                          .groupby('transcript_id'))

                features = [self._data_to_transcript(grp) for _, grp in groups]
            else:
                features = []

        return features

    def _data_to_gene(self, grp):
        first = grp.iloc[0]
        return Gene(
            seqname=first.seqname,
            start=grp['start'].min(),
            end=grp['end'].max(),
            strand=first.strand,
            name=first.gene_id,
            facecolor=first.color,
            **self._plot_kws)

    def _data_to_transcript(self, grp):
        first = grp.iloc[0]
        return Transcript(
            grp,
            name=first.transcript_id,
            facecolor=first.color,
            **self._plot_kws)

    def _collapse_transcripts(self, data):
        # TODO: handle empty case.
        collapsed = (self._collapse_transcript(grp)
                     for gene_id, grp in data.groupby('gene_id'))
        return pd.concat(collapsed, ignore_index=True)

    def _collapse_transcript(self, data):
        # Subset transcript to exons.
        data = data.ix[data['feature'] == 'exon']

        if len(data) == 0:
            # Return empty dataframe if no exons in transcript.
            collapsed = pd.DataFrame([], columns=data.columns)
        else:
            # Otherwise, merge overlapping transcript regions.
            regions = zip(data['start'], data['end'])
            start, end = zip(*merge_intervals(regions))

            # And build a new frame, using the first row as reference.
            first = data.iloc[0]

            collapsed = (pd.DataFrame([first] * len(start)).assign(
                start=start, end=end, transcript_id=first.gene_id))

        return collapsed


class GeneTrack(BaseGeneTrack):
    def __init__(self, data, **kwargs):
        super().__init__(**kwargs)
        self._data = data

    def _fetch_data(self, seqname, start, end):
        # Perform initial query.
        data = get_region(self._data, seqname, start, end)

        # Check if we have any features extending beyond range.
        min_ = data['start'].min()
        max_ = data['end'].max()

        if min_ < start or max_ > end:
            # We have features extending outside queried range,
            # so we need to peform a second query to get all features.

            # Perform second query in extended range, but limit
            # new features to the original transcripts.
            transcripts = set(data[self._transcript_id])

            data = get_region(data, seqname, min(min_, start), max(max_, end))
            data = data.ix[data[self._transcript_id].isin(transcripts)]

        return data


class GtfTrack(BaseGeneTrack):
    def __init__(self, gtf, **kwargs):
        super().__init__(**kwargs)
        self._gtf = gtf

    @classmethod
    def from_path(cls, gtf_path, **kwargs):
        return cls(GtfFile(gtf_path), **kwargs)

    def _fetch_data(self, seqname, start, end):
        # Perform initial query.
        data = self._gtf.fetch_frame(seqname, start, end)

        # Check if we have any features extending beyond range.
        min_ = data['start'].min()
        max_ = data['end'].max()

        if min_ < start or max_ > end:
            # We have features extending outside queried range,
            # so we need to peform a second query to get all features.

            # Perform second query in extended range, but limit
            # new features to the original transcripts.
            transcripts = set(data[self._transcript_id])

            data = self._gtf.fetch_frame(data, seqname, min(min_, start),
                                         max(max_, end))

            data = data.ix[data[self._transcript_id].isin(transcripts)]

        # Transform to correct format.
        data = (data
                .rename(columns={'contig': 'seqname'})
                .assign(strand=lambda df: df['strand'].map({'+': 1, '-': -1})))

        return data


class BiomartTrack(BaseGeneTrack):
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
                 **kwargs):
        super().__init__(**kwargs)

        if pybiomart is None:
            raise ValueError('Pybiomart must be installed to use '
                             'the BiomartGeneTrack ')

        # Fetch dataset from server.
        server = pybiomart.Server(host=host)
        self._dataset = (server.marts[mart].datasets[dataset])

    def _fetch_data(self, seqname, start, end):
        # TODO: Fetch transcript name instead of id (needs extra query).

        # Determine which transcripts are within region.
        transcripts = self._dataset.query(
            attributes=['ensembl_transcript_id'],
            filters={'chromosome_name': [seqname],
                     'start': [start],
                     'end': [end]})

        transcript_ids = list(transcripts[transcripts.columns[0]])

        # Retrieve exons for these transcripts.
        attrs = ['external_gene_name', 'ensembl_gene_id',
                 'ensembl_transcript_id', 'exon_chrom_start', 'exon_chrom_end',
                 'strand']

        filters = {'link_ensembl_transcript_stable_id': transcript_ids}

        data = self._dataset.query(attributes=attrs, filters=filters)
        data.columns = attrs

        # Reshape data to 'conform' to gtf format.
        data = pd.DataFrame({
            'seqname': seqname,
            'start': data['exon_chrom_start'],
            'end': data['exon_chrom_end'],
            'strand': data['strand'],
            'feature': 'exon',
            'gene_name': data['external_gene_name'],
            'gene_id': data['ensembl_gene_id'],
            'transcript_name': data['ensembl_transcript_id']
        })

        return data
