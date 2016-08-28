import itertools

import toolz
import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches

from .base import Actor, StackedTrack, BoundingBox
from .. import Track


class Feature(Actor):
    """Actor that represents a generic (single) genomic feature.

    Feature objects can be used to represent generic features. A
    feature object may be undirected (strand = None), in which case
    it is drawn as a Rectangle, or directed (strand = 1 or strand = -1),
    in which case the feature is drawn as a directed arrow.

    Args:
        seqname (str): Chromosome of feature.
        start (int): Start position of feature.
        end (int): End position of feature.
        strand (int, optional): Strand of feature, should either be
            1 if feature is on forward strand, or -1 for the reverse strand.
        label (str, optional): Label of the feature.
        arrow_kws (dict): Keywords specifying how the arrow should be drawn
            for a directed feature. Keywords are passed directly to
            matplotlibs FancyArrow constructor.
        rect_kws (dict): Keywords specifying how the rectangle should be
            drawn of an undirected feature. Keywords are passed directly
            to matplotlibs Rectangle constructor.
        shared_kws (dict): Keywords specifying how features should be
            drawn, but that are not specific to directed/undirected features.
            Examples are feature colors, etc.

    """

    _default_plot_kws = {'edgecolor': 'black', 'lw': 1}

    def __init__(self,
                 seqname,
                 start,
                 end,
                 strand=None,
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

    def draw(self, ax, y=1):
        if self.strand is None or np.isnan(self.strand):
            # Draw rectangle.
            # height = self._rect_kws['height']
            patch = mpl_patches.Rectangle(
                xy=(self.start, y),
                width=self.end - self.start,
                height=self._height,
                **self._plot_kws)
        else:
            # Draw directed arrow.
            if self.strand == 1:
                x, dx = self.start, self.end - self.start
            else:
                x, dx = self.end, self.start - self.end

            arrow_kws = dict(
                head_width=0.9 * self._height,
                head_length=abs(self.start - self.end) * 0.5,
                width=0.5 * self._height,
                length_includes_head=True)

            patch = mpl_patches.FancyArrow(
                x=x,
                dx=dx,
                y=y + (0.5 * self._height),
                dy=0,
                **toolz.merge(arrow_kws, self._plot_kws))

        ax.add_patch(patch)

        if self.name is not None:
            self._draw_label(ax, y=y)

    def _draw_label(self, ax, y=1):
        return ax.annotate(
            xy=(self.start, y + (0.5 * self._height)),
            xycoords='data',
            xytext=(-5, 0),
            textcoords='offset points',
            s=self.name,
            fontsize=16,
            horizontalalignment='right',
            verticalalignment='center',
            clip_on=True)


class FeatureTrack(Track):
    """Track for plotting generic genomic features.

    FeatureTrack draws a list of features on a given track axis,
    using the StackedTrack to stack overlapping features. The main
    utility of the class lises in the `from_frame` class method,
    which constructs a feature track directly from a dataframe.

    Args:
        data (pandas.Dataframe): DataFrame containing the features to
            draw. Should have the following columns: seqname, start,
            end, strand. Each row is taken as a separate feature.

    """

    # seqname, start, end, strand=None, label=None,
    # arrow_kws=None, rect_kws=None, **kwargs)

    def __init__(self,
                 data,
                 name=None,
                 hue=None,
                 palette=None,
                 stack_kws=None,
                 height=1):
        super().__init__()

        palette = palette or sns.color_palette()

        self._data = self._preprocess_data(
            data, name=name, hue=hue, palette=palette)
        self._stack_kws = stack_kws or {}

        self._height = height

    def _preprocess_data(self, data, name, hue, palette):
        plot_data = data[['seqname', 'start', 'end', 'strand']].copy()

        if hue is not None:
            if not isinstance(palette, dict):
                palette = dict(
                    zip(data[hue].unique(), itertools.cycle(palette)))
            plot_data['color'] = data[hue].map(palette)
        else:
            plot_data['color'] = None

        if name is not None:
            plot_data['name'] = data[name]
        else:
            plot_data['name'] = None

        return plot_data

    def draw(self, ax, seqname, start, end):
        features = self._get_features(seqname, start, end)
        track = StackedTrack(features, **self._stack_kws)
        return track.draw(ax, seqname, start, end)

    def get_height(self, ax, seqname, start, end):
        features = self._get_features(seqname, start, end)
        track = StackedTrack(features, **self._stack_kws)
        return track.get_height(ax, seqname, start, end)

    def _get_features(self, seqname, start, end):
        sel_data = self._data.query(('seqname == {!r} and '
                                     'end >= {} and start <= {}')
                                    .format(seqname, start, end))

        for tup in sel_data.itertuples():
            yield self._feature_from_row(tup)

    def _feature_from_row(self, row):
        return Feature(
            seqname=row.seqname,
            start=row.start,
            end=row.end,
            strand=row.strand,
            facecolor=row.color,
            name=row.name,
            height=self._height)
