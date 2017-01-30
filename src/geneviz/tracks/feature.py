from collections import namedtuple
import itertools
import operator

import pandas as pd
import numpy as np
import seaborn as sns
import toolz

from intervaltree import IntervalTree
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.collections import PatchCollection, LineCollection

from geneviz.util.colormap import build_colormap

from .base import Track


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

    def __init__(self,
                 data,
                 group=None,
                 label=None,
                 hue=None,
                 hue_order=None,
                 palette=None,
                 strand_junctions=False,
                 height=1,
                 spacing=0.05,
                 color='dimgrey',
                 patch_kws=None,
                 line_kws=None,
                 label_kws=None):
        super().__init__()

        # TODO: Add legend for hue.

        # Setup internal dataframe.
        color_map = build_colormap(
            data, hue=hue, palette=palette, order=hue_order)

        self._data = data
        self._color_map = color_map

        # Various visual parameters.
        self._hue = hue
        self._group = group
        self._label = label

        self._height = height
        self._spacing = spacing

        self._strand_junctions = strand_junctions

        # Detailed style kws for different plot aspects.
        default_patch_kws = {'facecolor': color, 'edgecolor': color, 'lw': 1}
        self._patch_kws = toolz.merge(default_patch_kws, patch_kws or {})

        default_line_kws = {'color': color, 'lw': 1.5}
        self._line_kws = toolz.merge(default_line_kws, line_kws or {})

        self._label_kws = label_kws or {}

    @classmethod
    def from_position(cls, data, width, **kwargs):
        """Constructs track from frame with positions instead of start/ends."""

        data = data.assign(
            start=data['position'] - (width // 2),
            end=data['position'] + (width // 2))

        return cls(data=data, **kwargs)

    def get_height(self, region, ax):
        """Returns the height of the track."""

        data = self._fetch_data(region).assign(height=self._height)
        stacked = stack(
            data,
            group=self._group,
            label=self._label,
            label_func=self._get_label_extent,
            ax=ax,
            spacing=self._spacing)

        return stacked['y'].max() + self._height + self._spacing

    def _fetch_data(self, region):
        """Fetches features within a given region."""

        return self._data.query(
            ('chromosome == {!r} and end >= {} and start <= {}')
            .format(*region))  # yapf: disable

    def draw(self, region, ax):
        """Draws the track."""

        # Fetch data within region.
        data = self._fetch_data(region).assign(height=self._height)
        stacked = stack(
            data,
            group=self._group,
            label=self._label,
            label_func=self._get_label_extent,
            ax=ax,
            spacing=self._spacing)

        # Draw features.
        if self._hue is None:
            patches = (self._feature_patch(tup)
                       for tup in stacked.itertuples())
            ax.add_collection(PatchCollection(patches, **self._patch_kws))
        else:
            for hue, grp in stacked.groupby(self._hue):
                patches = (self._feature_patch(tup)
                           for tup in grp.itertuples())
                patch_kws = toolz.merge(self._patch_kws,
                                        {'facecolor': self._color_map[hue]})
                ax.add_collection(PatchCollection(patches, **patch_kws))

        # Draw junctions/labels.
        if self._group is not None:
            segments = []
            for _, grp in stacked.groupby(self._group):
                # Collect junction segments.
                segments += list(self._junction_segments(grp))

                # Draw labels for groups.
                if self._label is not None:
                    self._draw_label_group(grp, ax=ax)

            ax.add_collection(LineCollection(segments, **self._line_kws))
        else:
            # Draw labels for single features.
            if self._label is not None:
                for tup in stacked.itertuples():
                    self._draw_label_single(tup, ax)

        # Set ylim and style axes.
        ax.set_ylim(0, stacked['y'].max() + self._height)
        ax.set_yticks([])

    def _feature_patch(self, tup):
        if self._strand_junctions or pd.isnull(tup.strand):
            # Return rectangle.
            return Rectangle(
                xy=(tup.start, tup.y),
                width=tup.end - tup.start,
                height=tup.height)
        else:
            # Return directed arrow.
            if tup.strand == 1:
                x, dx = tup.start, tup.end - tup.start
            else:
                x, dx = tup.end, tup.start - tup.end

            arrow_kws = dict(
                head_width=0.9 * tup.height,
                head_length=abs(tup.start - tup.end) * 0.5,
                width=0.5 * tup.height,
                length_includes_head=True)

            return FancyArrow(
                x=x, dx=dx, y=tup.y + (0.5 * tup.height), dy=0, **arrow_kws)

    def _junction_segments(self, grp):
        first = grp.iloc[0]
        y = first.y + (0.5 * first.height)

        if self._strand_junctions and not pd.isnull(first.strand):
            # Plot stranded junction with offset middle.
            mid_offset = 0.25 * first.height
            for start, end in zip(grp.end[:-1], grp.start[1:]):
                mid = (start + end) / 2
                y_mid = y + (mid_offset * first.strand)
                yield ((start, y), (mid, y_mid), (end, y))
        else:
            # Plot unstranded junction as straight line.
            for start, end in zip(grp.end[:-1], grp.start[1:]):
                yield ((start, y), (end, y))

    def _draw_label_group(self, grp, ax):
        # Check if we are drawing in reverse.
        xlim = ax.get_xlim()
        reverse = xlim[1] < xlim[0]

        # Use end as anchor for reverse situation.
        first = grp.iloc[0]
        x = grp.end.max() if reverse else grp.start.min()
        y = first.y + (0.5 * first.height)

        # Draw label.
        self._draw_label(x=x, y=y, label=first[self._label], ax=ax)

    def _draw_label_single(self, tup, ax):
        # Check if we are drawing in reverse.
        xlim = ax.get_xlim()
        reverse = xlim[1] < xlim[0]

        # Use end as anchor for reverse situation.
        x = tup.end if reverse else tup.start
        y = tup.y + (0.5 * tup.height)

        # Draw label.
        self._draw_label(x=x, y=y, label=getattr(tup, self._label), ax=ax)

    def _draw_label(self, x, y, label, ax):
        return ax.annotate(
            xy=(x, y),
            xycoords='data',
            xytext=(-5, 0),
            textcoords='offset points',
            s=label,
            horizontalalignment='right',
            verticalalignment='center',
            clip_on=True,
            **self._label_kws)

    def _get_label_extent(self, x, y, label, ax):
        txt = self._draw_label(x=x, y=y, label=label, ax=ax)
        plt.draw()

        bbox = txt.get_window_extent()
        dbox = bbox.transformed(ax.transData.inverted())

        txt.remove()

        return dbox


class RugTrack(Track):
    """Track that plots density ticks for features."""

    def __init__(self, data, hue=None, palette=None, height=1.0):
        super().__init__()
        self._data = self._preprocess_data(data, hue, palette)
        self._height = height

    def _preprocess_data(self, data, hue, palette):
        plot_data = data[['chromosome', 'position', 'strand']].copy()

        if hue is not None:
            if not isinstance(palette, dict):
                palette = dict(
                    zip(data[hue].unique(), itertools.cycle(palette)))
            plot_data['color'] = data[hue].map(palette)
        else:
            plot_data['color'] = None

        return plot_data

    def get_height(self, region, ax):
        return self._height

    def draw(self, region, ax):
        data = self._data.query(
            'chromosome == {!r} and position > {} and position < {}'
            .format(*region))  # yapf: disable

        for row in data.itertuples():
            ax.axvline(row.position, color=row.color)

        ax.yaxis.set_visible(False)


BoundingBox = namedtuple('BoundingBox', ['start', 'end', 'height'])


def pack_ffdh(bounding_boxes):
    """Implementation of the First-Fit Decreasing Height algorithm.

    See http://cgi.csc.liv.ac.uk/~epa/surveyhtml.html for a description
    of the algorithm and other options for packing algorithms.

    """

    # Sort boxes by decreasing height and size.
    bounding_boxes = sorted(
        bounding_boxes, key=lambda box: (box.height, box.end - box.start))
    bounding_boxes = bounding_boxes[::-1]

    levels = []
    level_heights = []

    for box in bounding_boxes:
        placed = False

        for level in levels:
            if not level.overlaps(box.start, box.end):
                level.addi(box.start, box.end, box)
                placed = True
                break

        if not placed:
            new_level = IntervalTree.from_tuples([(box.start, box.end, box)])
            levels.append(new_level)
            level_heights.append(box.height)

    levels_list = [[interval[2] for interval in level] for level in levels]

    return levels_list, level_heights


def stack(data, group=None, label=None, label_func=None, ax=None,
          spacing=0.05):
    """Stacks features in data frame."""

    # TODO: Refactor out label function? (Not a core feature of stack).

    if group is not None:
        agg_funcs = {'start': min, 'end': max, 'height': max}

        if label is not None:
            agg_funcs[label] = 'first'

        grouped = (data.groupby(group).agg(agg_funcs))

        heights = _stack(
            grouped,
            label=label,
            label_func=label_func,
            ax=ax,
            spacing=spacing)
        heights.index.name = group
        heights = heights.reset_index()

        return pd.merge(data, heights[[group, 'y']], on=group, how='outer')
    else:
        heights = _stack(data, label=label, label_func=label_func, ax=ax)
        return pd.concat([data, heights], axis=1)


def _stack(data, label=None, label_func=None, ax=None, spacing=0.05):
    if label is not None:
        data = _augment_with_labels(data, label, label_func, ax)

    levels, level_heights = pack_ffdh(data.itertuples())

    level_offsets = np.cumsum(level_heights, dtype=np.float) - level_heights[0]
    level_offsets += np.arange(1, len(level_offsets) + 1) * spacing

    rows = []
    for objs, offset in zip(levels, level_offsets):
        for obj in objs:
            rows.append((obj.Index, offset))

    index, y_values = zip(*rows)
    heights = pd.DataFrame({'y': y_values}, index=index)

    return heights


def _augment_with_labels(data, label, label_func, ax):
    # Augment positions.
    anchor_col = 'end' if reversed_axis(ax) else 'start'

    positions = [
        label_func(
            x=x, y=0, label=lab, ax=ax).x0
        for lab, x in zip(data[label], data[anchor_col])
    ]

    return data.assign(**{anchor_col: positions})


def reversed_axis(ax):
    """Checks if x-axis is reversed."""
    xlim = ax.get_xlim()
    return xlim[1] < xlim[0]
