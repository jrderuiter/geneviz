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
    """Track for plotting generic genomic features in a stacked fashion.

    The FeatureTrack is the main track that is used to plot genomic features.
    Overlapping features are plotted in a stacked fashion to avoid
    overplotting. Optionally, features can be grouped together if they belong
    to the same group, such as exons for a given gene/transcript. Grouped
    features are connected by junctions to indicate that they belong to a
    single group. Features can also be colored by their properties, using
    the hue and palette arugments.

    Parameters
    ----------
        data : pandas.Dataframe
            Dataset for plotting. Each row in the DataFrame is expected to
            correspond with a single feature. The DataFrame should have the
            following columns: seqname, start, end, strand; which together
            specify the location and orientation of the corresponding feature.
        group : str
            Column (categorical) by which features should be grouped. Grouped
            features are drawn together (at the same height) and are visually
            connected by junctions.
        label : str
            Column containing text labels that should be used to name features
            in the plot. If group is specified, the label of the first feature
            in the group is used for the group.
        hue : str
            Column (categorical) that should be used to determine the color
            of a given feature.
        hue_order : List[str]
            Order to plot the categorical hue levels in, otherwise the levels
            are inferred from the data objects.
        palette : List[Union[str, Tuple[float, float, float]]]
            Colors to use for the different levels of the hue variable.
            Should be specified as a list of colors (strs) or a list of
            tuples with RGB values (similar to Seaborn color palettes).
        strand_junctions : bool
            Boolean that indicates if group junctions should be drawn in
            a stranded fashion for groups. If False, the features are drawn
            with their own respective strands. If True, the strand of a group
            is determined from the first element of the group and the junctions
            are drawn upwards (for groups on the forward strand) or downwards
            (for groups on the reverse strand).
        height : float
            The height of an individual feature.
        spacing : float
            The (vertical) spacing to use between features.
        color : str
            Color for all of the elements. Overridded by hue.
        patch_kws : dict[str, Any]
            Dict of keyword arguments to pass to RectangleCollection or
            FancyArrow collection when drawing feature patches. Used to
            precisely specify the aesthetics of features.
        line_kws : dict[str, Any]
            Dict of keyword arguments to pass to LineCollection when drawing
            the junctions between groups. Used to specify modify the
            aesthetics of the junctions.
        label_kws : dict[str, Any]
            Dict of keyword arguments to pass to ax.annotate when drawing
            feature/group labels. Used to precisely specify the aesthetics
            of the labels.

    """

    def __init__(self,
                 data,
                 group=None,
                 label=None,
                 hue=None,
                 hue_order=None,
                 palette=None,
                 strand_junctions=False,
                 height=1.0,
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
        default_patch_kws = {
            'facecolor': color,
            'edgecolor': color,
            'linewidth': 1
        }
        self._patch_kws = toolz.merge(default_patch_kws, patch_kws or {})

        default_line_kws = {'color': color, 'linewidth': 1.5}
        self._line_kws = toolz.merge(default_line_kws, line_kws or {})

        self._label_kws = label_kws or {}

    @classmethod
    def from_position(cls, data, width, **kwargs):
        """Constructs instance from frame with positions instead of start/ends.

        Assumes that the DataFrame contains a 'position' column, which defines
        the exact position of the given feature. This function expands features
        so that the width of the plotted feature is equal to the given width
        and features are centered around their position.

        Parameters
        ----------
        data : pandas.Dataframe
            Dataset for plotting. Assumed to be the same format as for the
            main constructor, apart from containing a 'position' column,
            rather than start/end columns.
        width : int
            The width to use for the expanded features.
        **kwargs
            Any kwargs are passed to the main constructor.

        """

        data = data.assign(
            start=data['position'] - (width // 2),
            end=data['position'] + (width // 2))

        return cls(data=data, **kwargs)

    def get_height(self, region, ax):
        """Returns the height of the dummy track.

        Parameters
        ----------
        region : Tuple[str, int, int]
            The genomic region that will be drawn. Specified as a tuple of
            (chromosome, start, end).
        ax : matplotlib.Axes
            Axis that the track will be drawn on.

        Returns
        -------
        height : int
            Height of the track.

        """

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
        """Draws the track on the given axis.

        Parameters
        ----------
        region : Tuple[str, int, int]
            Genomic region to draw.
        ax : matplotlib.Axes
            Axis to draw track on.

        """

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
        ax.set_ylim(0, stacked['y'].max() + self._height + self._spacing)
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
        grp = grp.sort_values(by='start')

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


def stack(data, group=None, label=None, label_func=None, ax=None,
          spacing=0.05):
    """Stacks features in given dataframe using the FFDH algorithm.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame containing the features to be stacked.
    group : str
        Name of categorical column used to group features.
    label : str
        Name of categorical column containing feature labels.
    label_func : Function
        Function that is used to draw feature labels. This is used to assess
        the size of each label in data coordinates, which is required when
        accounting for the label sizes when stacking features.
    ax : matplotlib.Axes
        Axes on which feature labels will be drawn. Also used to assess
        feature label size.
    spacing : float
        Amount of vertical spacing to use between features.

    Returns
    -------
    pandas.DataFrame
        Annotated version of the input dataframe, containing an extra
        column 'y' that denotes the y-offset assigned to each feature.

    """

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

    levels, level_heights = _pack_ffdh(data.itertuples())

    level_offsets = np.cumsum(level_heights, dtype=np.float) - level_heights[0]
    level_offsets += np.arange(1, len(level_offsets) + 1) * spacing

    rows = []
    for objs, offset in zip(levels, level_offsets):
        for obj in objs:
            rows.append((obj.Index, offset))

    index, y_values = zip(*rows)
    heights = pd.DataFrame({'y': y_values}, index=index)

    return heights


def _pack_ffdh(objects):
    """Implementation of the First-Fit Decreasing Height packing algorithm.

    Packs objects into multiple levels, whilst trying to minimize the total
    used height. Uses a greedy approach that sorts objects by decreasing height,
    and then packs them accordingly.

    See http://cgi.csc.liv.ac.uk/~epa/surveyhtml.html for a description
    of the algorithm and other options for packing algorithms.

    Parameters
    ----------
    objects : List[Any]
        Objects to pack. Any object can be used, as long as the objects have
        'height', 'start' and 'end' properties, which define the height
        and range of the object.

    Returns
    -------
    Tuple[List[List[Any]], List[int]]
        Returns a tuple of level assignments (containing the objects) assigned
        to each level and a list containing the heights of each level.

    """

    # Sort boxes by decreasing height and size.
    objects = sorted(
        objects, key=lambda obj: (obj.height, obj.end - obj.start))
    objects = objects[::-1]

    levels = []
    level_heights = []

    for obj in objects:
        placed = False

        for level in levels:
            if not level.overlaps(obj.start, obj.end):
                level.addi(obj.start, obj.end, obj)
                placed = True
                break

        if not placed:
            new_level = IntervalTree.from_tuples([(obj.start, obj.end, obj)])
            levels.append(new_level)
            level_heights.append(obj.height)

    levels_list = [[interval[2] for interval in level] for level in levels]

    return levels_list, level_heights


def _augment_with_labels(data, label, label_func, ax):
    # Augment positions.
    anchor_col = 'end' if _reversed_axis(ax) else 'start'

    positions = [
        label_func(
            x=x, y=0, label=lab, ax=ax).x0
        for lab, x in zip(data[label], data[anchor_col])
    ]

    return data.assign(**{anchor_col: positions})


def _reversed_axis(ax):
    """Checks if x-axis is reversed."""
    xlim = ax.get_xlim()
    return xlim[1] < xlim[0]


class RugTrack(Track):
    """Track that plots density ticks for features.

    Parameters
    ----------
    data : pandas.Dataframe
        Dataset for plotting. Each row in the DataFrame is expected to
        correspond with a single feature. The DataFrame should have the
        following columns: seqname, position, strand; which together
        specify the location and orientation of the corresponding feature.
    hue : str
        Column (categorical) that should be used to determine the color
        of a given feature.
    hue_order : List[str]
        Order to plot the categorical hue levels in, otherwise the levels
        are inferred from the data objects.
    palette : List[Union[str, Tuple[float, float, float]]]
        Colors to use for the different levels of the hue variable.
        Should be specified as a list of colors (strs) or a list of
        tuples with RGB values (similar to Seaborn color palettes).
    height : float
        The height of the track.
    line_kws : dict[str, Any]
        Dict of keyword arguments to pass to LineCollection when drawing
        the ticks. Used to specify modify the aesthetics of the ticks.

    """

    def __init__(self,
                 data,
                 hue=None,
                 hue_order=None,
                 palette=None,
                 height=1.0,
                 line_kws=None):
        super().__init__()

        self._data = data
        self._height = height

        self._hue = hue
        self._color_map = build_colormap(
            data, hue=hue, palette=palette, order=hue_order)
        self._line_kws = line_kws or {}

    def get_height(self, region, ax):
        """Returns the height of the track.

        Parameters
        ----------
        region : Tuple[str, int, int]
            The genomic region that will be drawn. Specified as a tuple of
            (chromosome, start, end).
        ax : matplotlib.Axes
            Axis that the track will be drawn on.

        Returns
        -------
        height : int
            Height of the track.

        """
        return self._height

    def draw(self, region, ax):
        """Draws the track on the given axis.

        Parameters
        ----------
        region : Tuple[str, int, int]
            Genomic region to draw.
        ax : matplotlib.Axes
            Axis to draw track on.

        """

        data = self._data.query(
            'chromosome == {!r} and position > {} and position < {}'
            .format(*region))  # yapf: disable

        if self._hue is not None:
            for hue, grp in data.groupby(self._hue):
                self._draw_lines(grp, ax, color=self._color_map[hue])
        else:
            self._draw_lines(data, ax)

        ax.yaxis.set_visible(False)

    def _draw_lines(self, data, ax, color=None):
        segments = (((tup.position, 0), (tup.position, 1))
                    for tup in data.itertuples())
        lines = LineCollection(segments, color=color, **self._line_kws)
        ax.add_collection(lines)
