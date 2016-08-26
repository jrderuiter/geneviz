from __future__ import absolute_import, division, print_function

#pylint: disable=W0622,W0611
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import itertools

import seaborn as sns
import toolz

from geneviz.tracks import Track


class DataTrack(Track):

    _default_kws = dict(marker='o')

    def __init__(self, data, y=None, hue=None, palette=None,
                 height=1, legend_kws=None, plot_kws=None):
        super().__init__(height=height)

        y = y or 'y'
        palette = palette or sns.color_palette()

        self._data = self._preprocess_data(data, y, hue, palette)

        self._legend_kws = legend_kws or {}
        self._plot_kws = toolz.merge(self._default_kws, plot_kws or {})

    @staticmethod
    def _preprocess_data(data, y, hue, palette):
        plot_data = data[['seqname', 'position']].copy()
        plot_data['y'] = data[y]

        if hue is not None:
            if not isinstance(palette, dict):
                palette = dict(zip(data[hue].unique(),
                                   itertools.cycle(palette)))
            plot_data['hue'] = data[hue]
            plot_data['color'] = data[hue].map(palette)

        return plot_data

    def draw(self, ax, seqname, start, end):
        # Subset data for range.
        data = self._data.query(('seqname == {!r} and '
                                 '{} <= position <= {}')
                                .format(seqname, start, end))

        # Gather plotting kwargs.
        if 'hue' in data and 'color' in data:
            for (hue, color), grp in data.groupby(['hue', 'color']):
                grp = grp.sort_values(by='position')
                ax.plot(grp['position'], grp['y'], '.', label=hue,
                        color=color, **self._plot_kws)
            ax.legend(**self._legend_kws)

        else:
            data = data.sort_values(by='position')
            ax.plot(data['position'], data['y'], **self._plot_kws)


class BoxplotTrack(Track):

    # TODO: Seaborn styling?

    def __init__(self, data, seqname='seqname', x='x', y='y',
                 height=1, **kwargs):
        super().__init__(height=height)

        self._data = data

        self._seqname = seqname
        self._x = x
        self._y = y

        self._kwargs = kwargs

    def draw(self, ax, seqname, start, end):
        # Subset data for range.
        mask = ((self._data[self._seqname] == seqname) &
                (self._data[self._x] >= start) &
                (self._data[self._x] <= end))
        data = self._data[mask]

        grouped = ((pos, grp[self._y].values)
                   for pos, grp in data.groupby(self._x))
        positions, values = zip(*grouped)

        ax.boxplot(values, positions=positions, **self._kwargs)
