from __future__ import absolute_import, division, print_function

#pylint: disable=W0622,W0611
from builtins import (ascii, bytes, chr, dict, filter, hex, input, int, map,
                      next, oct, open, pow, range, round, str, super, zip)
from future.utils import native_str

import numpy as np
import pysam
import toolz

from matplotlib import (path as mpath, collections as mcollections, patches as
                        mpatches)

from geneviz.tracks import Track
from geneviz.util.tabix import BedIterator


class SpliceTrack(Track):

    _default_kws = dict(facecolor=None, lw=1)

    def __init__(self, data, height=1, color=None, patch_kws=None):
        super().__init__()

        self._data = data
        self._height = height

        patch_kws = toolz.merge(self._default_kws, {'edgecolor': color},
                                patch_kws or {})
        self._patch_kws = patch_kws

    def get_height(self, region, ax):
        return self._height

    def _fetch_data(self, region):
        return self._data.query(
            ('chromosome == {!r} and end >= {} and start <= {}')
            .format(*region))  # yapf: disable

    def draw(self, region, ax):
        data = self._fetch_data(region)

        arcs = (self._splice_arc(start, end, score)
                for start, end, score in
                zip(data['start'], data['end'], data['score'])) # yapf: disable

        patches = mcollections.PatchCollection(arcs, **self._patch_kws)

        ax.add_collection(patches)

        ax.set_ylim(0, data['score'].max())

    @staticmethod
    def _splice_arc(start, end, height):
        bezier_height = height / 0.75
        vertices = [(start, 0.), (start, bezier_height), (end, bezier_height),
                    (end, 0)]

        codes = [
            mpath.Path.MOVETO, mpath.Path.CURVE4, mpath.Path.CURVE4,
            mpath.Path.CURVE4
        ]

        path = mpath.Path(vertices, codes)
        patch = mpatches.PathPatch(path)

        return patch


class CoverageTrack(Track):
    def __init__(self,
                 bam_path,
                 height=1,
                 fill=True,
                 stepper='all',
                 plot_kwargs=None,
                 resample_interval=None):
        super().__init__()

        # Bam file parameters.
        self._bam_path = bam_path
        self._stepper = stepper

        # Draw parameters.
        self._height = height
        self._fill = fill
        self._resample_interval = resample_interval

        self._plot_kwargs = {} if plot_kwargs is None else plot_kwargs

    def get_height(self, region, ax):
        return self._height

    def _get_coverage(self, region):
        seqname, start, end = region
        hist = np.zeros(end - start)

        with pysam.AlignmentFile(native_str(self._bam_path), 'rb') as file_:
            # Truncate = True truncates pileups at start-end
            # positions, avoids looking outside range.
            pileups = file_.pileup(
                reference=seqname,
                start=start,
                end=end,
                stepper=self._stepper,
                truncate=True)

            for pileup in pileups:
                # Count reads at each position.
                coverage = 0
                for read in pileup.pileups:
                    # Is_del == 0 if read is positioned here, i.e.,
                    # this ignores reads spliced over position.
                    if read.is_del == 0:
                        coverage += 1

                hist[pileup.pos - start] = coverage

        return hist

    def draw(self, region, ax):
        _, start, end = region

        # Determine coverage.
        x_range = np.arange(start, end)
        coverage = self._get_coverage(region)

        # Resample if needed.
        if self._resample_interval is not None:
            x_new = np.arange(start, end, step=self._resample_interval)
            coverage = np.interp(x_new, xp=x_range, fp=coverage)
            x_range = x_new

        # Plot coverage line.
        ax.plot(x_range, coverage, **self._plot_kwargs)

        if self._fill:
            ax.fill_between(x_range, 0, coverage)
