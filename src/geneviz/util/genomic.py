import operator

import itertools
import operator
from typing import Any, Iterable

import numpy as np
import pandas as pd
from intervaltree import IntervalTree


class GenomicDataFrame(pd.DataFrame):
    """DataFrame with fast indexing by genomic position.

    Requires columns 'chromosome', 'start' and 'end' to be present in the
    DataFrame, as these columns are used for indexing.
    """

    def __init__(self,
                 *args,
                 chrom_col='chromosome',
                 start_col='start',
                 end_col='end',
                 **kwargs):
        super().__init__(*args, **kwargs)

        self._chrom_col = chrom_col
        self._start_col = start_col
        self._end_col = end_col
        self._trees = None

    @property
    def _constructor(self):
        return GenomicDataFrame

    @property
    def trees(self):
        """Returns trees used for indexing the DataFrame."""
        if self._trees is None:
            self._trees = self._build_trees()
        return self._trees

    def _build_trees(self):
        # Subset frame to positions and rename columns to defaults.
        position_df = self[[self._chrom_col, self._start_col, self._end_col]]
        position_df = position_df.rename(columns={
            self._chrom_col: 'chromosome',
            self._start_col: 'start',
            self._end_col: 'end'
        })

        # Add index and sort by chromosome (for grouping).
        position_df = position_df.assign(index=np.arange(len(self)))
        position_df = position_df.sort_values(by='chromosome')

        # Convert to tuples.
        tuples = ((tup.chromosome, tup.start, tup.end, tup.index)
                  for tup in position_df.itertuples())

        return GenomicIntervalTree.from_tuples(tuples)

    def search(self, chromosome, begin, end):
        """Subsets the DataFrame for rows within given range."""
        overlap = self.trees.search(chromosome, begin, end)
        indices = [interval[2] for interval in overlap]
        return self.iloc[indices].sort_index()


class GenomicIntervalTree(object):
    """Datastructure for efficiently accessing genomic objects by position."""

    def __init__(self, trees):
        # type: (Dict[str, IntervalTree]) -> None
        self._trees = trees

    def __getitem__(self, i):
        # type: (str) -> IntervalTree
        """Returns tree with given chromosome name."""
        return self._trees[i]

    @classmethod
    def from_tuples(cls, tuples):
        """Builds an instance from tuples.

        Assumes tuples are sorted by chromosome.
        """

        # Group by chromosome.
        groups = itertools.groupby(tuples, key=operator.itemgetter(0))

        # Build trees.
        trees = {}
        for chrom, group in groups:
            trees[chrom] = IntervalTree.from_tuples(
                (start, end, obj) for _, start, end, obj in group)

        return cls(trees)

    def search(self, chromosome, begin, end=None):
        # type: (str, int, int) -> Iterable[object]
        """Searches the tree for objects within given range."""
        return self._trees[chromosome].search(begin, end)


def merge_intervals(intervals):
    """Merges overlapping intervals.

    Parameters:
        intervals:
        is_sorted:

    Return:
        Generator containing merged intervals.

    """

    # Sort intervals.
    intervals = sorted(intervals, key=operator.itemgetter(0))

    if not intervals:
        # No intervals to merge.
        return

    # Low and high represent the bounds of
    # the current run of merges.
    low, high = intervals[0]

    for iv in intervals[1:]:
        if iv[0] <= high:
            # New interval overlaps current run,
            # merge with the current run.
            high = max(high, iv[1])
        else:
            # Current run is over, yield accumulated
            # interval and start new run.
            yield low, high
            low, high = iv

    yield low, high
