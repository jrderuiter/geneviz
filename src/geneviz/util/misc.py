from __future__ import absolute_import, division, print_function

#pylint: disable=W0622,W0611
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import operator


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


def get_region(data, seqname, start, end):
    """Returns region from data frame with seqname, start and end columns."""
    return data.query('seqname == {!r} and end >= {} and start <= {}'
                      .format(seqname, start, end))
