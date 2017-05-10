# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

from geneviz.util import misc

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


class TestMergeIntervals(object):

    def test_basic(self):
        intervals = [(10, 20), (30, 40), (15, 25)]
        merged = list(misc.merge_intervals(intervals))
        assert merged == [(10, 25), (30, 40)]

    def test_non_overlapping(self):
        intervals = [(10, 20), (30, 40), (5, 9)]
        merged = list(misc.merge_intervals(intervals))
        assert merged == [(5, 9), (10, 20), (30, 40)]

    def test_empty(self):
        merged = list(misc.merge_intervals([]))
        assert merged == []
