"""This module provides basic tracks and supporting actor classes
for plotting generic genomic features. The module provides the
following tracks:

    - StackedTrack: A generic Track that takes any objects
        implementing the Actor interface and plots these
        in a stacked fashion along the genomic axis.
    - FeatureTrack:
        A more specific track that draws (un)directed
        Feature actors on the genomic axis. The track can
        optionally be generated from a dataframe containing
        genomic features.

"""

from __future__ import absolute_import, division, print_function

#pylint: disable=W0622,W0611
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from operator import attrgetter
from collections import namedtuple

from intervaltree import IntervalTree
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
import numpy as np
import toolz

from geneviz.tracks import Track


_BoundingBox = namedtuple('BoundingBox', ['start', 'end', 'height', 'actor'])


class BoundingBox(_BoundingBox):
    """BoundingBox representing the extents of actor classes.

    The attributes start and end respectively represent the genomic start and
    end location of an actor. The height attribute represents the maximum
    height of the actor on the axis (in data coordinates).

    """


class Actor(object):
    """Abstract base class for objects drawn on a StackedTrack.

    Specifies two methods (get_bbox and draw) that should be overridden
    in each Actor subclass. `get_bbox` is used to determine the bounds
    of each actor, which are taken into account when stacking actors
    within the plox. `draw `is called by the StackedTrack to actually
    draw the respective actors.

    The seqname, start, end attributes of an actor are used to determine
    whether an actor falls within the plotted genomic range. `get_bbox`
    is used when actually drawing an actor to identify overlaps between
    different actors (and elements such as their labels).

    Attributes:
        seqname (str): Chromosome of the actor.
        start (int): Genomic start of the actor.
        end (int): Genomic end of the actor.

    """

    def __init__(self, seqname, start, end):
        self.seqname = seqname
        self.start = start
        self.end = end

    def get_height(self):
        """Height of the actor."""
        raise NotImplementedError()

    def get_bounding_box(self, ax):
        """Calculates bounding box on the given axis."""
        raise NotImplementedError()

    def draw(self, ax, y=1):
        """Draws the actor on the given axis.

        Args:
            ax (matplotlib.Axes): Axis on which actor should be drawn.
            y (int): Height at which to draw the actor.

        """
        raise NotImplementedError()


class StackedTrack(Track):
    """Track for plotting actors in a stacked fashion.

    Track class that accepts any generic actor class and plots these
    together on an axis. Actors are stacked if they have overlapping
    positions.

    Args:
        actors (List[Actor]): List of actors to draw.
        reverse (bool): Wether to reverse the y-axis.
        scale (float): Relative scale at which the y-axis is drawn.

    """

    def __init__(self, actors, reverse=False, spacing=0.2):
        super().__init__()
        self._actors = actors
        self._reverse = reverse
        self._spacing = spacing

    def _stack(self, actors, ax):
        """Determines stacked positions for given actors.

        Arguments:
            actors (List[Actor]): List of actors to stack.
            ax (matplotlib.Axes): Axis to be used in bounding box
                calculation for the actors (if needed).

        """

        bounding_boxes = [a.get_bounding_box(ax) for a in actors]
        levels, level_heights = pack_ffdh(bounding_boxes)

        level_cumsum = np.hstack([0, np.cumsum(level_heights)[:-1]])
        level_spacing = (np.arange(len(level_cumsum)) + 1) * self._spacing

        level_offsets = level_cumsum + level_spacing

        for level, offset in zip(levels, level_offsets):
            for box in level:
                yield box.actor, offset

    def _get_actors(self, seqname, start, end):
        return (a for a in self._actors
                if a.seqname == seqname and
                a.end > start and a.start < end)

    def draw(self, ax, seqname, start, end):
        actors = self._get_actors(seqname, start, end)
        stacked = self._stack(actors, ax=ax)

        y_max = 0
        for actor, y in stacked:
            actor.draw(ax, y=y)
            y_max = max(y_max, y + actor.get_height())

        y_max += self._spacing

        # Set ylim.
        if not self._reverse:
            ax.set_ylim(0, y_max)
        else:
            ax.set_ylim(y_max, 0)

        # Style axes.
        ax.yaxis.set_visible(False)

    def get_height(self, ax, seqname, start, end):
        actors = self._get_actors(seqname, start, end)
        bounding_boxes = (a.get_bounding_box(ax) for a in actors)

        _, level_heights = pack_ffdh(bounding_boxes)
        total_spacing = (len(level_heights) + 1) * self._spacing

        return sum(level_heights) + total_spacing


def pack_ffdh(bounding_boxes):
    """Implementation of the First-Fit Decreasing Height algorithm.

    See http://cgi.csc.liv.ac.uk/~epa/surveyhtml.html for a description
    of the algorithm and other options for packing algorithms.

    """

    height_sorted_boxes = sorted(bounding_boxes,
                                 key=attrgetter('height'))[::-1]

    levels = []
    level_heights = []

    for box in height_sorted_boxes:
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

    levels_list = [[interval[2] for interval in level]
                   for level in levels]

    return levels_list, level_heights
