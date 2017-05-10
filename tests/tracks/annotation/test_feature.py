# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches

import pandas as pd
import seaborn as sns

import pytest

from geneviz.tracks.annotation import feature as ftrack

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


@pytest.fixture
def feature_kws():
    return dict(seqname='1', start=200, end=300,
                strand=-1, name='test', height=2)


class TestFeature(object):

    def test_construction(self, feature_kws):
        """Tests construction of a basic feature object."""

        feature = ftrack.Feature(**feature_kws)

        assert feature.seqname == feature_kws['seqname']
        assert feature.start == feature_kws['start']
        assert feature.end == feature_kws['end']
        assert feature.strand == feature_kws['strand']
        assert feature.name == feature_kws['name']
        assert feature.get_height() == feature_kws['height']

    def test_bounding_box(self, feature_kws):
        """Tests bounding box for feature without a label."""

        _, ax = plt.subplots()
        ax.set_xlim(150, 350)

        feature_kws.update({'name': None})
        feature = ftrack.Feature(**feature_kws)

        bbox = feature.get_bounding_box(ax=ax)
        assert bbox.start == feature_kws['start']
        assert bbox.end == feature_kws['end']
        assert bbox.height == feature_kws['height']
        assert bbox.actor == feature

    def test_bounding_box_label(self, feature_kws):
        """Tests bounding box for feature with a label."""

        # TODO: devise a better test for the start position.
        # Position seems to be fixed now, but this may be due to
        # used backend which isn't actually rendering in a window.

        _, ax = plt.subplots()
        ax.set_xlim(150, 350)

        feature = ftrack.Feature(**feature_kws)

        bbox = feature.get_bounding_box(ax=ax)
        assert bbox.start < feature_kws['start']
        assert bbox.end == feature_kws['end']
        assert bbox.height == feature_kws['height']

    def test_draw(self, feature_kws):
        """Tests feature draw."""

        _, ax = plt.subplots()

        feature = ftrack.Feature(**feature_kws)
        feature.draw(ax)

        assert len(ax.patches) == 1
        assert len(ax.texts) == 1
        assert isinstance(ax.patches[0], mpl_patches.FancyArrow)

    def test_draw_unstranded(self, feature_kws):
        """Tests unstranded feature draw."""

        _, ax = plt.subplots()

        feature_kws.update({'strand': None})
        feature = ftrack.Feature(**feature_kws)
        feature.draw(ax)

        assert len(ax.patches) == 1
        assert len(ax.texts) == 1
        assert isinstance(ax.patches[0], mpl_patches.Rectangle)

    def test_draw_without_label(self, feature_kws):
        _, ax = plt.subplots()

        feature_kws.update({'name': None})
        feature = ftrack.Feature(**feature_kws)
        feature.draw(ax)

        assert len(ax.patches) == 1
        assert len(ax.texts) == 0


@pytest.fixture
def test_data():
    return pd.DataFrame.from_records(
        [('1', 10, 20, -1, 'feat_a'),
         ('1', 30, 40, 1, 'feat_b'),
         ('2', 10, 20, -1, 'feat_c')],
        columns=['seqname', 'start', 'end', 'strand', 'name'])


@pytest.fixture
def test_region():
    return {'seqname': '1', 'start': 0, 'end': 50}


class TestFeatureTrack(object):

    def test_features(self, test_data, test_region):
        """Test features from track."""

        track = ftrack.FeatureTrack(test_data)
        features = list(track._get_features(**test_region))
        assert len(features) == 2

    def test_get_height(self, test_data, test_region):
        """Test get_height from track."""

        _, ax = plt.subplots()

        track = ftrack.FeatureTrack(test_data)
        assert track.get_height(ax, **test_region) == 1.4

    def test_draw(self, test_data, test_region):
        """Test basic draw."""

        _, ax = plt.subplots()

        track = ftrack.FeatureTrack(test_data)
        track.draw(ax, **test_region)

        assert len(ax.patches) == 2
        assert len(ax.texts) == 0

    def test_draw_labels(self, test_data, test_region):
        """Test basic draw with labels."""

        _, ax = plt.subplots()

        track = ftrack.FeatureTrack(test_data, name='name')
        track.draw(ax, **test_region)

        assert len(ax.patches) == 2
        assert len(ax.texts) == 2

    def test_hue(self, test_data, test_region):
        """Test applying hue to features."""

        palette = sns.color_palette()
        track = ftrack.FeatureTrack(test_data, name='name',
                                    hue='strand', palette=palette)

        _, ax = plt.subplots()
        track.draw(ax, **test_region)

        assert ax.patches[0].get_facecolor()[:3] in palette
        assert ax.patches[1].get_facecolor()[:3] in palette

    def test_hue_dict(self, test_data, test_region):
        """Test applying hue (as dict) to features."""

        palette = {1: sns.color_palette()[0],
                   -1: sns.color_palette()[1]}

        track = ftrack.FeatureTrack(test_data, name='name',
                                    hue='strand', palette=palette)

        _, ax = plt.subplots()
        track.draw(ax, **test_region)

        assert ax.patches[0].get_facecolor()[:3] == palette[1]
        assert ax.patches[1].get_facecolor()[:3] == palette[-1]
