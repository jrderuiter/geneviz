import pandas as pd
from matplotlib import pyplot as plt

from geneviz.tracks import plot_tracks
from geneviz.tracks.feature import FeatureTrack

import seaborn as sns
sns.set_style('white')

# Setup features.
features = pd.DataFrame.from_records([
        ('1', 10, 25, 1, 'a', 'a'),
        ('1', 30, 40, None, 'a', 'b'),
        ('1', 50, 70, None, 'a', 'c')
    ], columns=['chromosome', 'start', 'end',
                'strand', 'group', 'type'])  # yapf: disable

# Plot un-grouped.
feat_track = FeatureTrack(data=features, hue='type', label='type')
plot_tracks([feat_track], region=('1', 0, 80), figsize=(5, None), despine=True)

plt.tight_layout()
