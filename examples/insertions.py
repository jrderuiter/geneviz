import pandas as pd
import seaborn as sns
from geneviz.tracks import plot_tracks, FeatureTrack, BiomartTrack

sns.set_style('white')

features = pd.DataFrame.from_records(
    [('1', 182448297, 1), ('1', 182448100, 1), ('1', 182419592, -1)],
    columns=['chromosome', 'position', 'strand'])

feature_track = FeatureTrack.from_position(
    features, width=3500, height=0.5, hue='strand')

gene_track = BiomartTrack(
    collapse='transcript',
    dataset='mmusculus_gene_ensembl',
    gene_id='gene_name',
    height=0.5)

plot_tracks(
    [feature_track, gene_track],
    region=('1', 182399172, 182462432),
    figsize=(10, None),
    despine=True)
