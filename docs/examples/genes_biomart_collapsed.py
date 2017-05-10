from matplotlib import pyplot as plt
import seaborn as sns

from geneviz.tracks import BiomartTrack, plot_tracks

sns.set_style('white')

gene_track = BiomartTrack(
    dataset='mmusculus_gene_ensembl',
    hue='strand',
    gene_id='gene_name',
    transcript_id='transcript_id',
    height=0.5,
    spacing=0.1,
    collapse='transcript')

plot_tracks(
    [gene_track],
    region=('15', 77618614, 77843156),
    despine=True,
    figsize=(8, None))
