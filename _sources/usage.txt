=====
Usage
=====

Drawing genomic features
------------------------

.. code-block:: python

    # Example features.
    features = pd.DataFrame.from_records([
            ('1', 10, 25, 1, 'a', 'a'),
            ('1', 30, 40, None, 'a', 'b'),
            ('1', 50, 70, None, 'a', 'c')
        ],
        columns=['chromosome', 'start', 'end',
                 'strand', 'group', 'type'])

    # Plot track.
    feat_track = FeatureTrack(data=features, hue='type')
    plot_tracks([feat_track], region=('1', 0, 80),
                figsize=(5, None), despine=True)


.. plot:: examples/features.py

.. code-block:: python

    FeatureTrack(data=features, hue='type', label='type')

.. plot:: examples/features_hue.py

.. code-block:: python

    FeatureTrack(data=features, hue='type', group='group', label='group')

.. plot:: examples/features_grouped.py

.. code-block:: python

    FeatureTrack(data=features, hue='type', group='group',
                 label='group', strand_junctions=True)

.. plot:: examples/features_grouped_stranded.py


Drawing genes/transcripts
-------------------------

.. code-block:: python

    from geneviz.tracks import BiomartTrack, plot_tracks

    gene_track = BiomartTrack(
        dataset='mmusculus_gene_ensembl',
        hue='strand',
        gene_id='gene_name',
        transcript_id='transcript_id',
        height=0.5,
        spacing=0.1)

    plot_tracks(
        [gene_track],
        region=('15', 77618614, 77843156),
        despine=True,
        figsize=(8, None))

.. plot:: examples/genes_biomart.py

.. code-block:: python

    BiomartTrack(
        dataset='mmusculus_gene_ensembl',
        hue='strand',
        gene_id='gene_name',
        transcript_id='transcript_id',
        height=0.5,
        spacing=0.1,
        collapse='transcript')

.. plot:: examples/genes_biomart_collapsed.py

.. code-block:: python

    BiomartTrack(
        dataset='mmusculus_gene_ensembl',
        hue='strand',
        gene_id='gene_name',
        transcript_id='transcript_id',
        height=0.5,
        spacing=0.1,
        collapse='gene')

.. plot:: examples/genes_biomart_gene.py


Combining tracks
----------------
