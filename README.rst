===============================
Geneviz
===============================

A Python library for creating visualizations of genomic data.

To improve interpretation, genomic data are frequently drawn together with a
variety of genomic annotation features using tools such as genome browsers.
Examples of such features include gene or transcript models from the reference
genome, CpG islands etc. Although many existing genome browsers are good at
visualizing genomic data, it is difficult or often impossible to interactively
use these tools from a scripting language such as Python. Geneviz aims to
bridge this gap by providing a basic framework for plotting various types
of data along genomic coordinates in Python.

Note that Geneviz is still under development. The API may therefore
likely be subject to change and documentation may be incomplete.

* Free software: MIT license
* Documentation: https://jrderuiter.github.io/geneviz.

Features
--------

* Plot various datatypes along a genomic axis.
* Annotate your data with genomic features such as genes/transcripts.
* Query genomic features from a variety of sources (GTF, Biomart, etc.).
* Visualize different types of NGS data (coverage/splice junctions).

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

