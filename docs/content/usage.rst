Usage
-----

Getting started
~~~~~~~~~~~~~~~

Linkapy is a Python package that is designed to facilitate the integrative analysis of single-cell multi-omics data, where multi-omics means multiple read outs of the same cell.
While an attempt is made to keep Linkapy as general as possible, for now it primarily focuses on data that includes methylation and transcription layers.
If available, transcriptome data should come as one (or more) featureCount tables, which will combined into a single matrix.
Methylation data can be provided in several formats (allcools, MethylDackel bedgraph, Bismark coverage, Bismark CpG report, or BedMethyl); see :doc:`notes` for format-specific details.

Usage
~~~~~

Linkapy can be used through the command line, or via the API in Python. For the latter, have a look at the `API Reference <../autoapi/index.html>`_.
To get started, example data from the original `scNMT-seq paper <https://www.nature.com/articles/s41467-018-03149-4>`_ can be downloaded:

.. code-block:: console

    linkapy example -h

Upon successfull download, an example command will be printed that you can use to get started and familiarize yourself with the data structures.

Output
~~~~~~

Running ``linkapy parsing`` populates the output directory (``--output``) with:

- ``matrices/`` - intermediate per-pattern matrices (arrow format for RNA, mtx format for methylation) that feed into the final object.
- ``<project>.h5mu`` - the final integrated object, in `MuData <https://mudata.readthedocs.io/>`_ format.
- ``<project>.log`` - a log of the run.
- ``cell_renaming.tsv`` - written only when multiple modalities are present; records how cell barcodes were matched/renamed across modalities.

The ``<project>.h5mu`` file is a **MuData** object, not a plain AnnData ``.h5ad`` file. Load it with the ``mudata`` package rather than ``anndata``:

.. code-block:: python

    import mudata as md
    mdata = md.read_h5mu("<project>.h5mu")

Each modality is stored under ``mdata.mod`` using the key ``RNA_<pattern>`` or ``METH_<pattern>``, where ``<pattern>`` is the (name-resolved) value passed to ``--transcriptome_pattern``/``--methylation_pattern`` (or their ``_names`` counterparts, or the two NOMe-derived patterns ``GCHN``/``WCGN`` when ``--NOMe`` is set). For example:

.. code-block:: python

    mdata.mod["METH_WCGN"]  # AnnData for the WCGN (CpG) NOMe modality
    mdata.mod["RNA_txt"]    # AnnData for transcriptome data matched on pattern "txt"

Each modality's ``.var`` index is additionally prefixed with its pattern (e.g. ``METH_WCGN:<region>``) to keep variable names unique across modalities.

Methylation modalities additionally carry per-cell QC in ``.obs`` and per-region QC in ``.var``, computed during aggregation; see :doc:`notes` for what each column means.

.. click:: linkapy.CLI:linkapy
   :prog: linkapy CLI
   :nested: full
