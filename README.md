[![documentation](https://github.com/WardDeb/linkapy/actions/workflows/docs.yml/badge.svg)](https://github.com/WardDeb/linkapy/actions/workflows/docs.yml)
[![installation](https://github.com/WardDeb/linkapy/actions/workflows/pip.yml/badge.svg)](https://github.com/WardDeb/linkapy/actions/workflows/pip.yml)
[![lint](https://github.com/WardDeb/linkapy/actions/workflows/lint.yml/badge.svg)](https://github.com/WardDeb/linkapy/actions/workflows/lint.yml)
[![PyPI version](https://img.shields.io/pypi/v/linkapy)](https://pypi.org/project/linkapy/)
[![real data test](https://github.com/FunctionalEpigeneticsLab/linkapy/actions/workflows/test_example.yml/badge.svg)](https://github.com/FunctionalEpigeneticsLab/linkapy/actions/workflows/test_example.yml)
[![pytests](https://github.com/FunctionalEpigeneticsLab/linkapy/actions/workflows/test_python.yml/badge.svg)](https://github.com/FunctionalEpigeneticsLab/linkapy/actions/workflows/test_python.yml)
[![rusttests](https://github.com/FunctionalEpigeneticsLab/linkapy/actions/workflows/test_rust.yml/badge.svg)](https://github.com/FunctionalEpigeneticsLab/linkapy/actions/workflows/test_rust.yml)
![badge](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/WardDeb/d0010eb142b962632f94c164c502b506/raw/coverage.json)
![badge](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/WardDeb/f4d532defe4f2caecec457d6653d933e/raw/coverage.json)
[![DOI](https://zenodo.org/badge/893384655.svg)](https://doi.org/10.5281/zenodo.23000297)


# linkapy
A framework to process and analyse (multi-modal) methylation  single-cell data.
Linkapy is mainly used to process single-cell methylation data, and is focused on creating an integrated (muData) object of multimodal single-cell data. In this context it means that multiple modalities (e.g. accessibility, DNA methylation, RNA expression) are available from the same cell. Note that if only methylation data is available, this package can still be used to process and analyse it.

Methylation input can be provided as allcools, MethylDackel bedgraph, Bismark coverage, Bismark CpG report, or BedMethyl files, including NOMe-seq data (`--NOMe`). Transcriptome input is provided as featureCounts tables. See [notes](docs/content/notes.rst) for format details, and [usage](docs/content/usage.rst) for the full CLI reference and output description.

# Documentation

Note that this package is still in development and is prone to change without notice. 
A rudimentary [documentation](https://linkapy.readthedocs.io/en/latest/) is available.

# Quickstart

  > pip install linkapy  
  > linkapy --help
