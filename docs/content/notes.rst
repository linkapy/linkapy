Notes
-----

Input files
===========

A couple of things need to be kept in mind with regards to methylation data.
At this point, five different file types for methylation data are supported:

- Allcools files
- MethylDackel bedgraph files
- Bismark coverage files
- Bismark CpG report files
- BedMethyl files

Note that MethylDackel bedgraph files, Bismark CpG report files and BedMethyl files are assumed to be 0-based start encoded (and 1-based end, if applicable).
The Allcools files and Bismark coverage files are assumed to be 1-based encoded. 
For the Bismark coverage files, keep in mind that `bismark_methylation_extractor` has a flag to output 0-based files, so pay attention that this is correct.
For BedMethyl files, please note that _no_ checks are performed to ensure that only one modification type is present.
If you have multiple modification types (i.e. more then one 'name' in column 4), please split htem into separate files (you can include them in the output by specific multiple methylation_path/methylation_pattern combinations).

For the RNA-part, for now only featureCounts tables are supported.

Region, blacklist and chromsizes files
=======================================

``--regions`` and ``--blacklist`` expect tab-separated, BED-like files (gzip-compressed files are also accepted) with at least 3 columns: ``chrom``, ``start``, ``end``.
A 4th column, if present, is used as the region name; otherwise the name defaults to ``<chrom>:<start>-<end>``. No header line should be present.

``--chromsizes`` expects a tab-separated, headerless file with 2 columns: ``chrom`` and ``size`` (in bp). It is only required when no ``--regions`` file is given, in which case the genome is tiled into consecutive bins of ``--binsize`` bp per chromosome.

Combining multiple featureCounts tables
========================================

When multiple files are matched by a single ``--transcriptome_pattern``, they are combined into a single count matrix. This requires all files to share identical ``Geneid``, ``Chr``, ``Start``, ``End``, ``Strand`` and ``Length`` columns (i.e. they were generated against the same annotation/feature set); an assertion error is raised otherwise. Sample names are derived from each count column's header, taking the substring before the first ``.``.

Modality naming and cell matching
==================================

Each pattern (methylation or transcriptome) yields its own AnnData object, stored in the final MuData object under the key ``METH_<pattern>`` or ``RNA_<pattern>``. When multiple modalities are produced, Linkapy attempts to match cell barcodes across them; the applied renaming (if any) is written to ``cell_renaming.tsv`` in the output directory. See :doc:`usage` for details on the output layout and how to load the resulting ``.h5mu`` file.