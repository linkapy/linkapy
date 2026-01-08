import pytest
import gzip
import numpy as np
import logging

'''
Create fixtures for testing.
The genome here is considered to be 1 chromosome (chr1) with 500 bps in total.
There are 5 methylation (WCGN context) sites:
 - pos 20 TCGG
 - pos 50 ACGT
 - pos 221 ACGA
 - pos 244 ACGA
 - pos 440 ACGA
and 4 methylation (GCHN context) sites:
 - pos 30 GCAA
 - pos 40 GCTC
 - pos 230 GCCC
 - pos 401 GCAG
There are 2 genes in total:
 - gene1: 0-100 bps
 - gene2: 200-250 bps

Also includes some helper function to parse mudata objects.
'''

@pytest.fixture
def dynamic_methylation_path(request):
    '''
    helps to return fixture to parametrize over.
    '''
    return request.getfixturevalue(request.param)

@pytest.fixture
def allcools_path(tmp_path):
    '''
    simulate allcools.tsv.gz files
    tsv with columns chr, pos, strand, context, meth, cov and sign.
    sign denotes if significantly methylated (1 = no test).
    1-based
    '''
    methpath = tmp_path / "allcools"
    methpath.mkdir()

    # Cell 1 - avg. meth in gene 1 = 50%
    with gzip.open(methpath / "cell1.WCGN.tsv.gz", "wt") as f:
        f.write('chr1\t21\t+\tTCGG\t25\t100\t1\n')
        f.write('chr1\t51\t+\tACGT\t75\t100\t1\n')

    # Cell 2 - avg. meth in gene 2 = 62.5%
    with gzip.open(methpath / "cell2.WCGN.tsv.gz", "wt") as f:
        f.write('chr1\t222\t+\tACGA\t2\t4\t1\n')
        f.write('chr1\t245\t+\tACGA\t12\t16\t1\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.WCGN.tsv.gz", "wt") as f:
        f.write('chr1\t21\t+\tTCGG\t21\t21\t1\n')
        f.write('chr1\t51\t+\tACGT\t35\t35\t1\n')
        f.write('chr1\t222\t+\tACGA\t0\t100\t1\n')
    
    # Cell 1 - avg. meth in gene 1 = 25%
    with gzip.open(methpath / "cell1.GCHN.tsv.gz", "wt") as f:
        f.write('chr1\t31\t+\tGCAA\t25\t100\t1\n')
        f.write('chr1\t41\t+\tGCTC\t1\t4\t1\n')

    # Cell 2 - avg. meth in gene 2 = 44%
    with gzip.open(methpath / "cell2.GCHN.tsv.gz", "wt") as f:
        f.write('chr1\t231\t+\tGCCC\t44\t100\t1\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.GCHN.tsv.gz", "wt") as f:
        f.write('chr1\t31\t+\tGCAA\t100\t100\t1\n')
        f.write('chr1\t41\t+\tGCTC\t4\t4\t1\n')
        f.write('chr1\t231\t+\tGCCC\t0\t100\t1\n')

    return methpath

@pytest.fixture
def methyldackel_path(tmp_path):
    '''
    simulate methyldackel files
    tsv with columns chrom, start, end, perc, methylcount, coverage
    0-based
    '''
    methpath = tmp_path / "methyldackel"
    methpath.mkdir()

    # Cell 1 - avg. meth in gene 1 = 50%
    with open(methpath / "cell1.WCGN.bedgraph", "w") as f:
        f.write('chr1\t20\t21\t25\t25\t100\n')
        f.write('chr1\t50\t51\t75\t75\t100\n')

    # Cell 2 - avg. meth in gene 2 = 62.5%
    with open(methpath / "cell2.WCGN.bedgraph", "w") as f:
        f.write('chr1\t221\t222\t50\t2\t4\n')
        f.write('chr1\t244\t245\t75\t12\t16\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with open(methpath / "cell3.WCGN.bedgraph", "w") as f:
        f.write('chr1\t20\t21\t100\t21\t21\n')
        f.write('chr1\t50\t51\t100\t35\t35\n')
        f.write('chr1\t221\t222\t0\t0\t100\n')
    
    # Cell 1 - avg. meth in gene 1 = 25%
    with open(methpath / "cell1.GCHN.bedgraph", "w") as f:
        f.write('chr1\t30\t31\t25\t25\t100\n')
        f.write('chr1\t40\t41\t25\t1\t4\n')

    # Cell 2 - avg. meth in gene 2 = 44%
    with open(methpath / "cell2.GCHN.bedgraph", "w") as f:
        f.write('chr1\t230\t231\t44\t44\t100\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with open(methpath / "cell3.GCHN.bedgraph", "w") as f:
        f.write('chr1\t30\t31\t100\t100\t100\n')
        f.write('chr1\t40\t41\t100\t4\t4\n')
        f.write('chr1\t230\t231\t0\t0\t100\n')
    return methpath

@pytest.fixture
def bismarkcov_path(tmp_path):
    '''
    simulate bismark coverage files
    tsv with columns chrom, start, end, perc, methylcount, non methylcount
    these are 1-based.
    '''
    methpath = tmp_path / "bismarkcov"
    methpath.mkdir()

    # Cell 1 - avg. meth in gene 1 = 50%
    with gzip.open(methpath / "cell1.WCGN.bismark.cov", "wt") as f:
        f.write('chr1\t21\t22\t25\t25\t75\n')
        f.write('chr1\t51\t52\t75\t75\t25\n')

    # Cell 2 - avg. meth in gene 2 = 62.5%
    with gzip.open(methpath / "cell2.WCGN.bismark.cov", "wt") as f:
        f.write('chr1\t222\t223\t50\t2\t2\n')
        f.write('chr1\t245\t246\t75\t12\t4\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.WCGN.bismark.cov", "wt") as f:
        f.write('chr1\t21\t22\t100\t21\t0\n')
        f.write('chr1\t51\t52\t100\t35\t0\n')
        f.write('chr1\t222\t223\t0\t0\t100\n')
    
    # Cell 1 - avg. meth in gene 1 = 25%
    with gzip.open(methpath / "cell1.GCHN.bismark.cov", "wt") as f:
        f.write('chr1\t31\t32\t25\t25\t75\n')
        f.write('chr1\t41\t42\t25\t1\t3\n')

    # Cell 2 - avg. meth in gene 2 = 44%
    with gzip.open(methpath / "cell2.GCHN.bismark.cov", "wt") as f:
        f.write('chr1\t231\t232\t44\t44\t56\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.GCHN.bismark.cov", "wt") as f:
        f.write('chr1\t31\t32\t100\t100\t0\n')
        f.write('chr1\t41\t42\t100\t4\t0\n')
        f.write('chr1\t231\t232\t0\t0\t100\n')
    return methpath

@pytest.fixture
def bismarkrep_path(tmp_path):
    '''
    simulate bismark CpG report files
    tsv with chrom, pos, strand, methylcount, non methylcount, c-context, trinucleotide context
    these are 0-based.
    '''
    methpath = tmp_path / "bismarkrep"
    methpath.mkdir()

    # Cell 1 - avg. meth in gene 1 = 50%
    with gzip.open(methpath / "cell1.WCGN.bismark.rep.gz", "wt") as f:
        f.write('chr1\t20\t+\t25\t75\tCG\tCGN\n')
        f.write('chr1\t50\t+\t75\t25\tCG\tCGN\n')

    # Cell 2 - avg. meth in gene 2 = 62.5%
    with gzip.open(methpath / "cell2.WCGN.bismark.rep.gz", "wt") as f:
        f.write('chr1\t221\t+\t2\t2\tCG\tCGN\n')
        f.write('chr1\t244\t+\t12\t4\tCG\tCGN\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.WCGN.bismark.rep.gz", "wt") as f:
        f.write('chr1\t20\t+\t21\t0\tCG\tCGN\n')
        f.write('chr1\t50\t+\t35\t0\tCG\tCGN\n')
        f.write('chr1\t221\t+\t0\t100\tCG\tCGN\n')
    
    # Cell 1 - avg. meth in gene 1 = 25%
    with gzip.open(methpath / "cell1.GCHN.bismark.rep.gz", "wt") as f:
        f.write('chr1\t30\t+\t25\t75\tCG\tCGN\n')
        f.write('chr1\t40\t+\t1\t3\tCG\tCGN\n')

    # Cell 2 - avg. meth in gene 2 = 44%
    with gzip.open(methpath / "cell2.GCHN.bismark.rep.gz", "wt") as f:
        f.write('chr1\t230\t+\t44\t56\tCG\tCGN\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.GCHN.bismark.rep.gz", "wt") as f:
        f.write('chr1\t30\t+\t100\t0\tCG\tCGN\n')
        f.write('chr1\t40\t+\t4\t0\tCG\tCGN\n')
        f.write('chr1\t230\t+\t0\t100\tCG\tCGN\n')
    return methpath


@pytest.fixture
def bedmethyl_path(tmp_path):
    '''
    simulate bedmethyl files.
    tsv with chrom, start, end, name, score strand,
    thickstart, thickend, color, validcov, %mod, count_mod
    count_can, count_other, count_delete, count_fail, count_diff, count_nocall
    columns.
    these are 0-based.
    '''
    methpath = tmp_path / "bedmethyl"
    methpath.mkdir()

    # Cell 1 - avg. meth in gene 1 = 50%
    with gzip.open(methpath / "cell1.WCGN.bed.gz", "wt") as f:
        f.write('chr1\t20\t21\tm\t100\t+\t20\t21\t255,0,0\t100\t25.00\t25\t75\t0\t0\t0\t0\t0\n')
        f.write('chr1\t50\t51\tm\t100\t+\t50\t51\t255,0,0\t100\t75.00\t75\t25\t0\t0\t0\t0\t0\n')

    # Cell 2 - avg. meth in gene 2 = 62.5%
    with gzip.open(methpath / "cell2.WCGN.bed.gz", "wt") as f:
        f.write('chr1\t221\t222\tm\t4\t+\t221\t222\t255,0,0\t4\t50.00\t2\t2\t0\t0\t0\t0\t0\n')
        f.write('chr1\t244\t245\tm\t16\t+\t224\t225\t255,0,0\t16\t75.00\t12\t4\t0\t0\t0\t0\t0\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.WCGN.bed.gz", "wt") as f:
        f.write('chr1\t20\t21\tm\t21\t+\t20\t21\t255,0,0\t21\t100.00\t21\t0\t0\t0\t0\t0\t0\n')
        f.write('chr1\t50\t51\tm\t35\t+\t50\t51\t255,0,0\t35\t100.00\t35\t0\t0\t0\t0\t0\t0\n')
        f.write('chr1\t221\t222\tm\t100\t+\t221\t222\t255,0,0\t100\t0.00\t0\t100\t0\t0\t0\t0\t0\n')
    
    # Cell 1 - avg. meth in gene 1 = 25%
    with gzip.open(methpath / "cell1.GCHN.bed.gz", "wt") as f:
        f.write('chr1\t30\t31\tm\t100\t+\t30\t31\t255,0,0\t100\t25.00\t25\t75\t0\t0\t0\t0\t0\n')
        f.write('chr1\t40\t41\tm\t4\t+\t40\t41\t255,0,0\t4\t25.00\t1\t3\t0\t0\t0\t0\t0\n')

    # Cell 2 - avg. meth in gene 2 = 44%
    with gzip.open(methpath / "cell2.GCHN.bed.gz", "wt") as f:
        f.write('chr1\t230\t231\tm\t100\t+\t230\t231\t255,0,0\t100\t44.00\t44\t56\t0\t0\t0\t0\t0\n')
    
    # Cell 3 - avg. meth in gene 1 = 100%, avg. meth in gene 2 = 0%
    with gzip.open(methpath / "cell3.GCHN.bed.gz", "wt") as f:
        f.write('chr1\t30\t31\tm\t100\t+\t30\t31\t255,0,0\t100\t100.00\t100\t0\t0\t0\t0\t0\t0\n')
        f.write('chr1\t40\t41\tm\t4\t+\t40\t41\t255,0,0\t4\t100.00\t4\t0\t0\t0\t0\t0\t0\n')
        f.write('chr1\t230\t231\tm\t100\t+\t230\t231\t255,0,0\t100\t0.00\t0\t100\t0\t0\t0\t0\t0\n')
    return methpath


@pytest.fixture
def rna_path(tmp_path):
    '''
    simulate transcriptome.tsv files
    tsv with columns Geneid, Chr, Start, End, Strand, Length
    '''
    rnapath = tmp_path / "transcriptome"
    rnapath.mkdir()

    with open(rnapath / 'count1.tsv', 'w') as f:
        f.write("#Program featureCounts comment line.\n")
        f.write('Geneid\tChr\tStart\tEnd\tStrand\tLength\tcell1_rna\tcell2_rna\n')
        f.write('gene1\tchr1\t0\t100\t+\t101\t50\t100\n')
        f.write('gene2\tchr1\t200\t250\t+\t51\t20\t100\n')
    
    with open(rnapath / 'count2.tsv', 'w') as f:
        f.write("#Program featureCounts comment line.\n")
        f.write('Geneid\tChr\tStart\tEnd\tStrand\tLength\tcell3_rna\n')
        f.write('gene1\tchr1\t0\t100\t+\t101\t33\n')
        f.write('gene2\tchr1\t200\t250\t+\t51\t33\n')
    
    return rnapath

@pytest.fixture
def bed_path(tmp_path):
    '''
    simulate regions.bed files
    '''
    bedpath = tmp_path / "regions"
    bedpath.mkdir()

    with open(bedpath / "gene1.bed", "w") as f:
        f.write('chr1\t0\t100\tgene1\t0\t+\n')

    with gzip.open(bedpath / "gene1.bed.gz", "wt") as f:
        f.write('chr1\t0\t100\tgene1\t0\t+\n')
    
    # gene2: 200-250 bps
    with open(bedpath / "gene2.bed", "w") as f:
        f.write('chr1\t200\t250\tgene2\t0\t+\n')

    with gzip.open(bedpath / "gene2.bed.gz", "wt") as f:
        f.write('chr1\t200\t250\tgene2\t0\t+\n')

    # blacklist regions
    with open(bedpath / "blacklist1.bed", 'w') as f:
        f.write('chr1\t40\t60\n')

    with open(bedpath / "blacklist2.bed", 'w') as f:
        f.write('chr1\t235\t250\n')

    with open(bedpath / "genome.chromsizes", 'w') as f:
        f.write('chr1\t500\n')

    return bedpath

def mu_to_dense(mu, var_name):
    """
    Takes a mu data object, and returns a dense matrix.
    Instead of directly calling . mu[var_name].X.todense(),
    this function retains distiction between NaN and 0 values.
    """
    if var_name not in mu.mod.keys():
        raise ValueError(f"'{var_name}' not found in muData object.")
    X = mu.mod[var_name].X
    X = X.tocoo()
    dense = np.full(X.shape, np.nan, dtype=X.dtype)
    dense[X.row, X.col] = X.data
    return dense

@pytest.fixture
def nome_path(tmp_path):
    '''
    Create empty chromsizes, empty GCHN and WCGN files for testing.
    '''

    chromsizes = tmp_path / "chromsizes.txt"
    chromsizes.touch()

    gchn = tmp_path / "GCHN"
    gchn.mkdir()
    wcgn = tmp_path / "WCGN"
    wcgn.mkdir()
    for ix in range(1, 4):
        (gchn / f"cell{ix}.GCHN.tsv.gz").touch()
        (wcgn / f"cell{ix}.WCGN.tsv.gz").touch()

    return tmp_path