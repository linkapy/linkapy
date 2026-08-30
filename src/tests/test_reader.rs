use crate::reader::{decide_methtype, is_gzipped, parse_chromsizes, parse_region, read_meth};
use crate::types::{MethFileType, MethRegion};
use std::io::BufRead;
use std::path::Path;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_region() {
        // get path to test data relative to test code.
        // This ensure cargo test can be ran from anywhere.
        let test_path = Path::new(file!());
        let bedfile = test_path.parent().unwrap().join("data/region.bed");
        let bedgzfile = test_path.parent().unwrap().join("data/region.bed.gz");

        // Parse bed file
        let bedregions = parse_region(bedfile.to_string_lossy().into_owned(), "bed".to_string());
        let bedgzregions = parse_region(
            bedgzfile.to_string_lossy().into_owned(),
            "bedgz".to_string(),
        );
        // bed file.
        assert_eq!(bedregions.len(), 1);
        assert_eq!(bedregions[0].chrom, "chr1");
        assert_eq!(bedregions[0].start, 100);
        assert_eq!(bedregions[0].end, 200);
        assert_eq!(bedregions[0].name, "chr1:100-200");
        assert_eq!(bedregions[0].class, "bed");
        // bed gz file.
        assert_eq!(bedgzregions.len(), 1);
        assert_eq!(bedgzregions[0].chrom, "chr1");
        assert_eq!(bedgzregions[0].start, 100);
        assert_eq!(bedgzregions[0].end, 200);
        assert_eq!(bedgzregions[0].name, "chr1:100-200");
        assert_eq!(bedgzregions[0].class, "bedgz");
    }

    #[test]
    fn test_is_gzipped() {
        let test_path = Path::new(file!());
        let bedfile = test_path.parent().unwrap().join("data/region.bed");
        let bedgzfile = test_path.parent().unwrap().join("data/region.bed.gz");
        assert_eq!(is_gzipped(&bedfile.to_string_lossy()).unwrap(), false);
        assert_eq!(is_gzipped(&bedgzfile.to_string_lossy()).unwrap(), true);
    }

    #[test]
    fn test_decide_methtype() {
        let test_path = Path::new(file!());
        let allcoolsf = test_path.parent().unwrap().join("data/methf_allcools");
        let bismarkcov = test_path.parent().unwrap().join("data/methf_bismarkcov");
        let cpgrep = test_path.parent().unwrap().join("data/methf_cpgrep");
        let methyldackel = test_path.parent().unwrap().join("data/methf_methyldackel");
        let bedmethyl = test_path.parent().unwrap().join("data/methf_bedmethyl");
        let pairs: Vec<(&str, MethFileType)> = vec![
            (allcoolsf.to_str().unwrap(), MethFileType::AllCools),
            (bismarkcov.to_str().unwrap(), MethFileType::BismarkCov),
            (cpgrep.to_str().unwrap(), MethFileType::BismarkCpGReport),
            (methyldackel.to_str().unwrap(), MethFileType::MethylDackel),
            (bedmethyl.to_str().unwrap(), MethFileType::BedMethyl),
        ];
        for (f, expected) in pairs {
            let reader = std::io::BufReader::new(std::fs::File::open(f).unwrap());
            let firstline = reader
                .lines()
                .filter_map(|l| l.ok())
                .find(|line| !line.trim_start().starts_with('#'));
            let methtype = decide_methtype(firstline);
            assert_eq!(methtype, expected);
        }
    }

    #[test]
    fn test_read_meth() {
        let test_path = Path::new(file!());
        let allcoolsf = test_path.parent().unwrap().join("data/methf_allcools");
        let bismarkcov = test_path.parent().unwrap().join("data/methf_bismarkcov");
        let cpgrep = test_path.parent().unwrap().join("data/methf_cpgrep");
        let methyldackel = test_path.parent().unwrap().join("data/methf_methyldackel");
        let bedmethyl = test_path.parent().unwrap().join("data/methf_bedmethyl");
        let exp_mr = vec![
            MethRegion {
                chrom: "chr1".to_string(),
                pos: 0,
                meth: 1,
                total: 1,
            },
            MethRegion {
                chrom: "chr1".to_string(),
                pos: 2,
                meth: 0,
                total: 1,
            },
        ];

        for f in vec![
            allcoolsf.to_str().unwrap(),
            bismarkcov.to_str().unwrap(),
            cpgrep.to_str().unwrap(),
            methyldackel.to_str().unwrap(),
            bedmethyl.to_str().unwrap(),
        ] {
            let methregions = read_meth(f);
            assert_eq!(methregions, exp_mr);
        }
    }

    #[test]
    fn test_parse_chromsizes() {
        let test_path = Path::new(file!());
        let chromsizef = test_path.parent().unwrap().join("data/chromsizes.txt");

        let regions = parse_chromsizes(&chromsizef.to_string_lossy(), 1000);
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[1].chrom, "chr1");
        assert_eq!(regions[1].start, 1000);
        assert_eq!(regions[1].end, 2000);
        assert_eq!(regions[1].name, "chr1:1000-2000");
        assert_eq!(regions[1].class, "bin");
    }

    #[test]
    #[should_panic(expected = "empty or only contains comments")]
    fn test_decide_methtype_none_panics() {
        decide_methtype(None);
    }

    #[test]
    #[should_panic(expected = "Could not decide between MethylDackel")]
    fn test_decide_methtype_ambiguous_six_columns_panics() {
        // meth=1, prop_cov=1 => MethylDackel-style percent would be 100, BismarkCov-style 50.
        // 10.0 matches neither, so the filetype can't be decided.
        decide_methtype(Some("chr1\t1\t2\t10.0\t1\t1".to_string()));
    }

    #[test]
    #[should_panic(expected = "Could not decide methylation filetype, as it has 5 columns")]
    fn test_decide_methtype_unsupported_column_count_panics() {
        decide_methtype(Some("chr1\t1\t2\t3\t4".to_string()));
    }

    #[test]
    #[should_panic(expected = "Could not parse BedMethyl coverage and modification counts")]
    fn test_decide_methtype_bedmethyl_unparseable_counts_panics() {
        let line = "chr1\t0\t1\tm\t1\t+\t0\t1\t255,0,0\tNA\t100.00\tNA\t0\t0\t0\t0\t0\t0";
        decide_methtype(Some(line.to_string()));
    }

    #[test]
    fn test_read_meth_gzipped() {
        let test_path = Path::new(file!());
        let plain = test_path.parent().unwrap().join("data/methf_methyldackel");
        let gzipped = test_path
            .parent()
            .unwrap()
            .join("data/methf_methyldackel.gz");
        assert_eq!(
            read_meth(&plain.to_string_lossy()),
            read_meth(&gzipped.to_string_lossy())
        );
    }

    #[test]
    fn test_read_meth_skips_track_and_comment_lines() {
        let test_path = Path::new(file!());
        let f = test_path
            .parent()
            .unwrap()
            .join("data/methf_methyldackel_trackskip");
        let methregions = read_meth(&f.to_string_lossy());
        assert_eq!(
            methregions,
            vec![
                MethRegion {
                    chrom: "chr1".to_string(),
                    pos: 0,
                    meth: 1,
                    total: 1
                },
                MethRegion {
                    chrom: "chr1".to_string(),
                    pos: 2,
                    meth: 0,
                    total: 1
                },
            ]
        );
    }

    #[test]
    #[should_panic(expected = "Error reading file")]
    fn test_read_meth_missing_file_panics() {
        read_meth("data/does_not_exist.tsv");
    }

    #[test]
    #[should_panic(expected = "Error reading file")]
    fn test_parse_region_missing_file_panics() {
        parse_region("data/does_not_exist.bed".to_string(), "bed".to_string());
    }

    #[test]
    fn test_parse_region_with_name_field() {
        let test_path = Path::new(file!());
        let f = test_path.parent().unwrap().join("data/region_named.bed");
        let regions = parse_region(f.to_string_lossy().into_owned(), "bed".to_string());
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].name, "regionA");
    }

    #[test]
    #[should_panic(expected = "Error reading file")]
    fn test_parse_chromsizes_invalid_utf8_panics() {
        let test_path = Path::new(file!());
        let f = test_path
            .parent()
            .unwrap()
            .join("data/chromsizes_badutf8.txt");
        parse_chromsizes(&f.to_string_lossy(), 1000);
    }

    #[test]
    #[should_panic(expected = "Error reading file")]
    fn test_parse_region_invalid_utf8_panics() {
        let test_path = Path::new(file!());
        let f = test_path.parent().unwrap().join("data/region_badutf8.bed");
        parse_region(f.to_string_lossy().into_owned(), "bed".to_string());
    }
}
