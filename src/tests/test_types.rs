use crate::types::{MethFileType, MethRegion};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_allcools_parse_line() {
        let region = MethFileType::AllCools
            .parse_line("chr1\t1\t+\tCG\t1\t1\t1")
            .unwrap()
            .unwrap();
        assert_eq!(
            region,
            MethRegion {
                chrom: "chr1".into(),
                pos: 0,
                meth: 1,
                total: 1
            }
        );
    }

    #[test]
    fn test_allcools_too_few_fields() {
        assert!(
            MethFileType::AllCools
                .parse_line("chr1\t1\t+\tCG\t1\t1")
                .is_err()
        );
    }

    #[test]
    fn test_allcools_too_many_fields() {
        assert!(
            MethFileType::AllCools
                .parse_line("chr1\t1\t+\tCG\t1\t1\t1\textra")
                .is_err()
        );
    }

    #[test]
    fn test_allcools_invalid_position() {
        assert!(
            MethFileType::AllCools
                .parse_line("chr1\tNA\t+\tCG\t1\t1\t1")
                .is_err()
        );
    }

    #[test]
    fn test_allcools_invalid_meth() {
        assert!(
            MethFileType::AllCools
                .parse_line("chr1\t1\t+\tCG\tNA\t1\t1")
                .is_err()
        );
    }

    #[test]
    fn test_allcools_invalid_total() {
        assert!(
            MethFileType::AllCools
                .parse_line("chr1\t1\t+\tCG\t1\tNA\t1")
                .is_err()
        );
    }

    #[test]
    fn test_methyldackel_parse_line() {
        let region = MethFileType::MethylDackel
            .parse_line("chr1\t0\t1\t100\t1\t1")
            .unwrap()
            .unwrap();
        // MethylDackel positions are already 0-based.
        assert_eq!(
            region,
            MethRegion {
                chrom: "chr1".into(),
                pos: 0,
                meth: 1,
                total: 1
            }
        );
    }

    #[test]
    fn test_methyldackel_too_few_fields() {
        assert!(
            MethFileType::MethylDackel
                .parse_line("chr1\t0\t1\t100\t1")
                .is_err()
        );
    }

    #[test]
    fn test_methyldackel_too_many_fields() {
        assert!(
            MethFileType::MethylDackel
                .parse_line("chr1\t0\t1\t100\t1\t1\textra")
                .is_err()
        );
    }

    #[test]
    fn test_methyldackel_invalid_position() {
        assert!(
            MethFileType::MethylDackel
                .parse_line("chr1\tNA\t1\t100\t1\t1")
                .is_err()
        );
    }

    #[test]
    fn test_methyldackel_invalid_meth() {
        assert!(
            MethFileType::MethylDackel
                .parse_line("chr1\t0\t1\t100\tNA\t1")
                .is_err()
        );
    }

    #[test]
    fn test_methyldackel_invalid_total() {
        assert!(
            MethFileType::MethylDackel
                .parse_line("chr1\t0\t1\t100\t1\tNA")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcov_parse_line() {
        let region = MethFileType::BismarkCov
            .parse_line("chr1\t1\t2\t100\t1\t0")
            .unwrap()
            .unwrap();
        // BismarkCov is 1-based; parse_line converts to 0-based. total = meth + unmeth.
        assert_eq!(
            region,
            MethRegion {
                chrom: "chr1".into(),
                pos: 0,
                meth: 1,
                total: 1
            }
        );
    }

    #[test]
    fn test_bismarkcov_too_few_fields() {
        assert!(
            MethFileType::BismarkCov
                .parse_line("chr1\t1\t2\t100\t1")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcov_too_many_fields() {
        assert!(
            MethFileType::BismarkCov
                .parse_line("chr1\t1\t2\t100\t1\t0\textra")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcov_invalid_position() {
        assert!(
            MethFileType::BismarkCov
                .parse_line("chr1\tNA\t2\t100\t1\t0")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcov_invalid_meth() {
        assert!(
            MethFileType::BismarkCov
                .parse_line("chr1\t1\t2\t100\tNA\t0")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcov_invalid_unmeth() {
        assert!(
            MethFileType::BismarkCov
                .parse_line("chr1\t1\t2\t100\t1\tNA")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcpgreport_parse_line() {
        let region = MethFileType::BismarkCpGReport
            .parse_line("chr1\t0\t+\t1\t0\tCG\tCGH")
            .unwrap()
            .unwrap();
        assert_eq!(
            region,
            MethRegion {
                chrom: "chr1".into(),
                pos: 0,
                meth: 1,
                total: 1
            }
        );
    }

    #[test]
    fn test_bismarkcpgreport_skips_zero_coverage() {
        // Positions with no methylated/unmethylated calls at all are skipped.
        let region = MethFileType::BismarkCpGReport
            .parse_line("chr1\t0\t+\t0\t0\tCG\tCGH")
            .unwrap();
        assert_eq!(region, None);
    }

    #[test]
    fn test_bismarkcpgreport_too_few_fields() {
        assert!(
            MethFileType::BismarkCpGReport
                .parse_line("chr1\t0\t+\t1\t0\tCG")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcpgreport_too_many_fields() {
        assert!(
            MethFileType::BismarkCpGReport
                .parse_line("chr1\t0\t+\t1\t0\tCG\tCGH\textra")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcpgreport_invalid_position() {
        assert!(
            MethFileType::BismarkCpGReport
                .parse_line("chr1\tNA\t+\t1\t0\tCG\tCGH")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcpgreport_invalid_meth() {
        assert!(
            MethFileType::BismarkCpGReport
                .parse_line("chr1\t0\t+\tNA\t0\tCG\tCGH")
                .is_err()
        );
    }

    #[test]
    fn test_bismarkcpgreport_invalid_unmeth() {
        assert!(
            MethFileType::BismarkCpGReport
                .parse_line("chr1\t0\t+\t1\tNA\tCG\tCGH")
                .is_err()
        );
    }

    #[test]
    fn test_bedmethyl_parse_line() {
        let line = "chr1\t0\t1\tm\t1\t+\t0\t1\t255,0,0\t1\t100.00\t1\t0\t0\t0\t0\t0\t0";
        let region = MethFileType::BedMethyl.parse_line(line).unwrap().unwrap();
        assert_eq!(
            region,
            MethRegion {
                chrom: "chr1".into(),
                pos: 0,
                meth: 1,
                total: 1
            }
        );
    }

    #[test]
    fn test_bedmethyl_too_few_fields() {
        assert!(
            MethFileType::BedMethyl
                .parse_line("chr1\t0\t1\tm\t1\t+")
                .is_err()
        );
    }

    #[test]
    fn test_bedmethyl_too_many_fields() {
        let line = "chr1\t0\t1\tm\t1\t+\t0\t1\t255,0,0\t1\t100.00\t1\t0\t0\t0\t0\t0\t0\textra";
        assert!(MethFileType::BedMethyl.parse_line(line).is_err());
    }

    #[test]
    fn test_bedmethyl_invalid_position() {
        let line = "chr1\tNA\t1\tm\t1\t+\t0\t1\t255,0,0\t1\t100.00\t1\t0\t0\t0\t0\t0\t0";
        assert!(MethFileType::BedMethyl.parse_line(line).is_err());
    }

    #[test]
    fn test_bedmethyl_invalid_meth() {
        // column 11 (modification count) is invalid, column 9 (coverage) stays valid.
        let line = "chr1\t0\t1\tm\t1\t+\t0\t1\t255,0,0\t1\t100.00\tNA\t0\t0\t0\t0\t0\t0";
        assert!(MethFileType::BedMethyl.parse_line(line).is_err());
    }

    #[test]
    fn test_bedmethyl_invalid_total() {
        // column 9 (coverage) is invalid, column 11 (modification count) stays valid.
        let line = "chr1\t0\t1\tm\t1\t+\t0\t1\t255,0,0\tNA\t100.00\t1\t0\t0\t0\t0\t0\t0";
        assert!(MethFileType::BedMethyl.parse_line(line).is_err());
    }

    #[test]
    fn test_parse_line_skips_comments_regardless_of_type() {
        for methtype in [
            MethFileType::AllCools,
            MethFileType::MethylDackel,
            MethFileType::BismarkCov,
            MethFileType::BismarkCpGReport,
            MethFileType::BedMethyl,
        ] {
            assert_eq!(methtype.parse_line("# a comment line").unwrap(), None);
        }
    }
}
