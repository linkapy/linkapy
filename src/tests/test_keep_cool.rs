use crate::keep_cool::{compute_cell_qc, finalize_region_qc, frac_to_sparse, tupvec_to_sparse};
use crate::types::MethRegion;

#[cfg(test)]
mod tests {
    use super::*;

    fn mr(meth: u32, total: u32) -> MethRegion {
        MethRegion {
            chrom: "chr1".into(),
            pos: 0,
            meth,
            total,
        }
    }

    #[test]
    fn test_frac_to_sparse() {
        let dense = vec![
            vec![1.0, 2.0, f32::NAN],
            vec![0.0, 2.0, 3.0],
            vec![f32::NAN, f32::NAN, 4.0],
        ];
        let sparse = frac_to_sparse(dense);
        assert_eq!(sparse.shape(), (3, 3));
        assert_eq!(sparse.nnz(), 6);
        assert_eq!(sparse.get(0, 0), Some(&1.0));
        assert_eq!(sparse.get(0, 1), Some(&2.0));
        assert_eq!(sparse.get(1, 0), Some(&0.0));
        assert_eq!(sparse.get(1, 1), Some(&2.0));
        assert_eq!(sparse.get(1, 2), Some(&3.0));
        assert_eq!(sparse.get(2, 2), Some(&4.0));
    }

    #[test]
    fn test_tupvec_to_sparse() {
        let dense = vec![
            vec![(1.0, 1.0, 1.0), (2.0, 2.0, 2.0), (f32::NAN, f32::NAN, f32::NAN)],
            vec![(0.0, 0.0, 0.0), (2.0, 2.0, 2.0), (3.0, 3.0 ,3.0)],
            vec![(f32::NAN, f32::NAN, f32::NAN), (f32::NAN, f32::NAN, f32::NAN), (4.0, 4.0, 4.0)],
        ];
        let (mat1, mat2, mat3) = tupvec_to_sparse(dense);
        for sparse in &[mat1, mat2, mat3] {
            assert_eq!(sparse.shape(), (3,3));
            assert_eq!(sparse.nnz(), 6);
            assert_eq!(sparse.get(0, 0), Some(&1.0));
            assert_eq!(sparse.get(0, 1), Some(&2.0));
            assert_eq!(sparse.get(1, 0), Some(&0.0));
            assert_eq!(sparse.get(1, 1), Some(&2.0));
            assert_eq!(sparse.get(1, 2), Some(&3.0));
            assert_eq!(sparse.get(2, 2), Some(&4.0));
        }
    }

    #[test]
    fn test_compute_cell_qc_normal() {
        // meth/total: 5/10, 15/20 -> total_meth=20, total_depth=30
        let methregions = vec![mr(5, 10), mr(15, 20)];
        let qc = compute_cell_qc(&methregions, 3, 8, 10);
        assert_eq!(qc.n_sites, 2);
        assert_eq!(qc.n_chroms, 3);
        assert_eq!(qc.total_depth, 30);
        assert_eq!(qc.total_meth, 20);
        assert_eq!(qc.mean_coverage, 15.0); // 30 / 2
        assert!((qc.global_meth_frac - (20.0 / 30.0)).abs() < 1e-6);
        assert_eq!(qc.n_regions_covered, 8);
        assert_eq!(qc.region_coverage_rate, 0.8); // 8 / 10
    }

    #[test]
    fn test_compute_cell_qc_no_sites() {
        let methregions: Vec<MethRegion> = Vec::new();
        let qc = compute_cell_qc(&methregions, 0, 0, 10);
        assert_eq!(qc.n_sites, 0);
        assert_eq!(qc.total_depth, 0);
        assert_eq!(qc.total_meth, 0);
        assert!(qc.mean_coverage.is_nan());
        assert!(qc.global_meth_frac.is_nan());
        assert_eq!(qc.region_coverage_rate, 0.0);
    }

    #[test]
    fn test_compute_cell_qc_no_regions() {
        let methregions = vec![mr(1, 1)];
        let qc = compute_cell_qc(&methregions, 1, 0, 0);
        assert!(qc.region_coverage_rate.is_nan());
    }

    #[test]
    fn test_finalize_region_qc_normal() {
        // 2 regions, 3 files. Region 0 covered by all 3 cells (fracs 0.0, 0.5, 1.0),
        // region 1 covered by none.
        let n_cells_covered = vec![3, 0];
        let frac_sum = vec![0.0 + 0.5 + 1.0, 0.0];
        let frac_sumsq = vec![0.0 * 0.0 + 0.5 * 0.5 + 1.0 * 1.0, 0.0];
        let qc = finalize_region_qc(n_cells_covered, &frac_sum, &frac_sumsq, 3);

        assert_eq!(qc.n_cells_covered, vec![3, 0]);
        assert_eq!(qc.missingness[0], 0.0);
        assert_eq!(qc.missingness[1], 1.0);

        assert!((qc.mean_frac[0] - 0.5).abs() < 1e-6);
        // population variance of [0.0, 0.5, 1.0] around mean 0.5 = (0.25 + 0 + 0.25) / 3
        let expected_std = ((0.25f64 + 0.0 + 0.25) / 3.0).sqrt() as f32;
        assert!((qc.std_frac[0] - expected_std).abs() < 1e-6);

        assert!(qc.mean_frac[1].is_nan());
        assert!(qc.std_frac[1].is_nan());
    }

    #[test]
    fn test_finalize_region_qc_single_cell_has_zero_variance() {
        let n_cells_covered = vec![1];
        let frac_sum = vec![0.7];
        let frac_sumsq = vec![0.7 * 0.7];
        let qc = finalize_region_qc(n_cells_covered, &frac_sum, &frac_sumsq, 4);

        assert!((qc.mean_frac[0] - 0.7).abs() < 1e-6);
        assert_eq!(qc.std_frac[0], 0.0);
        assert_eq!(qc.missingness[0], 0.75); // 1 - 1/4
    }

    #[test]
    fn test_finalize_region_qc_no_files() {
        let n_cells_covered = vec![0];
        let frac_sum = vec![0.0];
        let frac_sumsq = vec![0.0];
        let qc = finalize_region_qc(n_cells_covered, &frac_sum, &frac_sumsq, 0);
        assert!(qc.missingness[0].is_nan());
    }
}