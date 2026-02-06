use pyo3::prelude::*;
use pyo3::types::{PyList, PyModule};
use std::fs::File;
use std::io::{Write};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::cmp::Ordering;
use std::collections::HashMap;
use sprs::{CsMat, TriMat};
use sprs::io::write_matrix_market;
use crate::types::{Region};
use crate::reader::{read_meth, parse_chromsizes, parse_region};

#[allow(clippy::too_many_arguments)]
#[pyfunction]
pub fn parse_cools(
    py: Python<'_>,
    _coolfiles: Py<PyList>,
    _regions: Py<PyList>,
    _blacklist: Py<PyList>,
    _regionlabels: Py<PyList>,
    threads: usize,
    prefix: &str,
    chromsizes: &str,
    binsize: u32
) -> PyResult<()> {
    // Set up the logging from python
    let logging = PyModule::import(py, "logging")?;
    let logger = logging.call_method1("getLogger", ())?;

    let mut coolfiles: Vec<String> = Vec::new();
    let mut regions: Vec<String> = Vec::new();
    let mut regionlabels: Vec<String> = Vec::new();
    let mut blacklist: Vec<String> = Vec::new();

    Python::with_gil(|py| {
        coolfiles = _coolfiles.extract(py).expect("Failed to retrieve allcoolfiles.");
        regions = _regions.extract(py).expect("Failed to retrieve regions.");
        regionlabels = _regionlabels.extract(py).expect("Failed to retrieve region labels.");
        blacklist = _blacklist.extract(py).expect("Failed to retrieve blacklist regions.");
    });
    // regions and regionlabels should always be same length.
    assert_eq!(regions.len(), regionlabels.len());


    let blacklist_regions: Option<Vec<Region>> = if blacklist.is_empty() {
        logger.call_method1("info", ("\'keep_cool\': No blacklist provided.",))?;
        None
    } else {
        let mut blacklist_regions: Vec<Region> = Vec::new();
        for _b in blacklist.into_iter() {
            blacklist_regions.extend(parse_region(_b, "blacklist".to_string()));
        }
        logger.call_method1(
            "info",
            (format!("\'keep_cool\': Blacklist(s) parsed. {} regions.", blacklist_regions.len()),)
        )?;
        Some(blacklist_regions)
    };

    let parsed_regions = if regions.is_empty() {
            logger.call_method1("info", ("\'keep_cool\': running in chromsize mode.",))?;
            parse_chromsizes(chromsizes, binsize)
        } else {
            logger.call_method1("info", ("\'keep_cool\': running in regions mode.",))?;
            // Parse regions.
            let mut parsed_regions: Vec<Region> = Vec::new();
            for (_r, _l) in regions.into_iter().zip(regionlabels.into_iter()) {
                parsed_regions.extend(parse_region(_r, _l));
            }
            // Sort per chromosome and start position.
            parsed_regions.sort_by(|a, b| {
                // First, compare by `chrom`
                let chrom_order = a.chrom.cmp(&b.chrom);
                if chrom_order != Ordering::Equal {
                    return chrom_order;
                }
                a.start[0].cmp(&b.start[0])
            });
            parsed_regions
        };

    logger.call_method1(
        "info",
        (format!("\'keep_cool\': Found {} regions.", parsed_regions.len()),)
    )?;
    let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    logger.call_method1(
        "info",
        (format!("\'keep_cool\': Starting pool with {} threads.", threads),)
    )?;
    let blacklist_by_chrom: HashMap<String, Vec<(u32, u32)>> = match &blacklist_regions {
        Some(blacklist) => {
            let mut map: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
            for bl in blacklist.iter() {
                let start = bl.start[0];
                let end = *bl.end.last().unwrap();
                map.entry(bl.chrom.clone()).or_default().push((start, end));
            }
            for intervals in map.values_mut() {
                intervals.sort_by_key(|(s, _)| *s);
            }
            map
        }
        None => HashMap::new(),
    };
    // Metrics
    let aggregated_metrics: Vec<Vec<((f32, f32, f32), f32)>>  = pool.install(|| {
        coolfiles
            .par_iter()
            .map(|methfile| {
                let mut methregions = read_meth(methfile);
                
                // Sort per chromosome and position (cannot really assume this at this point, though most times it'll be.)
                methregions.sort_by(|a, b| {
                    let chrom_order = a.chrom.cmp(&b.chrom);
                    if chrom_order != Ordering::Equal {
                        return chrom_order;
                    }
                    a.pos.cmp(&b.pos)
                });

                // Index regions by chromosome
                let mut by_chrom: HashMap<String, (usize, usize)> = HashMap::new();
                let mut start = 0;
                for i in 1..=methregions.len() {
                    if i == methregions.len() || methregions[i].chrom != methregions[start].chrom {
                        by_chrom.insert(methregions[start].chrom.clone(), (start, i));
                        start = i;
                    }
                }

                parsed_regions
                    .par_iter()
                    .map(|region| {
                        let (meth_sum, total_sum, sites, frac_sum, frac_count) = if let Some((s, e)) = by_chrom.get(&region.chrom) {
                            let chrom_regions = &methregions[*s..*e];
                            let start_idx = chrom_regions.partition_point(|x| x.pos < region.start[0]);
                            let end_idx = chrom_regions.partition_point(|x| x.pos < *region.end.last().unwrap());
                            let mut meth_sum = f32::NAN;
                            let mut total_sum = f32::NAN;
                            let mut sites = f32::NAN;
                            let mut frac_sum = 0.0f32;
                            let mut frac_count = 0.0f32;

                            let blacklist_intervals = blacklist_by_chrom.get(&region.chrom);
                            for x in &chrom_regions[start_idx..end_idx] {
                                let is_blacklisted = if let Some(intervals) = blacklist_intervals {
                                    let pos = x.pos;
                                    let idx = intervals.partition_point(|(s, _)| *s <= pos);
                                    if idx == 0 {
                                        false
                                    } else {
                                        let (s, e) = intervals[idx - 1];
                                        pos >= s && pos < e
                                    }
                                } else {
                                    false
                                };
                                if is_blacklisted {
                                    continue;
                                }
                                let meth = x.meth as f32;
                                let total = x.total as f32;
                                let frac = meth / total; // total will never be zero.
                                meth_sum = if meth_sum.is_nan() { meth } else { meth_sum + meth };
                                total_sum = if total_sum.is_nan() { total } else { total_sum + total };
                                sites = if sites.is_nan() { 1.0 } else { sites + 1.0 };
                                if !frac.is_nan() {
                                    frac_sum += frac;
                                    frac_count += 1.0;
                                }
                            }
                            (meth_sum, total_sum, sites, frac_sum, frac_count)
                        } else {
                            (f32::NAN, f32::NAN, f32::NAN, 0.0, 0.0)
                        };

                        let mean_fraction = if frac_count > 0.0 {
                            frac_sum / frac_count
                        } else {
                            f32::NAN
                        };

                        ((meth_sum, total_sum, sites), mean_fraction)
                    })
            .collect()
            })
    .collect()
    });
    
    let regvals: Vec<Vec<(f32, f32, f32)>> = aggregated_metrics.iter().map(|v| v.iter().map(|(vals, _)| *vals).collect()).collect();
    let fractions_vec: Vec<Vec<f32>> = aggregated_metrics.iter().map(|v| v.iter().map(|(_, fracs)| fracs.clone()).collect()).collect();
    let fracm = frac_to_sparse(fractions_vec);
    let (methm, covm, sitem) = tupvec_to_sparse(regvals);
    logger.call_method1(
        "info",
        (format!("\'keep_cool\': Finished parsing {} files.", coolfiles.len()),)
    )?;
    // Define output files taken the prefix.
    let ometh = format!("{}.meth.mtx", prefix);
    let ocov = format!("{}.cov.mtx", prefix);
    let osite = format!("{}.site.mtx", prefix);
    let ofrac = format!("{}.frac.mtx", prefix);

    let oregionfile: String = format!("{}.regions.tsv", prefix);
    let ocellfile: String = format!("{}.cells.tsv", prefix);
    write_matrix_market(ometh, &methm).unwrap();
    write_matrix_market(ocov, &covm).unwrap();
    write_matrix_market(osite, &sitem).unwrap();
    write_matrix_market(ofrac, &fracm).unwrap();

    logger.call_method1(
        "info",
        (format!("\'keep_cool\': Finished writing matrices with prefix {}.", prefix),)
    )?;

    let mut ofile = File::create(oregionfile).unwrap();
    writeln!(ofile, "chrom\tstart\tend\tname\tclass").unwrap();
    for region in parsed_regions {
        writeln!(ofile, "{}\t{}\t{}\t{}\t{}", region.chrom, region.start[0], *region.end.last().unwrap(), region.name, region.class).unwrap();
    }
    let mut ofile = File::create(ocellfile).unwrap();
    for coolfile in coolfiles {
        writeln!(ofile, "{}", coolfile).unwrap();
    }
    logger.call_method1(
        "info",
        (format!("\'keep_cool\': Finished writing metadata with prefix {}.", prefix),)
    )?;

    Ok(())
}

pub fn frac_to_sparse(dense: Vec<Vec<f32>>) -> CsMat<f32> {
    let max_row = dense.len();
    let max_col = dense.iter().map(|row| row.len()).max().unwrap_or(0);

    let mut mat = TriMat::new((max_row, max_col));

    for (i, row) in dense.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            if !v.is_nan() {
                mat.add_triplet(i, j, v);
            }
        }
    }
    mat.to_csr()
}

pub fn tupvec_to_sparse(dense: Vec<Vec<(f32, f32, f32)>>) -> (CsMat<f32>, CsMat<f32>, CsMat<f32>) {
    let max_row = dense.len();
    let max_col = dense.iter().map(|row| row.len()).max().unwrap_or(0);

    let mut mat1 = TriMat::new((max_row, max_col));
    let mut mat2 = TriMat::new((max_row, max_col));
    let mut mat3 = TriMat::new((max_row, max_col));

    for (i, row) in dense.iter().enumerate() {
        for (j, &(v1, v2, v3)) in row.iter().enumerate() {
            if !v1.is_nan() {
                mat1.add_triplet(i, j, v1);
            }
            if !v2.is_nan() {
                mat2.add_triplet(i, j, v2);
            }
            if !v3.is_nan() {
                mat3.add_triplet(i, j, v3);
            }
        }
    }
    (mat1.to_csr(), mat2.to_csr(), mat3.to_csr())
}