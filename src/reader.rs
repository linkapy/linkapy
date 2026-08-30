use crate::types::{MethFileType, MethRegion, Region};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

fn open_meth_reader(_f: &str) -> Box<dyn BufRead> {
    match is_gzipped(_f) {
        Ok(true) => Box::new(BufReader::new(MultiGzDecoder::new(File::open(_f).unwrap()))),
        Ok(false) => Box::new(BufReader::new(File::open(_f).unwrap())),
        Err(e) => panic!("Error reading file {}: {}", _f, e),
    }
}

pub fn read_meth(_f: &str) -> Vec<MethRegion> {
    let mut methregions: Vec<MethRegion> = Vec::with_capacity(1 << 16);

    let mut lines = open_meth_reader(_f).lines().map_while(|l| l.ok());
    // Decide methtype
    let mut methtype: Option<MethFileType> = None;
    for line in &mut lines {
        let trimmed = line.trim_start();
        if trimmed.starts_with('#') || trimmed.starts_with("track") {
            continue;
        }
        let mt = decide_methtype(Some(line.clone()));
        if let Some(region) = mt.parse_line(&line).unwrap() {
            methregions.push(region);
        }
        methtype = Some(mt);
        break;
    }
    let methtype = methtype.unwrap_or_else(|| decide_methtype(None));

    for line in lines {
        // Skip track lines
        if line.trim_start().starts_with("track") {
            continue;
        }

        match methtype.parse_line(&line).unwrap() {
            Some(region) => methregions.push(region),
            None => { /* e.g., skip empty or invalid CpGReport lines */ }
        }
    }
    methregions
}

pub fn parse_chromsizes(file: &str, binsize: u32) -> Vec<Region> {
    let mut regions: Vec<Region> = Vec::new();
    let reader = BufReader::new(File::open(file).unwrap());
    for line in reader.lines() {
        match line {
            Ok(line) => {
                let fields: Vec<&str> = line.split('\t').collect();
                let chrom = fields[0].to_string();
                let chromsize = fields[1].parse::<u32>().unwrap();
                let mut start = 0;
                let mut end = start + binsize;
                while end < chromsize {
                    regions.push(Region {
                        chrom: chrom.clone(),
                        start,
                        end,
                        name: format!("{}:{}-{}", chrom, start, end),
                        class: "bin".to_string(),
                    });
                    start = end;
                    end += binsize;
                    // Capture chromosome end
                    if end >= chromsize {
                        regions.push(Region {
                            chrom: chrom.clone(),
                            start,
                            end: chromsize,
                            name: format!("{}:{}-{}", chrom, start, end),
                            class: "bin".to_string(),
                        });
                    }
                }
            }
            Err(_e) => {
                panic!("Error reading file {}", file);
            }
        }
    }
    regions
}

pub fn parse_region(reg: String, class: String) -> Vec<Region> {
    let mut regions = Vec::new();
    let sample = reg.clone();

    let reader: Box<dyn BufRead> = match is_gzipped(&reg) {
        Ok(true) => Box::new(BufReader::new(MultiGzDecoder::new(
            File::open(reg).unwrap(),
        ))),
        Ok(false) => Box::new(BufReader::new(File::open(reg).unwrap())),
        Err(e) => {
            panic!("Error reading file {}: {}", reg, e);
        }
    };

    let lines = reader.lines();
    for line in lines {
        match line {
            Ok(line) => {
                let fields: Vec<&str> = line.split('\t').collect();
                let chrom = fields[0].to_string();
                let start = fields[1].to_string();
                let end = fields[2].to_string();
                let name: String = if fields.len() > 3 {
                    fields[3].to_string()
                } else {
                    format!("{}:{}-{}", chrom, start, end)
                };
                let start = start.parse::<u32>().unwrap();
                let end = end.parse::<u32>().unwrap();
                regions.push(Region {
                    chrom,
                    start,
                    end,
                    name,
                    class: class.to_string(),
                });
            }
            Err(_e) => {
                panic!("Error reading file {}", sample);
            }
        }
    }
    regions
}

pub fn is_gzipped(path: &str) -> std::io::Result<bool> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    file.read_exact(&mut magic)?;
    Ok(magic == [0x1F, 0x8B])
}

pub fn decide_methtype(line: Option<String>) -> MethFileType {
    if let Some(l) = line {
        let fields: Vec<&str> = l.split('\t').collect();
        match fields.len() {
            6 => {
                // Either MethylDackel bedgraph or BismarkCov
                // The difference is small:
                // MethylDackel: chrom, start, end, percent, methylcount, coverage
                // BismarkCov: chrom, start, end, percent, methylcount, nonmethylcount
                let perc = fields[3].parse::<f32>().unwrap();
                let meth: u32 = fields[4].parse::<u32>().unwrap();
                let prop_cov: u32 = fields[5].parse::<u32>().unwrap();
                let methyldackel_perc = (meth as f32 / prop_cov as f32) * 100.0;
                let bismarkcov_perc = (meth as f32 / (meth as f32 + prop_cov as f32)) * 100.0;
                if (perc - methyldackel_perc).abs() < 0.01 {
                    return MethFileType::MethylDackel;
                }
                if (perc - bismarkcov_perc).abs() < 0.01 {
                    return MethFileType::BismarkCov;
                }
                panic!(
                    "Could not decide between MethylDackel {} != {}, and BismarkCov {} != {}",
                    perc, methyldackel_perc, perc, bismarkcov_perc
                );
            }
            7 => {
                // Either Allcools tsv of Bismark CpG report
                // The difference here is bit more difficult to detect
                // Allcools: chrom, pos, strand, context, methylcount, coverage, sigtest
                // Bismark CpG report: chrom, pos, strand, methylcount, nonmethylcount, C-context, Trinucleotide context
                // If fields 4 and 5 are both integers, it's allcools file.
                let meth = fields[4].parse::<u32>();
                let cov = fields[5].parse::<u32>();
                if meth.is_ok() && cov.is_ok() {
                    MethFileType::AllCools
                } else {
                    MethFileType::BismarkCpGReport
                }
            }
            18 => {
                // BedMethyl file
                let meth = fields[11].parse::<u32>();
                let cov = fields[9].parse::<u32>();
                if meth.is_ok() && cov.is_ok() {
                    MethFileType::BedMethyl
                } else {
                    panic!(
                        "Could not parse BedMethyl coverage and modification counts (columns 9/11) from: {}",
                        l
                    );
                }
            }
            _ => panic!(
                "Could not decide methylation filetype, as it has {} columns.",
                fields.len()
            ),
        }
    } else {
        panic!(
            "Could not decide methylation filetype, as the file is empty or only contains comments."
        );
    }
}
