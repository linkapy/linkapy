pub struct Region {
    pub chrom: String,
    pub start: Vec<u32>,
    pub end: Vec<u32>,
    pub name: String,
    pub class: String,
}

#[derive(Debug, PartialEq)]
pub struct MethRegion {
    pub chrom: String,
    pub pos: u32,
    pub meth: u32,
    pub total: u32,
}

#[derive(Debug, PartialEq)]
pub enum MethFileType {
    AllCools,
    MethylDackel,
    BedMethyl,
    BismarkCov,
    BismarkCpGReport,
}

fn split_exact<const N: usize>(line: &str) -> Option<[&str; N]> {
    let mut iter = line.split('\t');
    let mut fields = [""; N];
    for field in fields.iter_mut() {
        *field = iter.next()?;
    }
    if iter.next().is_some() {
        return None;
    }
    Some(fields)
}

impl MethFileType {
    pub fn parse_line(&self, line: &str) -> Result<Option<MethRegion>, String> {
        // Ignore any comment lines (#)
        if line.starts_with('#') {
            return Ok(None);
        };
        match self {
            MethFileType::AllCools => {
                let fields = split_exact::<7>(line)
                    .ok_or_else(|| format!("Invalid AllCools line: {}", line))?;
                let chrom = fields[0].to_string();
                // Allcool files is 1-based. Convert to 0-based.
                let pos = fields[1]
                    .parse::<u32>()
                    .map_err(|e| format!("Invalid position in AllCools line: {}: {}", line, e))?
                    - 1;
                let meth = fields[4].parse::<u32>().map_err(|e| {
                    format!("Invalid methylated count in AllCools line: {}: {}", line, e)
                })?;
                let total = fields[5].parse::<u32>().map_err(|e| {
                    format!("Invalid total count in AllCools line: {}: {}", line, e)
                })?;
                Ok(Some(MethRegion {
                    chrom,
                    pos,
                    meth,
                    total,
                }))
            }
            MethFileType::MethylDackel => {
                let fields = split_exact::<6>(line)
                    .ok_or_else(|| format!("Invalid MethylDackel line: {}", line))?;
                let chrom = fields[0].to_string();
                let pos = fields[1].parse::<u32>().map_err(|e| {
                    format!("Invalid position in MethylDackel line: {}: {}", line, e)
                })?;
                let meth = fields[4].parse::<u32>().map_err(|e| {
                    format!(
                        "Invalid methylated count in MethylDackel line: {}: {}",
                        line, e
                    )
                })?;
                let total = fields[5].parse::<u32>().map_err(|e| {
                    format!("Invalid total count in MethylDackel line: {}: {}", line, e)
                })?;
                Ok(Some(MethRegion {
                    chrom,
                    pos,
                    meth,
                    total,
                }))
            }
            //MethFileType::BismarkBedgraph => {}
            MethFileType::BismarkCov => {
                let fields = split_exact::<6>(line)
                    .ok_or_else(|| format!("Invalid BismarkCov line: {}", line))?;
                let chrom = fields[0].to_string();
                // Bismark Coverage is 1-based by default. Even though 0-based can be specified.
                // We assume it's 1-based.
                let pos = fields[1]
                    .parse::<u32>()
                    .map_err(|e| format!("Invalid position in BismarkCov line: {}: {}", line, e))?
                    - 1;
                let meth = fields[4].parse::<u32>().map_err(|e| {
                    format!(
                        "Invalid methylated count in BismarkCov line: {}: {}",
                        line, e
                    )
                })?;
                let unmeth = fields[5].parse::<u32>().map_err(|e| {
                    format!(
                        "Invalid unmethylated count in BismarkCov line: {}: {}",
                        line, e
                    )
                })?;
                let total = meth + unmeth;
                Ok(Some(MethRegion {
                    chrom,
                    pos,
                    meth,
                    total,
                }))
            }
            MethFileType::BismarkCpGReport => {
                let fields = split_exact::<7>(line)
                    .ok_or_else(|| format!("Invalid BismarkCpGReport line: {}", line))?;
                let chrom = fields[0].to_string();
                // Bismark CpG Report is 0-based.
                let pos = fields[1].parse::<u32>().map_err(|e| {
                    format!("Invalid position in BismarkCpGReport line: {}: {}", line, e)
                })?;
                let meth = fields[3].parse::<u32>().map_err(|e| {
                    format!(
                        "Invalid methylated count in BismarkCpGReport line: {}: {}",
                        line, e
                    )
                })?;
                let unmeth = fields[4].parse::<u32>().map_err(|e| {
                    format!(
                        "Invalid unmethylated count in BismarkCpGReport line: {}: {}",
                        line, e
                    )
                })?;
                if meth != 0 || unmeth != 0 {
                    let total = meth + unmeth;
                    Ok(Some(MethRegion {
                        chrom,
                        pos,
                        meth,
                        total,
                    }))
                } else {
                    Ok(None)
                }
            }
            MethFileType::BedMethyl => {
                let fields = split_exact::<18>(line)
                    .ok_or_else(|| format!("Invalid BedMethyl line: {}", line))?;
                let chrom = fields[0].to_string();
                // BedMethyl is 0-based.
                let pos = fields[1]
                    .parse::<u32>()
                    .map_err(|e| format!("Invalid position in BedMethyl line: {}: {}", line, e))?;
                let meth = fields[11].parse::<u32>().map_err(|e| {
                    format!(
                        "Invalid methylated count in BedMethyl line: {}: {}",
                        line, e
                    )
                })?;
                let total = fields[9].parse::<u32>().map_err(|e| {
                    format!(
                        "Invalid unmethylated count in BedMethyl line: {}: {}",
                        line, e
                    )
                })?;
                Ok(Some(MethRegion {
                    chrom,
                    pos,
                    meth,
                    total,
                }))
            }
        }
    }
}
