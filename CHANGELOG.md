# Changelog

## [1.10.0](https://github.com/linkapy/linkapy/compare/v1.9.0...v1.10.0) (2026-09-06)


### Features

* logger version dump ([#77](https://github.com/linkapy/linkapy/issues/77)) ([971487e](https://github.com/linkapy/linkapy/commit/971487e79f6663550cb3b8af7b4590d9cb5a54b9))
* mqc compatible QC ([18d97e5](https://github.com/linkapy/linkapy/commit/18d97e51b9fed90649a37df648cb75b333e184b9))
* performance gain, drop dense setting ([#80](https://github.com/linkapy/linkapy/issues/80)) ([18d97e5](https://github.com/linkapy/linkapy/commit/18d97e51b9fed90649a37df648cb75b333e184b9))

## [1.9.0](https://github.com/linkapy/linkapy/compare/v1.8.0...v1.9.0) (2026-08-15)


### Features

* exact split for methtypes ([0a0433f](https://github.com/linkapy/linkapy/commit/0a0433f78e8525acfb37b747bc627907e3fa3c1e))
* single gz pass for methtype ([0a0433f](https://github.com/linkapy/linkapy/commit/0a0433f78e8525acfb37b747bc627907e3fa3c1e))


### Bug Fixes

* drop global lock in latest pyo3 ([0a0433f](https://github.com/linkapy/linkapy/commit/0a0433f78e8525acfb37b747bc627907e3fa3c1e))

## [1.8.0](https://github.com/linkapy/linkapy/compare/v1.7.0...v1.8.0) (2026-06-19)


### Features

* unique var names to accomodate pull_on_update behavior in mudata ([#67](https://github.com/linkapy/linkapy/issues/67)) ([d93f4aa](https://github.com/linkapy/linkapy/commit/d93f4aaba9b38d118f1cbad36bc047d537ac9823))

## [1.7.0](https://github.com/linkapy/linkapy/compare/v1.6.0...v1.7.0) (2026-02-06)


### Features

* skip 'track' headers ([43a17a6](https://github.com/linkapy/linkapy/commit/43a17a6d474d7c010d36f234abe2c36c01215f5d))


### Bug Fixes

* full gz file parsing ([#60](https://github.com/linkapy/linkapy/issues/60)) ([43a17a6](https://github.com/linkapy/linkapy/commit/43a17a6d474d7c010d36f234abe2c36c01215f5d))

## [1.6.0](https://github.com/linkapy/linkapy/compare/v1.5.0...v1.6.0) (2026-02-04)


### Features

* more clear logs on decision making for matching cells ([6333e72](https://github.com/linkapy/linkapy/commit/6333e723429253227429e75be08c4fc468e5726a))
* obs_index written as object dtypes instead of strings ([6333e72](https://github.com/linkapy/linkapy/commit/6333e723429253227429e75be08c4fc468e5726a))

## [1.5.0](https://github.com/FunctionalEpigeneticsLab/linkapy/compare/v1.4.0...v1.5.0) (2025-09-19)


### Features

* Bedmethyl support ([#43](https://github.com/FunctionalEpigeneticsLab/linkapy/issues/43)) ([263c64c](https://github.com/FunctionalEpigeneticsLab/linkapy/commit/263c64cc1e8b2f4d45b5af4159d14957ac622fb7))

## [1.4.0](https://github.com/FunctionalEpigeneticsLab/linkapy/compare/v1.3.0...v1.4.0) (2025-09-15)


### Features

* example CLI to get example usage. ([#38](https://github.com/FunctionalEpigeneticsLab/linkapy/issues/38)) ([daee609](https://github.com/FunctionalEpigeneticsLab/linkapy/commit/daee6094341d989fc8749c5ef1636977c9fb5b9a))

## [1.3.0](https://github.com/FunctionalEpigeneticsLab/linkapy/compare/v1.2.1...v1.3.0) (2025-09-12)


### Features

* meth file support for methyldackel, allcools, bismark cov and bismark CpG. Both gz and regular. ([#32](https://github.com/FunctionalEpigeneticsLab/linkapy/issues/32)) ([88efc83](https://github.com/FunctionalEpigeneticsLab/linkapy/commit/88efc83d9372f0475cda6ad663194607751f5c68))

## [1.2.1](https://github.com/WardDeb/linkapy/compare/v1.2.0...v1.2.1) (2025-08-22)


### Bug Fixes

* drop dense pass in mudata dump to avoid loss of NA information ([a94300e](https://github.com/WardDeb/linkapy/commit/a94300e3b45bab70e7dc9c90fa89f52680237fa1))
* frac moved to rust code  ([a94300e](https://github.com/WardDeb/linkapy/commit/a94300e3b45bab70e7dc9c90fa89f52680237fa1))
* Retain zero-nan discrimination in mudata creation ([#23](https://github.com/WardDeb/linkapy/issues/23)) ([a94300e](https://github.com/WardDeb/linkapy/commit/a94300e3b45bab70e7dc9c90fa89f52680237fa1))

## [1.2.0](https://github.com/WardDeb/linkapy/compare/v1.1.0...v1.2.0) (2025-08-21)


### Features

* chromsizes / bins mode included ([e87602a](https://github.com/WardDeb/linkapy/commit/e87602a99a3a7c03c05254dcbdbbdeb765fd18d6))
* CLI ([e87602a](https://github.com/WardDeb/linkapy/commit/e87602a99a3a7c03c05254dcbdbbdeb765fd18d6))
* drop scNMT specificity for more general RNA/Methylation assays ([e87602a](https://github.com/WardDeb/linkapy/commit/e87602a99a3a7c03c05254dcbdbbdeb765fd18d6))
* obs matching feature to identify common cells between assays based on edit distance ([e87602a](https://github.com/WardDeb/linkapy/commit/e87602a99a3a7c03c05254dcbdbbdeb765fd18d6))
* pattern names to rename assays ([e87602a](https://github.com/WardDeb/linkapy/commit/e87602a99a3a7c03c05254dcbdbbdeb765fd18d6))
* verbose mode ([e87602a](https://github.com/WardDeb/linkapy/commit/e87602a99a3a7c03c05254dcbdbbdeb765fd18d6))


### Bug Fixes

* summarization on half-open intervals ([e87602a](https://github.com/WardDeb/linkapy/commit/e87602a99a3a7c03c05254dcbdbbdeb765fd18d6))

## [1.1.0](https://github.com/WardDeb/linkapy/compare/v1.0.0...v1.1.0) (2025-07-29)


### Features

* automated pypi upload on release ([#10](https://github.com/WardDeb/linkapy/issues/10)) ([d190211](https://github.com/WardDeb/linkapy/commit/d190211286d8d5fb0a1c10f86f16b62fee6acda1))
* versioning, libname, templates ([#7](https://github.com/WardDeb/linkapy/issues/7)) ([da6bc41](https://github.com/WardDeb/linkapy/commit/da6bc41b0845a0ab1116cb52dfa4db701d129a96))

## 1.0.0 (2025-07-26)


### Features

* release please ([da22150](https://github.com/WardDeb/linkapy/commit/da221505a6e27dcd75102640e8a4b201d9e1e862))
