# nf-cmgg/report: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [../07/2026]

### `Added`
#### New Workflows
- **Targeted** (`workflows/targeted.nf`): Targeted variant analysis workflow with gene-specific configurations
- **Pacvar reporting**: enabled (but no `workflows/pacvar.nf` yet)

#### New Modules
- **cat/fastq**: Fastq concatenation module (replaced merge_reads)

### `Changed`
- Module **hotcount**: Optimized process for counting
- Module **samtools**: Update version from 1.22.1 to 1.23.1
- Added pixi environment for reproducible Python setup
- Updated documentation

### `Removed`
- Deprecated scripts and legacy code paths
- Replaced `merge_reads` with `cat_fastq` for better compatibility

## v1.0.0 - [19/03/2026]

Initial release of nf-cmgg/report, created with the [nf-core](https://nf-co.re/) template.
