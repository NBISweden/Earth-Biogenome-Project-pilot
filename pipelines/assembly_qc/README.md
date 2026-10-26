# Assembly QC

K-mer consistency (MerquryFK, Merqury), gene space completeness (BUSCO) and contiguity (gfastats,
QUAST) of one or more assemblies of a sample. It runs on its own, and the same subworkflows are the
evaluation step of the EBP assembly workflow.

## Run

Write a parameter file and pass it with `-params-file`:

```yaml
sample:
  name: 'Species name'
  genome_size: 176040000   # haploid, in bases
  ploidy: 2
assembly:
  - assembler: hifiasm
    stage: scaffolded      # raw, decontaminated, purged, polished, scaffolded or curated
    id: opts01
    pri_fasta: '/path/to/primary.fasta'
    alt_fasta: '/path/to/alternate.fasta'   # optional
hifi:                      # HiFi reads, for the k-mer databases
  - reads: '/path/to/reads.fastq.gz'
busco:
  lineages: 'metazoa_odb10,eukaryota_odb10'
  # lineages_db_path: '/path/to/busco/lineages'
# fastk:
#   kmer_size: 31
# use_phased: true         # evaluate the primary and alternate assemblies separately
```

```bash
nextflow run pipelines/assembly_qc -profile singularity -params-file qc_parameters.yml --outdir results
```

The BUSCO lineages are not looked up from the sample: `busco.lineages` is required.
Results are published under `<outdir>/assembly_qc` and, per tool, under the stage folders
(`<outdir>/<stage>/busco`, `gfastats`, `merquryfk`).

## Test

From the repository root, a stub run on the tiny data set:

```bash
nextflow run pipelines/assembly_qc -profile docker,test -stub-run -params-file pipelines/assembly_qc/tests/params.yml
```
