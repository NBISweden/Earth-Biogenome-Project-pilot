# Test data set provenance

## Authors

- Guilherme Dias
- Martin Pippel
- Cormac Kinsella

## Provenance

* **PB_HIFI/dmel_2Mb.fasta.gz**: Subsampled PacBio Hifi data containing nuclear and mitochondrial reads
* **HIC/*.fastq.gz**: Subsampled Illumina Hi-C data containing only nuclear reads
* **assembly/*.fasta**: Hifiasm nuclear genome assemblies from the above test data. 

### Nuclear genome data

- species: _Drosophila melanogaster_
- PacBio HiFi data was downloaded from NBCI's SRA [SRR10238607](https://www.ncbi.nlm.nih.gov/sra/SRR10238607)
- Illumina HiC data was downloaded from NBCI's SRA [SRR10512944](https://www.ncbi.nlm.nih.gov/sra/SRR10512944)

- PacBio HiFi data were randomly subsampled to 15X read coverage and assembled with [hifiasm version 0.19.8](https://github.com/chhylp123/hifiasm)
- PacBio reads were mapped back to the genome assembly with minimap2
- Illumina HiC reads were mapped back to the genome assembly with bwa mem

- PacBio and HiC reads were extracted from the alignment bam files:
   - from a 2Mb region within a 28Mb contig
   - a 20Kb gap was introduces at Position 1Mb and all PacBio and HiC reads that mapped into that region were filtered out
   - HiC reads were further filtered down to 50%, representing a 38X coverage

### Mitochondrial genome data

- A subset of PacBio HiFi reads (10k) from SRA project [SRR10238607](https://www.ncbi.nlm.nih.gov/sra/SRR10238607) were downloaded
- These were aligned with minimap2 to the _Drosophila melanogaster_ mitochondrial genome, reference sequence PP764103.1
- The first thirty aligned reads were extracted and run through the mitohifi workflow against reference PP764103.1, which retained 6 reads passing all filters
- Filtered reads were added to the existing PacBio HiFi read set and the organelle assembly subworkflow was successfully tested in both contig- and read-based assembly input modes
