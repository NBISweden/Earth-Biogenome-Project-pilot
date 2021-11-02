# EBP subworkflows

This folder contains subworkflows which can be run either
stand-alone or as part of a larger workflow.

## Available subworkflows

- Assembly Validation: A workflow to check draft assemblies for
various metrics of accuracy.

## Using subworkflows

An example using the assembly validation workflow.
```bash
# nextflow run <repository> -main-script <path/to/script>
nextflow run NBISweden/Earth-Biogenome-Project-pilot \
    -main-script subworkflows/assembly_validation/assembly_validation.nf
```