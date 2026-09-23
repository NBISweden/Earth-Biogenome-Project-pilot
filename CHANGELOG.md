# Changelog

All notable changes to this pipeline will be documented in this file, going
forward from this release. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/).

## [1.0.0] - 2026-09-23

First tagged release. This formalizes the pipeline's current state rather
than marking a fresh rewrite: NBIS/SciLifeLab already runs this pipeline in
production for de novo genome assembly and curation handoff, ahead of it
having a version number attached.

### Included

- Staged assembly workflow driven by `params.steps`: read inspection,
  assembly (HiFi, with optional Hi-C-guided hifiasm and organelle assembly),
  contamination screening, duplicate purging, error polishing, scaffolding,
  and rapid curation handoff (Pretext/CurationPretext).
- RNAseq alignment as an optional side branch.
- Assembly evaluation (BUSCO, Merqury/MerquryFK, k-mer spectra) and a
  Quarto-based assembly report, with MultiQC-collected logs.
- Resume-friendly design: any stage's output can be supplied directly as
  input, letting a run start partway through.

### Known limitations, tracked for post-1.0 work

- Stage gating is a flat, unordered `params.steps` membership list with
  repeated per-stage boilerplate, and some prerequisites (Hi-C to CRAM
  conversion, FastK/Meryl database builds) run unconditionally regardless of
  which stages are requested - see the architecture refactor epic (#381) and
  ADR 0001 (#378) for the planned replacement.
- The assembly data model doesn't yet cleanly support multiple
  assemblies/haplotypes per sample (#246, #338) - a record-type redesign is
  proposed on #338.
- No CI is currently configured for this repository.

See the [architecture refactor epic](https://github.com/NBISweden/Earth-Biogenome-Project-pilot/issues/381)
for the tracked plan going forward.
