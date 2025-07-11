# Configuration

This folder contains workflow configuration files that affect various
aspects of running this pipeline.

These are all in some way included by the `nextflow.config`.

## Test Profile

The test profile is included through a profile called `test`, activated
using `-profile test`. It runs most of the pipeline on a tiny test dataset,
with reduced resources.

## Modules config

The `modules.config` supplies additional parameters to each module, such
as file prefixes, additional tool options, and which files should be published
where. This file extends the `nextflow.config`, supplying module specific
configuration. Resource requirements should not be added here as they will
also affect the test profile.

## Tool Resources config

The `tool_resources.config` file contains the resource requirements for
each tool in the pipeline. It tries to customize the resources as best as
possible to the genome size.