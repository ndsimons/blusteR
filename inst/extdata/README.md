# Reference data for blusteR

This directory holds the **pre-built reference objects** consumed by the
package at run time:

- `bluster_reference_human.rds`
- `bluster_reference_mouse.rds`

Each file is a serialized object of class `bluster_reference` (OAS
background sequences, IEDB epitopes, SAbDab structures, and the
pre-computed background k-mer frequency tables).

## How these files are generated

They are produced by the stand-alone developer script:

```sh
Rscript data-raw/build_references.R
```

See `data-raw/build_references.R` for options (`--species`, `--n`,
`--force`, `--cache`).

The package loads these files automatically via
`system.file("extdata", "bluster_reference_<species>.rds", package = "blusteR")`,
which avoids re-downloading source data and re-computing k-mer
frequencies on every run.
