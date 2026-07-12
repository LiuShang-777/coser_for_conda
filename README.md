# CoSEr

CoSEr calculates coupling scores between two traits from a phenotype table.
It provides two modes:

- `angle`: regression-based angle score, `1 - sin(theta)`.
- `distance`: distance from the shifted regression line, optionally normalized.

## Installation

Local conda build:

```bash
conda build .
conda install --use-local coser
```

## Usage

Angle-based score:

```bash
coser angle -i input.csv -x TraitA -y TraitB -o angle_cs.csv
```

Distance-based score:

```bash
coser distance -i input.csv -x TraitA -y TraitB -o distance_cs.csv
```

For tab-delimited input:

```bash
coser angle -i input.tsv -x TraitA -y TraitB -o angle_cs.csv --sep "\t"
```

If the input file has no sample index column:

```bash
coser angle -i input.csv -x TraitA -y TraitB -o angle_cs.csv --index-col none
```
