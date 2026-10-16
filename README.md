# bioinf-utils

A small package of utilities for working with bioinformatic files.
Currently it provides a FASTQ read filter by GC content, length and quality.

Uses only the Python standard library.

## Installation

Requires Python 3.10 or newer.

```bash
git clone https://github.com/odoaseio/bioinf-utils.git
cd bioinf-utils
```

Import the package from the repository folder:

```python
import bioinf_utils
```

## filter_fastq

Filters reads from a FASTQ file and saves the reads that pass into the
`filtered/` folder (created automatically). Reads are processed one by one,
so no more than one read is kept in memory at a time — suitable for large files.

```python
bioinf_utils.filter_fastq(
    input_fastq='reads.fastq',
    output_fastq='filtered_reads.fastq',
    gc_bounds=(20, 80),
    length_bounds=(50, 100),
    quality_threshold=30,
)
# result: filtered/filtered_reads.fastq
```

| Argument | Description | Default |
|---|---|---|
| `input_fastq` | path to the input FASTQ file | — |
| `output_fastq` | output file name inside `filtered/` | — |
| `gc_bounds` | GC content interval, in % | `(0, 100)` |
| `length_bounds` | read length interval | `(0, 2**32)` |
| `quality_threshold` | minimal mean read quality (phred33 scale) | `0` |

- All intervals include both bounds.
- If a single number is passed to `gc_bounds` or `length_bounds`, it is treated
  as the upper bound: `gc_bounds=44.4` is the same as `(0, 44.4)`.
- Reads with mean quality below `quality_threshold` are discarded.

## Structure

```
bioinf-utils/
|- bioinf_utils.py   # utilities (filter_fastq)
|- helpers/
    |- fastq.py      # FASTQ reading/writing and read checks
```

## Contacts

Kate — [@odoaseio](https://github.com/odoaseio)
