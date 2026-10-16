"""
Bioinformatics utilities for working with FASTQ files.
"""

import helpers.fastq


def filter_fastq(
        input_fastq: str,
        output_fastq: str,
        gc_bounds: tuple | float = (0, 100),
        length_bounds: tuple | int = (0, 2**32),
        quality_threshold: float = 0) -> None:
    """
    Filters FASTQ reads by GC content, length and mean quality.
    Reads are processed one by one, so only one read is kept in memory.

    Arguments:
    - input_fastq (str): path to the input FASTQ file
    - output_fastq (str): output file name; saved to the 'filtered' folder
    - gc_bounds (tuple | float): GC interval in %, (0, 100) by default;
      a single number is treated as the upper bound
    - length_bounds (tuple | int): length interval, (0, 2**32) by default;
      a single number is treated as the upper bound
    - quality_threshold (float): minimal mean quality (phred33), 0 by default

    All intervals include both bounds. Returns None.
    """
    output_path = helpers.fastq.prepare_output_path(output_fastq)

    with open(input_fastq) as input_file:
        with open(output_path, mode='w') as output_file:
            while True:
                record = helpers.fastq.read_fastq_record(input_file)
                if record is None:
                    break
                _, seq, _, quality_string = record

                if (helpers.fastq.passes_gc_filter(seq, gc_bounds)
                        and helpers.fastq.passes_length_filter(
                            seq, length_bounds)
                        and helpers.fastq.passes_quality_filter(
                            quality_string, quality_threshold)):
                    helpers.fastq.write_fastq_record(output_file, record)
