"""
Helper functions for filter_fastq: reading, writing and read filters.
"""

import os


def is_in_bounds(value: float, bounds: tuple | float) -> bool:
    """
    Checks whether value lies within bounds (both ends included).

    Arguments:
    - value (float): number to check
    - bounds (tuple | float): (lower, upper) interval or a single number,
      which is treated as the upper bound (lower bound is 0)

    Returns True if lower <= value <= upper, else False.
    """
    if isinstance(bounds, tuple):
        lower_bound, upper_bound = bounds
    else:
        lower_bound, upper_bound = 0, bounds

    return lower_bound <= value <= upper_bound


def passes_gc_filter(seq: str, gc_bounds: tuple | float) -> bool:
    """
    Checks whether GC content of a sequence (in %) lies within gc_bounds.

    Arguments:
    - seq (str): nucleotide sequence (may end with a newline)
    - gc_bounds (tuple | float): GC interval in %, see is_in_bounds

    Returns bool. Empty sequence never passes.
    """
    seq = seq.strip().upper()
    if not seq:
        return False
    gc_percent = (seq.count('C') + seq.count('G')) / len(seq) * 100
    return is_in_bounds(gc_percent, gc_bounds)


def passes_length_filter(seq: str, length_bounds: tuple | int) -> bool:
    """
    Checks whether sequence length lies within length_bounds.

    Arguments:
    - seq (str): nucleotide sequence (may end with a newline)
    - length_bounds (tuple | int): length interval, see is_in_bounds

    Returns bool.
    """
    return is_in_bounds(len(seq.strip()), length_bounds)


def passes_quality_filter(quality_string: str,
                          quality_threshold: float) -> bool:
    """
    Checks whether mean read quality (phred33) is not below the threshold.

    Arguments:
    - quality_string (str): quality line of a FASTQ record
      (may end with a newline)
    - quality_threshold (float): minimal allowed mean quality

    Returns bool. Empty quality string never passes.
    """
    quality_string = quality_string.strip()
    if not quality_string:
        return False

    quality_sum = 0
    for char in quality_string:
        quality_sum += ord(char) - 33
    quality_mean = quality_sum / len(quality_string)

    return quality_mean >= quality_threshold


def read_fastq_record(input_file) -> tuple[str, str, str, str] | None:
    """
    Reads one FASTQ record (4 lines) from an opened file.

    Arguments:
    - input_file: file object opened for reading

    Returns tuple (header, seq, plus_line, quality_string) with lines
    kept as is (with newlines), or None if the end of file is reached.
    """
    header = input_file.readline()
    if not header:
        return None
    seq = input_file.readline()
    plus_line = input_file.readline()
    quality_string = input_file.readline()

    return header, seq, plus_line, quality_string


def write_fastq_record(output_file, record: tuple[str, str, str, str]) -> None:
    """
    Writes one FASTQ record to an opened file. Lines are written as is.

    Arguments:
    - output_file: file object opened for writing
    - record (tuple): (header, seq, plus_line, quality_string)
    """
    for line in record:
        output_file.write(line)


def prepare_output_path(output_fastq: str,
                        output_dir: str = 'filtered') -> str:
    """
    Creates output_dir if it does not exist and builds path to the output file.

    Arguments:
    - output_fastq (str): output file name
    - output_dir (str): output directory, 'filtered' by default

    Returns str: path to the output file.
    """
    os.makedirs(output_dir, exist_ok=True)
    return os.path.join(output_dir, output_fastq)
