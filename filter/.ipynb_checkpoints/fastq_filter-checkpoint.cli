import argparse
import logging
from .filter import filter_fastq
from typing import Tuple, Union

def parse_gc_bounds(value: str) -> Union[Tuple[float, float], float]:
    """Parse GC bounds from string input."""
    if "," in value:
        lower, upper = map(float, value.split(","))
        return (lower, upper)
    return float(value)

def parse_length_bounds(value: str) -> Union[Tuple[int, int], int]:
    """Parse length bounds from string input."""
    if "," in value:
        lower, upper = map(int, value.split(","))
        return (lower, upper)
    return int(value)

def setup_logging(log_file: str = "fastq_filter.log") -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def main():
    """Command-line interface for FastQ filtering."""
    parser = argparse.ArgumentParser(
        description="Filter FastQ files based on GC content, length, and quality."
    )
    
    parser.add_argument(
        "-i", "--input", 
        required=True,
        help="Input FastQ file path"
    )
    parser.add_argument(
        "-o", "--output", 
        default="",
        help="Output FastQ file name (default: filtered.fastq)"
    )
    parser.add_argument(
        "--gc-bounds",
        type=parse_gc_bounds,
        default="0,100",
        help="GC content bounds as either a single value (upper bound) or comma-separated pair (default: 0,100)"
    )
    parser.add_argument(
        "--length-bounds",
        type=parse_length_bounds,
        default="0,4294967296",  # 2^32
        help="Length bounds as either a single value (upper bound) or comma-separated pair (default: 0,4294967296)"
    )
    parser.add_argument(
        "-q", "--quality-threshold",
        type=float,
        default=0,
        help="Minimum average quality score (default: 0)"
    )
    parser.add_argument(
        "--log-file",
        default="fastq_filter.log",
        help="Log file path (default: fastq_filter.log)"
    )
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log_file)
    
    # Run filtering
    filter_fastq(
        input_fastq=args.input,
        output_fastq=args.output,
        gc_bounds=args.gc_bounds,
        length_bounds=args.length_bounds,
        quality_threshold=args.quality_threshold
    )

if __name__ == "__main__":
    main()
