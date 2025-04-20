import os
import sys
import unittest
import tempfile
import logging

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
                
from filter.filter import filter_fastq

class TestFastQFilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Create test FASTQ files."""
        cls.test_dir = tempfile.mkdtemp()
        cls.good_fastq = os.path.join(cls.test_dir, "test_good.fastq")
        cls.bad_fastq = os.path.join(cls.test_dir, "test_bad.fastq")
        
        # Create good test file with known sequences
        with open(cls.good_fastq, "w") as f:
            f.write("@good_seq\nACGT\n+\nIIII\n")
            f.write("@bad_gc\nAAAA\n+\nIIII\n")
            f.write("@bad_len\nACGTACGT\n+\nIIIIIIII\n")
            f.write("@bad_qual\nACGT\n+\n!!!!\n")
        
        # Create bad test file (empty)
        with open(cls.bad_fastq, "w") as f:
            f.write("")

    def test_gc_filtering(self):
        """Test GC content filtering."""
        output = filter_fastq(
            self.good_fastq,
            gc_bounds=(50, 100),
            output_fastq="gc_filtered.fastq"
        )
        self.assertTrue(os.path.exists(output))
        with open(output) as f:
            content = f.read()
        self.assertIn("@good_seq", content)
        self.assertNotIn("@bad_gc", content)

    def test_length_filtering(self):
        """Test length filtering."""
        output = filter_fastq(
            self.good_fastq,
            length_bounds=(1, 4),
            output_fastq="length_filtered.fastq"
        )
        with open(output) as f:
            content = f.read()
        self.assertIn("@good_seq", content)
        self.assertNotIn("@bad_len", content)

    def test_quality_filtering(self):
        """Test quality filtering."""
        output = filter_fastq(
            self.good_fastq,
            quality_threshold=30,
            output_fastq="qual_filtered.fastq"
        )
        with open(output) as f:
            content = f.read()
        self.assertIn("@good_seq", content)
        self.assertNotIn("@bad_qual", content)

    def test_file_not_found(self):
        """Test handling of missing input file."""
        with self.assertRaises(FileNotFoundError):
            filter_fastq("nonexistent.fastq")

    def test_invalid_bounds(self):
        """Test invalid bounds handling."""
        with self.assertRaises(ValueError):
            filter_fastq(
                self.good_fastq,
                gc_bounds=(100, 0)
            )

    def test_empty_file(self):
        """Test handling of empty input file."""
        output = filter_fastq(
            self.bad_fastq,
            output_fastq="empty_filtered.fastq"
        )
        self.assertTrue(os.path.exists(output))
        with open(output) as f:
            content = f.read()
        self.assertEqual(content, "")

    def test_log_file_creation(self):
        """Test that log file is created."""
        log_file = os.path.join(self.test_dir, "test.log")
        if os.path.exists(log_file):
            os.remove(log_file)
            
        filter_fastq(
            self.good_fastq,
            output_fastq="log_test.fastq"
        )
        self.assertTrue(os.path.exists("fastq_filter.log"))

    def test_single_bound_values(self):
        """Test single value bounds (upper bound only)."""
        output = filter_fastq(
            self.good_fastq,
            gc_bounds=25,  # Only upper bound
            output_fastq="single_bound.fastq"
        )
        with open(output) as f:
            content = f.read()
        self.assertIn("@bad_gc", content)  # GC=0 should pass
        self.assertNotIn("@good_seq", content)  # GC=50 should fail

    @classmethod
    def tearDownClass(cls):
        """Clean up test files."""
        for f in os.listdir(cls.test_dir):
            os.remove(os.path.join(cls.test_dir, f))
        os.rmdir(cls.test_dir)
        # Clean up any created filtered directories
        filtered_dir = os.path.join(os.getcwd(), "filtered")
        if os.path.exists(filtered_dir):
            for f in os.listdir(filtered_dir):
                os.remove(os.path.join(filtered_dir, f))
            os.rmdir(filtered_dir)

if __name__ == "__main__":
    unittest.main()
