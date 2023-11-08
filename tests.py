import unittest
from typing import Dict, List, Tuple
import tempfile
import csv
import subprocess
import os

from recip_blast import (
    find_reciprocal_best_hits,
    create_blast_db,
    run_blast,
)

# Assume the following functions are defined as provided earlier:
# parse_blast_output, filter_hits, find_reciprocal_best_hits


class TestReciprocalBestHits(unittest.TestCase):
    def test_find_reciprocal_best_hits(self):
        # Simulate BLAST output data
        blast_results_file1_content = [
            [
                "geneA_strainA",
                "geneX_strainB",
                "95.0",
                "100",
                "5",
                "0",
                "1",
                "100",
                "1",
                "100",
                "0.0",
                "200",
                "100",
                "100",
            ],
            [
                "geneA_strainA",
                "geneY_strainB",
                "85.0",
                "100",
                "15",
                "0",
                "1",
                "100",
                "1",
                "100",
                "0.1",
                "180",
                "100",
                "100",
            ],  # This is a 'miss' due to low identity
        ]
        blast_results_file2_content = [
            [
                "geneX_strainB",
                "geneA_strainA",
                "95.0",
                "100",
                "5",
                "0",
                "1",
                "100",
                "1",
                "100",
                "0.0",
                "200",
                "100",
                "100",
            ],
            [
                "geneY_strainB",
                "geneZ_strainA",
                "85.0",
                "100",
                "15",
                "0",
                "1",
                "100",
                "1",
                "100",
                "0.1",
                "180",
                "100",
                "100",
            ],  # This is a 'miss' because it doesn't match geneA_strainA
        ]

        with tempfile.NamedTemporaryFile(
            mode="w+", delete=False
        ) as file1, tempfile.NamedTemporaryFile(mode="w+", delete=False) as file2:
            writer1 = csv.writer(file1, delimiter="\t")
            writer2 = csv.writer(file2, delimiter="\t")
            writer1.writerows(blast_results_file1_content)
            writer2.writerows(blast_results_file2_content)

            # Flush the content to ensure it's written before the function reads it
            file1.flush()
            file2.flush()

            # Call the function to test
            reciprocal_hits = find_reciprocal_best_hits(
                file1.name, file2.name, identity_threshold=90.0, coverage_threshold=0.8
            )

        # Check the results
        expected_hits = [("geneA_strainA", "geneX_strainB", 95.0, 1.0)]
        self.assertEqual(reciprocal_hits, expected_hits)

        # No need to explicitly delete the files - they are deleted after the context block


class TestEndToEndBLAST(unittest.TestCase):
    def setUp(self):
        # Create temporary FASTA files with known sequences
        self.fasta_content_strainA = """>geneA_strainA
ATGCGTACGTAGCTAGCTGACTGATGCTGACT
>geneB_strainA
ATGCGTACGTAGCTAGCTGACTGATGCTGACG
"""
        self.fasta_content_strainB = """>geneX_strainB
ATGCGTACGTAGCTAGCTGACTGATGCTGACT
>geneY_strainB
ATGCGTACGTAGCTAGCTGACTGATGCTGAAA
"""

    def test_end_to_end_reciprocal_hits(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            fastaA_path = os.path.join(tmp_dir, "strainA.fasta")
            fastaB_path = os.path.join(tmp_dir, "strainB.fasta")
            with open(fastaA_path, "w") as fastaA, open(fastaB_path, "w") as fastaB:
                fastaA.write(self.fasta_content_strainA)
                fastaB.write(self.fasta_content_strainB)

            # Create BLAST databases in the temporary directory
            db_name_A = os.path.basename(fastaA_path).replace(".fasta", "")
            db_name_B = os.path.basename(fastaB_path).replace(".fasta", "")
            create_blast_db(fastaA_path, db_name_A, "nucl", tmp_dir)
            create_blast_db(fastaB_path, db_name_B, "nucl", tmp_dir)

            # Run BLAST
            blast_resultA = os.path.join(
                tmp_dir, db_name_A + "_vs_" + db_name_B + ".blast"
            )
            blast_resultB = os.path.join(
                tmp_dir, db_name_B + "_vs_" + db_name_A + ".blast"
            )
            run_blast(
                fastaA_path, db_name_B, blast_resultA, "nucl", tmp_dir, num_threads=1
            )
            run_blast(
                fastaB_path, db_name_A, blast_resultB, "nucl", tmp_dir, num_threads=1
            )

            # Find reciprocal best hits
            reciprocal_hits = find_reciprocal_best_hits(
                blast_resultA, blast_resultB, 90, 0.8
            )

            # Check the results
            expected_hits = [
                ("geneA_strainA", "geneX_strainB", 100.0, 1.0)
            ]  # Assuming that these are reciprocal best hits
            self.assertEqual(reciprocal_hits, expected_hits)


if __name__ == "__main__":
    unittest.main()
