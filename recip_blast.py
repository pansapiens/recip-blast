#!/usr/bin/env python

import os
import subprocess
import csv
import argparse
from typing import List, Dict, Tuple


def infer_dbtype(fasta_path: str, line_limit: int = 1000) -> str:
    """
    Infer the dbtype (nucl or prot) from the first non-header lines of a FASTA file.
    """
    with open(fasta_path, "r") as fasta_file:
        non_header_lines = 0
        for line in fasta_file:
            if not line.startswith(">"):  # Ignore header lines
                # Check for amino acid-specific characters
                if any(c in line for c in "DEFHIKLMNPQRSVWY"):
                    return "prot"
                non_header_lines += 1
                if non_header_lines >= line_limit:
                    break
    return "nucl"


def create_blast_db(
    fasta_path: str, db_name: str, db_type: str, output_dir: str
) -> None:
    db_path = os.path.join(output_dir, db_name)
    subprocess.run(
        ["makeblastdb", "-in", fasta_path, "-dbtype", db_type, "-out", db_path],
        check=True,
    )


def run_blast(
    query: str,
    db: str,
    out_file: str,
    db_type: str,
    output_dir: str,
    num_threads: int = 1,
    outfmt_str: str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
) -> None:
    db_path = os.path.join(output_dir, db) if output_dir else db
    blast_cmd = "blastn" if db_type == "nucl" else "blastp"
    subprocess.run(
        [
            blast_cmd,
            "-query",
            query,
            "-db",
            db_path,
            "-outfmt",
            outfmt_str,
            "-out",
            out_file,
            "-num_threads",
            str(num_threads),
        ],
        check=True,
    )


def parse_blast_output(file_path: str) -> Dict[str, List[Tuple[str, float, float]]]:
    results = {}
    with open(file_path, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            (
                query,
                subject,
                identity,
                _,
                _,
                _,
                qstart,
                qend,
                sstart,
                send,
                _,
                _,
                qlen,
                slen,
            ) = row[:14]
            identity = float(identity)
            coverage = max(
                int(qend) - int(qstart) + 1, int(send) - int(sstart) + 1
            ) / float(qlen)
            results.setdefault(query, []).append((subject, identity, coverage))
    return results


def filter_hits(
    blast_results: Dict[str, List[Tuple[str, float, float]]],
    identity_threshold: float,
    coverage_threshold: float,
) -> Dict[str, str]:
    best_hits = {}
    for query, hits in blast_results.items():
        for subject, identity, coverage in hits:
            if identity >= identity_threshold and coverage >= coverage_threshold:
                if query not in best_hits or best_hits[query][1] < identity:
                    best_hits[query] = (subject, identity)
    return best_hits


def find_reciprocal_best_hits(
    file1: str, file2: str, identity_threshold: float, coverage_threshold: float
) -> List[Tuple[str, str, float, float]]:
    # Parse the BLAST outputs
    results_file1 = parse_blast_output(file1)
    results_file2 = parse_blast_output(file2)

    # Filter hits above the identity and coverage thresholds
    best_hits_file1 = filter_hits(results_file1, identity_threshold, coverage_threshold)
    best_hits_file2 = filter_hits(results_file2, identity_threshold, coverage_threshold)

    # Find reciprocal best hits
    reciprocal_hits = []
    for query, (subject, identity) in best_hits_file1.items():
        if subject in best_hits_file2 and best_hits_file2[subject][0] == query:
            # Retrieve the coverage for the hit from the original results
            coverage = next((cov for subj, iden, cov in results_file1[query] if subj == subject), 0)
            reciprocal_hits.append((query, subject, identity, coverage))

    return reciprocal_hits



def ensure_directory_exists(directory: str) -> None:
    if not os.path.exists(directory):
        os.makedirs(directory)


def main(args: argparse.Namespace) -> None:
    # Ensure the output directory exists
    ensure_directory_exists(args.output_dir)

    # Infer DB types if not provided
    dbtype_strainA = args.dbtype if args.dbtype else infer_dbtype(args.strainA)
    dbtype_strainB = args.dbtype if args.dbtype else infer_dbtype(args.strainB)

    # Create BLAST databases
    create_blast_db(args.strainA, "strainA_db", dbtype_strainA, args.output_dir)
    create_blast_db(args.strainB, "strainB_db", dbtype_strainB, args.output_dir)

    # Paths for BLAST output
    strainB_vs_strainA_path = os.path.join(args.output_dir, "strainB_vs_strainA.blast")
    strainA_vs_strainB_path = os.path.join(args.output_dir, "strainA_vs_strainB.blast")

    # Run BLAST searches and save the results in the output directory
    run_blast(
        args.strainB,
        "strainA_db",
        strainB_vs_strainA_path,
        dbtype_strainA,
        args.output_dir,
        num_threads=args.threads,
    )
    run_blast(
        args.strainA,
        "strainB_db",
        strainA_vs_strainB_path,
        dbtype_strainB,
        args.output_dir,
        num_threads=args.threads,
    )

    # Find reciprocal best hits
    reciprocal_best_hits = find_reciprocal_best_hits(
        strainB_vs_strainA_path, strainA_vs_strainB_path, args.identity, args.coverage
    )

    # Output results to the specified output directory
    output_file = os.path.join(args.output_dir, args.output)
    with open(output_file, "w") as out_file:
        writer = csv.writer(out_file, delimiter="\t")
        for hit in reciprocal_best_hits:
            writer.writerow(hit)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run reciprocal best hit BLAST search between two strains."
    )
    parser.add_argument("strainA", help="FASTA file of strain A")
    parser.add_argument("strainB", help="FASTA file of strain B")
    parser.add_argument(
        "--identity", type=float, default=90.0, help="Identity threshold for BLAST hits"
    )
    parser.add_argument(
        "--coverage", type=float, default=0.8, help="Coverage threshold for BLAST hits"
    )
    parser.add_argument(
        "--output",
        help="Output filename for reciprocal best hits",
        default="reciprocal_best_hits.tsv",
    )
    parser.add_argument("--output_dir", help="Output directory", default="output")
    parser.add_argument(
        "--dbtype",
        help="Manually specify the BLAST database type (nucl or prot). If not provided, it will be inferred.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for BLAST search (default: 1)",
    )

    args = parser.parse_args()
    main(args)
