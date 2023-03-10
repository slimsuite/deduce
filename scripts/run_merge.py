import argparse
import seaborn as sns
import numpy as np
from deduce_uces.io.sam import parse_sam_lazy
from deduce_uces.stages.candidates import do_merge, calculate_max_gaps


def main():
    parser = argparse.ArgumentParser(
        description="Run the merging algorithm on a sorted BAM file"
    )

    parser.add_argument("bam", metavar="BAM", type=str)
    parser.add_argument("--uce-min-length", type=int, default=200)
    parser.add_argument("--uce-min-homology", type=int, default=100)
    parser.add_argument("--non-unique-core-kmer-buffer", type=int, default=0)

    args = parser.parse_args()

    uce_max_gaps = calculate_max_gaps(
        args.uce_min_length, args.uce_min_homology, args.non_unique_core_kmer_buffer
    )

    print("uce_max_gaps", uce_max_gaps)

    generator = parse_sam_lazy(args.bam)

    n_merged = 0
    lengths = []
    for merged in do_merge(generator, args.uce_min_length, uce_max_gaps):
        n_merged += 1
        lengths.append(len(merged.sequence))

    print("# merged: ", n_merged)
    figure = sns.histplot(data=np.array(lengths)).get_figure()
    figure.savefig("run_lengths.png")


if __name__ == "__main__":
    main()
