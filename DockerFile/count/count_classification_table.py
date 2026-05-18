import argparse
import sys
import os

# Add shared module path for imports (handles both local and Docker environments)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'shared'))
sys.path.insert(0, '/opt')  # Docker environment path

import plot_pie_charts
from column_order import sort_columns_with_counts, COLUMN_ORDER, get_sort_key


def count_read_group(table_file):
    readtag = open(table_file,'r').read().split('\n')[:-1]
    readtag = [c.split('\t')[1] for c in readtag]

    read_count = {}
    for row in readtag:
        if row not in read_count:
            read_count[row] = 1
        else:
            read_count[row] += 1

    return read_count

def get_args():
    """ Get command-line arguments """
    parser = argparse.ArgumentParser(description="Count read classification table",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                     )

    parser.add_argument("-i", "--input", type=str, required=True,
                        help="input read classification file.")
    parser.add_argument("-o","--output", type=str, required=True,
                        help="output filename prefix.")
    return parser.parse_args()

def main():
    args = get_args()

    # Count reads in classification table
    sample_name = args.input.split('/')[-1].split('readname_')[1][:-4]
    read_count = count_read_group(args.input)

    # Sort columns in fixed order and add missing categories
    ordered_columns = sort_columns_with_counts(read_count)

    # Write into output table
    header = '\t'.join(col for col, _ in ordered_columns)
    counts = '\t'.join(str(cnt) for _, cnt in ordered_columns)

    f = open(args.output + '_read_classification.txt', 'w')
    f.write('#Sample\t' + header + '\n')
    f.write(sample_name + '\t' + counts + '\n')
    f.close()

    plot_pie_charts.plot_varint_count_analysis(args.output + '_read_classification.txt')

if __name__ == '__main__':
    main()

