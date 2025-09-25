import argparse
import plot_pie_charts

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

    # Write into output table
    header=''
    count=''
    for key in read_count:
        header += '\t'+key
        count += '\t'+str(read_count[key])

    f=open(args.output+'_read_classification.txt','w')
    f.write('#Sample')
    f.write(header+'\n')
    f.write(sample_name)
    f.write(count+'\n')
    f.close()

    plot_pie_charts.plot_varint_count_analysis(args.output+'_read_classification.txt')

if __name__ == '__main__':
    main()

