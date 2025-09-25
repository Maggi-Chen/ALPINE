import argparse
import plot_pie_charts

def read_inputs(table_file):
    readtag = open(table_file,'r').read().split('\n')[:-1]
    header = readtag[0].split('\t')
    count = readtag[1].split('\t')
    sample = count[0]

    read_count={}
    for i in range(len(header)-1):
        read_count[header[i+1]] = count[i+1]

    return (sample, read_count, header[1:])

def get_args():
    """ Get command-line arguments """
    parser = argparse.ArgumentParser(description="Merge read classification tables for individual samples",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                     )

    parser.add_argument("-i", "--input", type=str, required=True, nargs='+',
                        help="input 1 or more read count file")
    parser.add_argument("-o","--output", type=str, required=False, default="Merged",
                        help="output filename prefix. [Merged]")
    return parser.parse_args()

def main():
    args = get_args()

    # Count reads in each sample and write into output file
    allcounts = {}
    for input_file in args.input:
        (sample, read_count, header) = read_inputs(input_file)
        allcounts[sample] = read_count

    # Get all read tags for all samples
    alltag=[]
    for c in header:
        alltag+=[c]
    for sample in allcounts:
        for tag in allcounts[sample]:
            if tag not in alltag:
                alltag += [tag]

    # Open output file and write header line
    f=open(args.output+'_read_classification.txt','w')
    f.write('#Sample')
    for tag in alltag:
        f.write('\t'+tag)
    f.write('\n')

    for sample in allcounts:
        f.write(sample)
        for tag in alltag:
            if tag not in allcounts[sample]:
                f.write('\t0')
            else:
                f.write('\t'+str(allcounts[sample][tag]))
        f.write('\n')
    f.close()

    plot_pie_charts.plot_varint_count_analysis(args.output+'_read_classification.txt')

if __name__ == '__main__':
    main()

