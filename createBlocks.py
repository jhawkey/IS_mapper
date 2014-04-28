from argparse import (ArgumentParser, FileType)

def parse_args():

    parser = ArgumentParser(description="concatenate tables file to make blocks file")
    parser.add_argument('--table', nrags='+', type=str, required=True, help='table files to plot blocks from')
    parser.add_argument('--output', type=str, required=True, help='output file name')
    return parser.parse_args()

def main():

    args = parse_args()

    out = open(args.output, 'w')

    for i in args.table:

        header = i.split('_table.txt')[0]
        f = open(i, 'r')
        starts = []
        ends = []
        for columns in (raw.strip().split() for raw in f):
            try:
                starts.append(columns[3])
                ends.append(columns[4])
            except IndexError:
                pass

        for coord in range(0, len(starts)):
            out.write(header + '\t' + starts[coord] + '\t' + ends[coord] + '\n')

    out.close()

if __name__ == '__main__':
    main()