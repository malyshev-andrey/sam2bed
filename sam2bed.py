#!/usr/bin/env python3
import argparse
import csv
import re
import sys

# the names of sam columns that must be present
SAM_COLUMNS = [
    'QNAME',
    'FLAG',
    'RNAME',
    'POS',
    'MAPQ',
    'CIGAR',
    'RNEXT',
    'PNEXT',
    'TLEN',
    'SEQ',
    'QUAL'
]

# output bed file columns
OUTPUT_COLUMNS = [
    'RNAME',
    'start',
    'end',
    'QNAME',
    'MAPQ',
    'strand',
    'FLAG',
    'cigar',
    'RNEXT',
    'PNEXT',
    'TLEN',
    'seq',
    'qual',
    'attributes'
]

# sam numeric columns (can be converted to int)
NUMERIC_COLUMNS = {'FLAG', 'POS', 'MAPQ', 'PNEXT', 'TLEN'}

# cigar operations
OPERATIONS = {}
for cigar_operation in {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}:
    OPERATIONS[cigar_operation] = {
        'consumes_query': cigar_operation in {'M', 'I', 'S', '=', 'X'},
        'consumes_reference': cigar_operation in {'M', 'D', 'N', '=', 'X'}
    }
OPERATION_PATTERN = re.compile(r'(\d+)(\D)')


def bed_line_blank(sam_entry: dict, **kwargs) -> dict:
    result = {}
    for column in OUTPUT_COLUMNS:
        if column in sam_entry:
            result[column] = sam_entry[column]
        else:
            result[column] = [] if column in {'seq', 'qual', 'cigar'} else None
    result['strand'] = '-' if sam_entry['FLAG'] % 32 // 16 else '+'
    result.update(kwargs)
    return result


def print_bed_line(line: dict, reduced: bool, output_bed_file) -> None:
    for column in line:
        if column in {'seq', 'qual', 'cigar'}:
            line[column] = ''.join(line[column])
    line['attributes'] = '\t'.join(line['attributes'])
    ncol = 6 if reduced else len(OUTPUT_COLUMNS)
    result = (line[column] for column in OUTPUT_COLUMNS[:ncol])
    print(*result, sep='\t', file=output_bed_file)


def main():
    # CLI arguments parsing
    parser = argparse.ArgumentParser(
        prog='sam2bed',
        description="""
        Converts sam to bed, translating 1-based coordinates of reads to 0-based
        bed intervals.
        """
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input sam file path')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output bed file path')
    parser.add_argument('-s', '--split',
                        help='Split reads by Ns in CIGAR strings',
                        action='store_true')
    parser.add_argument('-S', '--split-with-deletions',
                        help='Split reads by N and D in CIGAR strings',
                        action='store_true')
    parser.add_argument('--reduced',
                        help='Only six output columns',
                        action='store_true')

    args = parser.parse_args()

    # set of operations that will be used to split CIGAR strings
    if args.split_with_deletions:
        split_operations = set('ND')
    elif args.split:
        split_operations = set('N')
    else:
        split_operations = set()

    # input
    if args.input == 'stdin':
        input_file = sys.stdin
    else:
        input_file = open(args.input, mode='rt', newline='')

    # output
    if args.output == 'stdout':
        output_file = sys.stdout
    else:
        output_file = open(args.output, mode='wt')

    # reading the lines of sam file to dictionaries with specified fields
    sam_reader = csv.DictReader(
        filter(lambda s: s[0] != '@', input_file),  # skipping header
        fieldnames=SAM_COLUMNS,
        delimiter='\t',
        restkey='attributes'
    )

    for entry in sam_reader:
        entry: dict
        for field in NUMERIC_COLUMNS:
            entry[field] = int(entry[field])

        # list of all CIGAR operations for current entry
        operations = re.findall(OPERATION_PATTERN, entry['CIGAR'])

        bed_line = bed_line_blank(
            entry,
            start=entry['POS'] - 1,
            end=entry['POS'] - 1
        )
        query_pointer = 0
        segment_n = 1

        for n, operation in operations:
            n = int(n)

            if OPERATIONS[operation]['consumes_query']:
                bed_line['seq'].append(
                    entry['SEQ'][query_pointer:query_pointer + n]
                )
                bed_line['qual'].append(
                    entry['QUAL'][query_pointer:query_pointer + n]
                )
                query_pointer += n

            if operation in split_operations:
                bed_line['QNAME'] += f'/{segment_n}'
                segment_n += 1
                print_bed_line(bed_line, args.reduced, output_file)
                bed_line = bed_line_blank(
                    entry,
                    start=bed_line['end'] + n,
                    end=bed_line['end'] + n
                )
            else:
                bed_line['cigar'].append(f'{n}{operation}')

                if OPERATIONS[operation]['consumes_reference']:
                    bed_line['end'] += n

        if bed_line['seq']:
            bed_line['QNAME'] += f'/{segment_n}'
            print_bed_line(bed_line, args.reduced, output_file)

    if args.input != 'stdin':
        input_file.close()
    if args.output != 'stdout':
        output_file.close()


if __name__ == '__main__':
    main()
