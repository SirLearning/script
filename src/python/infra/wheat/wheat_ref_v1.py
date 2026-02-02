import json
import sys
import csv
import argparse

parser = argparse.ArgumentParser(description="Convert chromosome names and positions based on a JSON map.")
parser.add_argument("-m", "--map_json", required=True, help="Path to the JSON mapping file")
parser.add_argument("-i", "--input_file", required=True, help="Path to the input data file")
parser.add_argument("-o", "--output_file", required=True, help="Path to the output file")
parser.add_argument("-c", "--chrom_col", required=True, help="Chromosome column name or 1-based index")
parser.add_argument("-p", "--pos_col", required=True, help="Position column name or 1-based index")

args = parser.parse_args()

map_file = args.map_json
input_file = args.input_file
output_file = args.output_file
chrom_col_arg = args.chrom_col
pos_col_arg = args.pos_col

with open(map_file, 'r') as f:
    mapping = json.load(f)

def get_col_idx(header, arg):
    # Try to match by name first
    if arg in header:
        return header.index(arg)
    # Try 1-based index
    if arg.isdigit():
        return int(arg) - 1 # 1-based arg
    return -1

with open(input_file, 'r', newline='') as fin, open(output_file, 'w', newline='') as fout:
    reader = csv.reader(fin, delimiter='\\t')
    writer = csv.writer(fout, delimiter='\\t')

    try:
        header = next(reader)
    except StopIteration:
        sys.exit(0)
    
    writer.writerow(header)

    c_idx = get_col_idx(header, chrom_col_arg)
    p_idx = get_col_idx(header, pos_col_arg)

    if c_idx == -1 or p_idx == -1:
            sys.stderr.write(f"Error: Columns {chrom_col_arg}, {pos_col_arg} not found in header: {header}\\n")
            sys.exit(1)

    for parts in reader:
        if len(parts) <= max(c_idx, p_idx):
            writer.writerow(parts)
            continue

        chrom = parts[c_idx]
        # Check if chrom is in mapping (e.g. "2")
        if chrom in mapping:
            rule = mapping[chrom]
            # Update Chromosome Name
            parts[c_idx] = rule.get('chrom', chrom)
            
            # Update Position
            try:
                pos = int(parts[p_idx])
                # Calculate offset: support 'offset' key
                offset = rule.get('offset', 0)
                
                parts[p_idx] = str(pos + offset)
            except ValueError:
                pass
        
        writer.writerow(parts)