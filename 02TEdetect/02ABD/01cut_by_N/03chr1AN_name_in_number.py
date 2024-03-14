import re
import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]

with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
    line_number = 1
    for line in input_file:
        if line.startswith('>'):
            line = re.sub(r":\w+-\d+", f':{line_number}', line)
            line_number += 1
        output_file.write(line)