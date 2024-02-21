with open('chr1AN.fa', 'r') as infile, open('chr1AN_modified.fa', 'w') as outfile:
    line_number = 1
    for line in infile:
        if line.startswith('>'):
            line = line.replace(':\w+-\w+', f':{line_number}')
            line_number += 1
        outfile.write(line)