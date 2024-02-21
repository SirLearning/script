with open('chr1AN.fa', 'r') as infile, open('chr1AN_modified.fa', 'w') as outfile:
    line_number = 1
    for line in infile:
        line = line.replace(r':\d+-\d+', f':{line_number}')
        line_number += 1
        outfile.write(line)