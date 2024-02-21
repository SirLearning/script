with open('chr1AN.fa', 'r') as input_file, open('chr1AN_modified.fa', 'w') as output_file:
    for line_number, line in enumerate(input_file, start=1):
        if line.startswith('>'):
            line = line.replace(r':\d+-\d+', f':{line_number}')
        output_file.write(line)