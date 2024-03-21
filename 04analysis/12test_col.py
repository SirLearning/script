with open('data/test.ln.gff3') as f:
    for line in f:
        part = line.split('\t')
        print(part[9])