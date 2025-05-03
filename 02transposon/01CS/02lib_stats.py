import sys

lib_filename = sys.argv[1]

with open(lib_filename, 'r') as lib_file:
    lib = {}
    for line in lib_file:
        element = line.split('\t')
        chrom = element[0]
        start = int(element[1])
        end = int(element[2])
        numN = end - start - 1
        lib[chrom] = []
        lib[chrom].append((start, end, numN))