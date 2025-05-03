import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]
N_length = int(sys.argv[3]) # 2000时，片段平均长度为2.46M，有241个片段

with open(input_filename) as f:
    elements = {}
    frag_list = []
    count = 0
    for line in f:
        element = line.split('\t')
        chrom = element[0]
        start = int(element[1])
        end = int(element[2])
        numN = end - start - 1
        if numN < N_length:
            continue
        count += 1
        elements[count] = {'chrom': chrom, 'start': start, 'end': end, 'numN': numN}
        if count == 1:
            frag = start
        else:
            frag = start - elements[count-1]['end']
        elements[count]['frag'] = frag
        frag_list.append(frag)
for key in elements:
    print(key, elements[key], sep='\t')
# average length of fragments
average_frag = sum(frag_list) / len(frag_list) if frag_list else 0
print("Average frag:", average_frag)

with open(output_filename, 'w') as f:
    for key in elements:
        f.write(f"{elements[key]['chrom']}\t{elements[key]['start']}\t{elements[key]['end']}\n")

# import sys
#
# input_filename = 'N.bed' # sys.argv[1]
# N_length = 1000 # int(sys.argv[2])
#
# with open(input_filename) as f:
#     elements = {}
#     chr0 = ''
#     count = 0
#     for line in f:
#         element = line.split('\t')
#         chrom = element[0]
#         start = int(element[1])
#         end = int(element[2])
#         numN = end - start - 1
#         if numN < N_length:
#             continue
#         count += 1
#         elements[chrom] = []
#         elements[chrom].append((start, end, numN))
#         if chr0 == '':
#             frag = start
#         else:
#             frag = start - elements[chr0][1]
#         elements[chrom].append(frag)
#         chr0 = chrom
# for chrom in elements:
#     print(chrom, elements[chrom], sep='\t')