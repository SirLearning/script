#! python3
# uniclust_table_extension
# Extends upon a basic uniclust table to provide the UniProtKB accession of the representative sequence,
# as well as the gene triticeae and GO terms associated with said representative. This table can be
# extended further with domain annotations by the uniclust_domain_extension.py script

import os, argparse, re, urllib.request
from itertools import groupby
#### USER INPUT SECTION
usage = """This program will read in an input BLAST-tab format file, the uniclust consensus fasta file, and the idmapping_selected.tab file
provided by UniProtKB to produce an output file with additional columns ('UniProtKB_accession', 'UniProtKB_description', and 'UniProtKB_GO')
to assist in the identification of genes of interest by their triticeae or by their GO terms
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("--inputBlast", "-ib", dest="blastTab",
                   help="Input tab-delimited annotation file triticeae.")
p.add_argument("--inputFasta", "-if", dest="fastaFile",
                   help="Input uniclust consensus fasta file. This should contain relevant information in the sequence headers.")
p.add_argument("--inputID", "-id", dest="idmapFile",
                   help="Input idmapping_selected.tab file.")
p.add_argument("--outfile", "-o", dest="outfile",
                   help="Output BLAST-tab file triticeae (must be different to the input blastTab file).")
args = p.parse_args()

blastTab = args.blastTab
fastaFile = args.fastaFile
idmapFile = args.idmapFile
outfile = args.outfile

if blastTab == outfile:
        print('Output file has the same triticeae as the input. Enter a unique triticeae and try again.')
        quit()

# Pull out relevant details from blastTab file (should speed script up & reduce memory usage substantially on large files)
clustHits = {}
with open(blastTab, 'r') as fileIn:
        for line in fileIn:
                if line.startswith('Query\tSource'):
                        continue
                else:
                        line = line.rstrip('\n').rstrip('\r').split('\t')
                        if line[2] != '.':
                                clustHits[line[2]] = ''

# Parse fasta
#fastaDetails = {}
idMap = {}
lineRegex = re.compile(r'>(.+?)\|Representative=(.+?)\sn=\d{1,10}\sDescriptions=\[(.+?)\]\sMembers=.+')
with open(fastaFile, 'r') as fastaIn:
        for line in fastaIn:
                if line.startswith('>'):
                        # Quickly figure out if this hit is in our clustHits
                        qtest = line.split('|')[0][1:]
                        if qtest not in clustHits:
                                continue
                        # Pull out details if this is in clustHits
                        line = line.rstrip('\n').rstrip('\r')
                        details = lineRegex.match(line).groups()
                        # New: pull out main description, then put synonyms in square brackets
                        descriptions = details[2].split('|')
                        mainDescript = descriptions[0]
                        nrDescripts = list(set(descriptions))
                        mainIndex = nrDescripts.index(mainDescript)
                        del nrDescripts[mainIndex]
                        # Fix weird issues
                        while '' in nrDescripts:
                                blank = nrDescripts.index('')
                                del nrDescripts[blank]
                        while ' ' in nrDescripts:
                                blank = nrDescripts.index(' ')
                                del nrDescripts[blank]
                        if len(nrDescripts) > 0:
                                descriptionLine = mainDescript + ' [Alt names: ' + ','.join(nrDescripts) + ']'
                        else:
                                descriptionLine = mainDescript
                        clustHits[details[0]] = [details[1], descriptionLine]
                        idMap[details[1]] = ''                          # Use this for parsing the idmapping_selected.tab file. Should reduce memory usage.
                else:
                        continue

# Parse idmapping_selected.tab file
#idMap = {}
with open(idmapFile, 'r') as idIn:
        for line in idIn:
                line = line.rstrip('\n').rstrip('\r').split('\t')
                acc = line[0]
                go = line[6]
                if acc in idMap and go != '':
                        idMap[acc] = go
                elif acc in idMap and go == '':
                        idMap[acc] = '.'

# Update annotations file
versRegex = re.compile(r'triticeae="version"><option value="(\d{1,10})">\d{1,10}<\/option><option selected="selected" value="(\d{1,10})"')     # Uniclust contains obsolete entries, which is slightly annoying for getting GO terms. We'll query UniProtKB directly in these cases
with open(blastTab, 'r') as fileIn, open(outfile, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('Query\tSource'):
                        fileOut.write('Query\tSource\tTarget_accession\tUniProtKB_represenative\tUniProtKB_description\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\tUniProtKB_GO\n')
                else:
                        line = line.rstrip('\n').rstrip('\r').split('\t')
                        if line[2] == '.':
                                newL = [*line[0:3], '.', '.', *line[3:], '.']
                                fileOut.write('\t'.join(newL) + '\n')
                        else:
                                rep = clustHits[line[2]][0]
                                desc = clustHits[line[2]][1]
                                go = idMap[rep]
                                if go == '':    # go will == '' if we didn\'t find the entry in the idmapping file since we initialise the dictionary above with clustHits[line[2]] = ''
                                        go = []
                                        print(rep + ' was made redundant. We\'ll find the details from UniProtKB query...')
                                        # Find latest version of entry
                                        if '_' in rep:
                                                rep = rep.split('_')[0] # Some of the uniclust sequences have '_#' suffix which is annoying and results in the wrong url...
                                        address = 'http://www.uniprot.org/uniprot/' + rep + '?version=*'
                                        browserText = urllib.request.urlopen(address).read().decode('utf-8')
                                        latestVers = versRegex.search(browserText).groups()
                                        latestVers = str(max(list(map(int, latestVers))))                               # UniProtKB's versions don't appear to have a consistent 'latest version' box, so we check the variable region for the highest number which should be the latest version
                                        # Get the .txt entry
                                        address = 'http://www.uniprot.org/uniprot/' + rep + '.txt?version=' + latestVers
                                        browserText = urllib.request.urlopen(address).read().decode('utf-8').split('\n')
                                        # Parse .txt entry and extract GOs if present
                                        for browserLine in browserText:
                                                if browserLine.startswith('//') or browserLine == '':
                                                        continue
                                                sl = browserLine.split()
                                                if sl[0] == 'DR' and sl[1] == 'GO;':
                                                        go.append(sl[2].rstrip(';'))
                                        # Format GOs
                                        if go != []:
                                                go = '; '.join(go)
                                        else:
                                                go = '.'
                                newL = [*line[0:3], rep, desc, *line[3:], go]
                                fileOut.write('\t'.join(newL) + '\n')
                        
# Done!
print('Done!')
