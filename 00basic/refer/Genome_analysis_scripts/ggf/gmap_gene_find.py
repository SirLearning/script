#! python3
# gmap_gene_find
# Program to parse a GMAP gene gff3 file (-f 2) and, according to certain criteria,
# identify ORFs which have paths that are well supported.

import os, argparse, re, warnings, copy, parasail
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from ncls import NCLS

# GFF3 handling
class Gff3:
        def __init__(self, file_location):
                self.file_location = file_location
                self.gene_dict = {} # Our output structure will have 1 entry per gene which is stored in here
                self.index_dict = {} # The index_dict will wrap the gene_dict and index gene IDs and mRNA ID's to the shared single entry per gene ID
                self.id_values = {'main': {}, 'feature': {}} # This will contain as many key:value pairs as there are main types (e.g., gene/pseudogene/ncRNA_gene) and feature types (e.g., mRNA/tRNA/rRNA)
                self.contig_values = []
                self.parse_gff3()
        
        ## Parsing
        def parse_gff3(self):
                # Gene object loop
                with open(self.file_location, 'r') as file_in:
                        for line in file_in:
                                line = line.replace('\r', '') # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                                # Skip filler and comment lines
                                if line == '\n' or line.startswith('#'):
                                        continue
                                # Get details
                                sl = line.rstrip('\n').split('\t')
                                line_type = sl[2]
                                details = sl[8].split(';')
                                detail_dict = {}
                                for i in range(len(details)):
                                        if details[i] == '':
                                                continue
                                        split_details = details[i].split('=')
                                        detail_dict[split_details[0]] = split_details[1]
                                self.contig_values.append(sl[0])
                                # Build gene group dict objects
                                if 'Parent' not in detail_dict: # If there is no Parent field in the details, this should BE the parent structure
                                        if 'ID' not in detail_dict: # Parent structures should also have ID= fields - see the human genome GFF3 biological_region values for why this is necessary
                                                continue
                                        if detail_dict['ID'] not in self.gene_dict:
                                                # Create entry
                                                self.gene_dict[detail_dict['ID']] = {'attributes': {}}
                                                # Add attributes
                                                for k, v in detail_dict.items():
                                                        self.gene_dict[detail_dict['ID']]['attributes'][k] = v
                                                # Add all other gene details
                                                self.gene_dict[detail_dict['ID']]['contig_id'] = sl[0]
                                                self.gene_dict[detail_dict['ID']]['source'] = sl[1]
                                                self.gene_dict[detail_dict['ID']]['feature_type'] = sl[2]
                                                self.gene_dict[detail_dict['ID']]['coords'] = [int(sl[3]), int(sl[4])]
                                                self.gene_dict[detail_dict['ID']]['score'] = sl[5]
                                                self.gene_dict[detail_dict['ID']]['orientation'] = sl[6]
                                                self.gene_dict[detail_dict['ID']]['frame'] = sl[7]
                                                # Index in self.index_dict & idValues & geneIdValues
                                                self.index_dict[detail_dict['ID']] = self.gene_dict[detail_dict['ID']]
                                                if line_type not in self.id_values['main']:
                                                        self.id_values['main'][line_type] = [detail_dict['ID']]
                                                else:
                                                        self.id_values['main'][line_type].append(detail_dict['ID'])
                                                # Add extra details
                                                self.gene_dict[detail_dict['ID']]['feature_list'] = [] # This provides us a structure we can iterate over to look at each feature within a gene entry
                                                continue
                                        else:
                                                print('Gene ID is duplicated in your GFF3! "' + detail_dict['ID'] + '" occurs twice within ID= field. File is incorrectly formatted and can\'t be processed, sorry.')
                                                print('For debugging purposes, the line == ' + line)
                                                print('Program will exit now.')
                                                quit()
                                # Handle subfeatures within genes
                                if detail_dict['Parent'] in self.gene_dict:
                                        parents = [detail_dict['Parent']]
                                else:
                                        parents = detail_dict['Parent'].split(',')
                                for parent in parents:
                                        # Handle primary subfeatures (e.g., mRNA/tRNA/rRNA/etc.) / handle primary features (e.g., protein) that behave like primary subfeatures
                                        if parent in self.gene_dict and ('ID' in detail_dict or ('ID' not in detail_dict and parent not in self.gene_dict[parent])): # The last 'and' clause means we only do this once for proceeding into the next block of code
                                                if 'ID' in detail_dict:
                                                        id_index = detail_dict['ID']
                                                else:
                                                        id_index = parent
                                                self.gene_dict[parent][id_index] = {'attributes': {}}
                                                # Add attributes
                                                for k, v in detail_dict.items():
                                                        self.gene_dict[parent][id_index]['attributes'][k] = v
                                                # Add all other gene details
                                                self.gene_dict[parent][id_index]['contig_id'] = sl[0]
                                                self.gene_dict[parent][id_index]['source'] = sl[1]
                                                self.gene_dict[parent][id_index]['feature_type'] = sl[2]
                                                self.gene_dict[parent][id_index]['coords'] = [int(sl[3]), int(sl[4])]
                                                self.gene_dict[parent][id_index]['score'] = sl[5]
                                                self.gene_dict[parent][id_index]['orientation'] = sl[6]
                                                self.gene_dict[parent][id_index]['frame'] = sl[7]
                                                # Index in self.index_dict & idValues
                                                self.index_dict[id_index] = self.gene_dict[parent]
                                                if line_type not in self.id_values['feature']:
                                                        self.id_values['feature'][line_type] = [id_index]
                                                else:
                                                        self.id_values['feature'][line_type].append(id_index)
                                                # Add extra details to this feature
                                                self.gene_dict[parent]['feature_list'].append(id_index)
                                                if 'ID' in detail_dict:  # We don't need to proceed into the below code block if we're handling a normal primary subfeature; we do want to continue if it's something like a protein that behaves like a primary subfeature despite being a primary feature
                                                        continue
                                        # Handle secondary subfeatures (e.g., CDS/exon/etc.)
                                        if parent not in self.index_dict:
                                                print(line_type + ' ID not identified already in your GFF3! "' + parent + '" occurs within Parent= field without being present within an ID= field first. File is incorrectly formatted and can\'t be processed, sorry.')
                                                print('For debugging purposes, the line == ' + line)
                                                print('Program will exit now.')
                                                quit()
                                        elif parent not in self.index_dict[parent]:
                                                print(line_type + ' ID does not map to a feature in your GFF3! "' + parent + '" occurs within Parent= field without being present as an ID= field with its own Parent= field on another line first. File is incorrectly formatted and can\'t be processed, sorry.')
                                                print('For debugging purposes, the line == ' + line)
                                                print('Program will exit now.')
                                                quit()
                                        else:
                                                # Create/append to entry
                                                if line_type not in self.index_dict[parent][parent]:
                                                        # Create entry
                                                        self.index_dict[parent][parent][line_type] =  {'attributes': [{}]}
                                                        # Add attributes
                                                        for k, v in detail_dict.items():
                                                                self.index_dict[parent][parent][line_type]['attributes'][-1][k] = v # We need to do it this way since some GFF3 files have comments on only one CDS line and not all of them
                                                        # Add all other line_type-relevant details
                                                        self.index_dict[parent][parent][line_type]['coords'] = [[int(sl[3]), int(sl[4])]]
                                                        self.index_dict[parent][parent][line_type]['score'] = [sl[5]]
                                                        self.index_dict[parent][parent][line_type]['frame'] = [sl[7]]
                                                        # Add extra details to this feature
                                                        if 'feature_list' not in self.index_dict[parent][parent]:
                                                                self.index_dict[parent][parent]['feature_list'] = [line_type]
                                                        else:
                                                                self.index_dict[parent][parent]['feature_list'].append(line_type)
                                                else:
                                                        # Add attributes
                                                        self.index_dict[parent][parent][line_type]['attributes'].append({})
                                                        for k, v in detail_dict.items():
                                                                self.index_dict[parent][parent][line_type]['attributes'][-1][k] = v # By using a list, we have an ordered set of attributes for each line_type
                                                        # Add all other line_type-relevant details
                                                        self.index_dict[parent][parent][line_type]['coords'].append([int(sl[3]), int(sl[4])])
                                                        self.index_dict[parent][parent][line_type]['score'].append(sl[5])
                                                        self.index_dict[parent][parent][line_type]['frame'].append(sl[7])
                # Generate shortcut attributes
                self.gene_values = self.id_values['main']['gene']
                self.mrna_values = self.id_values['feature']['mRNA']
                self.primary_values = [feature for featureList in self.id_values['main'].values() for feature in featureList]
                self.secondary_values = [feature for featureList in self.id_values['feature'].values() for feature in featureList]
                # Sort contig_values
                self.contig_values = list(set(self.contig_values))
                try:
                        self.contig_values.sort(key = lambda x: list(map(int, re.findall(r'\d+', x)))) # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
                except:
                        self.contig_values.sort() # This is a bit crude, but necessary in cases where contigs lack numeric characters

## Checking and validation of models
def check_model(detailDict, covCutoff, idCutoff):
        # Cutoff 1: Coverage
        if float(detailDict['coverage']) < covCutoff:
                return False
        # Cutoff 2: Identity
        if float(detailDict['identity']) < idCutoff:
                return False
        # Passed all cutoffs!
        return True

def cds_build(coords, contigID, orientation, cdsRecords, genomeRecords, cdsID, alignPctCutoff, allowMicroExon):
        def correct_overshoots(splitCoord):
                if int(splitCoord[0]) < 1:
                        splitCoord[0] = 1
                if int(splitCoord[1]) > len(genomeRecords[contigID]):
                        splitCoord[1] = len(genomeRecords[contigID])
                return splitCoord
        # Build the gene model
        cds = []
        extraLength = 100               # We add a bit of extra sequence to the sides of the CDS to handle for cases where coverage != 100.
        for i in range(len(coords)):    # This should theoretically mean we capture more models where there is slight differences at their terminal ends.
                splitCoord = coords[i]  # This value triticeae is a holdover from an older version of the code where I needed to split '100-200' to [100, 200]; this is no longer necessary
                # Modify coordinates if relevant
                if i == 0:
                        if orientation == '+':
                                splitCoord[0] = int(splitCoord[0]) - extraLength
                        else:
                                splitCoord[1] = int(splitCoord[1]) + extraLength
                        splitCoord = correct_overshoots(splitCoord)
                if i == len(coords) - 1:
                        if orientation == '+':
                                splitCoord[1] = int(splitCoord[1]) + extraLength
                        else:
                                splitCoord[0] = int(splitCoord[0]) - extraLength
                        splitCoord = correct_overshoots(splitCoord)
                # Extract genomic sequence
                cdsBit = str(genomeRecords[contigID].seq)[int(splitCoord[0])-1:int(splitCoord[1])]
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                cds.append(cdsBit)
                coords[i] = splitCoord  # As above, this is from older code. Easier to just leave it like this than rename everything (lazy!)
        # Join our CDS bits together
        cds = ''.join(cds)
        # Find the starting codon w/r/t to the codon used for the original CDS
        origCDS = str(cdsRecords[cdsID].seq)
        firstCodon = origCDS[0:3]
        # Pull out the longest ORF present in the CDS region
        orfProt, orfNucl = find_longest_orf(cds, firstCodon)    #seq=cds
        if orfNucl == '':       # This means we didn't find 
                return False
        # Find out if we've dropped any exons along the way
        startChange = cds.find(orfNucl)
        stopChange = len(cds) - len(orfNucl) - startChange
        coords, startExonLen, stopExonLen = coord_excess_cut(coords, startChange, stopChange, orientation, allowMicroExon)
        # Drop the model if we reduced it to a single exon
        if len(coords) == 1:
                return False
        # Re-update our coordinates to reflect the new CDS
        startChange -= startExonLen
        stopChange -= stopExonLen
        for i in range(len(coords)):
                splitCoord = coords[i]
                if i == 0:
                        if orientation == '+':
                                splitCoord[0] = int(splitCoord[0]) + startChange
                        else:
                                splitCoord[1] = int(splitCoord[1]) - startChange
                if i == len(coords) - 1:
                        if orientation == '+':
                                splitCoord[1] = int(splitCoord[1]) - stopChange
                        else:
                                splitCoord[0] = int(splitCoord[0]) + stopChange
                coords[i] = splitCoord
        # Extend the CDS where possible
        result = cds_extension(coords, contigID, orientation, genomeRecords)
        if result == False:
                return False    # If we encountered a problem with our CDS, drop the model
        coords, cds = result
        # Validate the ORF to see that it is sufficiently similar to its origin
        result = validate_translated_cds(cds, cdsRecords, cdsID, alignPctCutoff)
        if result == False:
                return False
        # Return coordinates and protein sequence otherwise
        return coords, result

def find_longest_orf(seq, firstCodon):
        # Translate into ORFs and grab the longest bits inbetween stop codons
        longest = ['', '']
        for frame in range(3):
                record = Seq(seq)
                # Get nucleotide for this frame
                nucl = str(record)[frame:]
                nucl = Seq(nucl)
                # Translate to protein
                with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        frameProt = str(nucl.translate(table=1))
                # Find the longest ORF
                prots = frameProt.split('*')
                tmpLongest = ['', '']
                for i in range(len(prots)-1):                           # Ignore the last section; this has no stop codon
                        if len(prots[i]) > len(tmpLongest[0]):
                                tmpLongest = [prots[i], i]
                if tmpLongest == ['', '']:
                        continue                                        # This means the sequence starts with a stop codon, and the ORF itself lacks a stop codon
                # Convert this ORF back into its nucleotide sequence
                beforeLength = len(''.join(prots[:tmpLongest[1]]))*3 + (3*tmpLongest[1])        # 3*tmpLongest adds back in the length of any stop codons
                nuclOrf = nucl[beforeLength:beforeLength + len(tmpLongest[0])*3 + 3]            # +3 for the last stop codon
                nucl = str(nuclOrf)
                # Find the starting codon
                codonIndex = -1
                codons = re.findall('..?.?', nucl)                      # Pulls out a list of codons from the nucleotide
                for codon in codons:
                        if codon == firstCodon or codon == 'ATG':       # By adding an ATG check we ensure we don't stupidly skip the start codon looking for a silly codon start predicted by PASA
                                codonIndex = codons.index(codon)
                                break
                if codonIndex == -1:
                        continue
                # Update the start position of this nucl
                nucl = nucl[codonIndex*3:]
                # Translate to amino acid
                record = Seq(nucl)
                with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        frameOrf = str(record.translate(table=1))
                if len(frameOrf) > len(longest[0]):
                        longest = [frameOrf, nucl]
        return longest

def validate_translated_cds(cdsNucl, cdsRecords, cdsID, alignPctCutoff):
        # Arbitrary preset values
        cdsLenCutoff = 30       # 30 AA makes sure the ORF is longer than most spurious ORFs that arise by chance in a genome; if it's transcribed then it should be OK.
        # Get the original sequence's details
        origCDS = str(cdsRecords[cdsID].seq)
        origProt = longest_orf(cdsRecords[cdsID])
        origLen = len(origCDS)
        # Convert the new sequence to protein
        record = Seq(cdsNucl)
        with warnings.catch_warnings():
                warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                cdsProt = str(record.translate(table=1))
        assert cdsProt[-1] == '*' and cdsProt.count('*') == 1   # Make sure things are all good
        # If the ORF isn't long enough, it doesn't pass validation
        if len(cdsProt) < cdsLenCutoff:
                return False
        # If the ORFs are identical, they automatically pass validation
        if cdsProt == origProt:
                return cdsProt
        # Check that the two sequences are roughly the same - our extensions could have resulted in the longest ORF being within an extension
        try:
                newAlign, origAlign = ssw_simple(cdsProt, origProt)
        except:
                return False                                            # SSW dies sometimes. This seems to happen with repetitive sequences. While annoying, these errors can serve as a way of identifying bad sequences - bright side!
        alignPct = len(newAlign) / len(cdsProt)
        if alignPct < alignPctCutoff:                                   # Identity cut-off is probably not necessary, just align percent and arbitrary value to ensure the alignment covers most of the original sequence.
                return False                                            # Originally I tried 0.60, then 0.75. While these seem okay, I think sticking to 0.85 or 0.90 is ideal since we want to ensure that the ORF is mostly supported.
        # If we have the same start codon and a stop codon, check length for consistency
        lowerBound = origLen - (origLen * 0.1)
        upperBound = origLen + (origLen * 0.1)
        if lowerBound <= len(cdsNucl) <= upperBound:
                return cdsProt
        else:
                return False

def coord_excess_cut(coords, startChange, stopChange, orientation, allowMicroExon):
        # Cull exons that aren't coding and chop into coding exons
        startReduction = startChange
        stopReduction = stopChange
        startExonLen = 0
        stopExonLen = 0
        if allowMicroExon == False:
            microExonSize = -3      # This value is an arbitrary measure where, if a terminal exon is less than this size, we consider it 'fake' and delete it
        else:
            '''microExonSize is being changed as of 1/03/21. In the original vision of the program, it's probable that I saw potential for gmap to "fish"
            for a start or stop for a gene. This is probably especially possible if using CDS for the alignment. However, when using transcripts, this risk should
            be greatly reduced. At this date, this behaviour was causing problems in the program, so I needed to make this change. It's a bit risky, but I've set it
            as an optional parameter so it should be okay?'''
            microExonSize = 0
        for i in range(2):
                while True:
                        if i == 0:
                                exon = coords[0]
                        else:
                                exon = coords[-1]
                        # Extract details
                        rightCoord = int(exon[1])
                        leftCoord = int(exon[0])
                        exonLen = rightCoord - leftCoord + 1
                        # Update our change values
                        if i == 0:
                                startReduction -= exonLen       # This helps us to keep track of how much we need to reduce startChange
                                if startReduction > 0:          # when we begin chopping into the first exon - if the original first exon 
                                        del coords[0]           # is the one we chop, we end up with reduction value == 0
                                        startExonLen += exonLen
                                # Handle microexons at gene terminal
                                elif startReduction > microExonSize:
                                        del coords[0]
                                        startExonLen += exonLen
                                else:
                                        break
                        else:
                                stopReduction -= exonLen
                                if stopReduction > 0:
                                        del coords[-1]
                                        stopExonLen += exonLen  # We hold onto exon lengths so we can calculate how much we chop into the new start exon
                                # Handle microexons at gene terminal
                                elif stopReduction > microExonSize:
                                        del coords[-1]
                                        stopExonLen += exonLen
                                else:                           # by calculating stopChange - stopExonLen... if we didn't remove an exon, stopExonLen == 0
                                        break
        return coords, startExonLen, stopExonLen

def cds_extension(coords, contigID, orientation, genomeRecords):
        # Crawl up the genome sequence looking for a way to extend the ORF to 1) an ATG, or 2) the same codon
        stopCodonsPos = ['tag', 'taa', 'tga']
        stopCodonsNeg = ['cta', 'tta', 'tca']
        extensionLength = 90    # This is arbitrary; converts to 30 AA; can alter or consider making available as an argument
        cds = make_cds(coords, genomeRecords, contigID, orientation)
        currentStart = cds[0:3]
        currPos = None
        atgPos = None
        if orientation == '+':
                # Crawl back looking for the first stop codon - this is our boundary
                startCoord = int(coords[0][0])
                genomeSeq = genomeRecords[contigID][0:startCoord-1]     # startCoord is 1-based so we -1 to counter that
                for i in range(len(genomeSeq)-1, -1, -3):
                        codon = str(genomeSeq[i-2:i+1].seq)
                        if codon.lower() in stopCodonsPos:
                                break
                if str(genomeSeq.seq) == '':                            # Handles scenario where the gene model starts at the first base of the contig
                        i = 0
                else:
                        i = i - 2                                       # This walks our coordinate value back to the start of the codon (Atg) since our index currently corresponds to (atG)
                # Crawl back up from the stop position looking for the first current start or ATG
                for x in range(i+3, len(genomeSeq), 3):                 # +3 to look at the next, non-stop codon
                        codon = str(genomeSeq[x:x+3].seq)
                        if codon.lower() == currentStart.lower() and currPos == None:
                                currPos = x                     # Note that this X represents the distance in from the stop codon boundary
                        if codon.lower() == 'atg' and atgPos == None:
                                atgPos = x
                # Compare this position to the original based on length and codon type
                if atgPos != None:
                        if currentStart.lower() != 'atg':
                                accepted = [atgPos, True]               # We'll accept any length replacement if it's replacing a non-ATG with an ATG
                        elif (startCoord - atgPos) >= extensionLength:  # We'll only replace an ATG with another ATG if it increases length by a predefined amount
                                accepted = [atgPos, True]
                        else:
                                accepted = [0, False]
                else:
                        accepted = [0, False]                           # The false tag lets us know that 0 is not an extension but the lack thereof
                if currPos != None:
                        if (startCoord - currPos) >= accepted[0] + extensionLength:     # We always set accepted as a value above, whether it be an ATG start or 0
                                accepted = [currPos, True]
                # Recompute the accepted start position with reference to the entire contig
                acceptedPos = accepted[0] + 1                           # +1 to reconvert this to 1-based
                # Update this in our coords value if relevant
                if accepted[1] == True:
                        coords[0] = [acceptedPos, coords[0][1]]
        else:
                # Crawl up looking for the first stop codon - this is our boundary
                startCoord = int(coords[0][1])
                genomeSeq = genomeRecords[contigID][startCoord:]        # startCoord is 1-based; we want just after it, so accepting it as-is is correct
                for i in range(0, len(genomeSeq), 3):
                        codon = str(genomeSeq[i:i+3].seq)
                        if codon.lower() in stopCodonsNeg:
                                break
                if str(genomeSeq.seq) == '':                            # Handles scenario where the gene model starts at the last base of the contig
                        i = startCoord
                # Crawl back down from the stop position looking for the first current start or ATG
                currentStart = reverse_comp(currentStart)               # '-' orientation means we need to reverse complement our cds' start codon
                for x in range(i-1, -1, -3):
                        codon = str(genomeSeq[x-2:x+1].seq)
                        if codon.lower() == currentStart.lower() and currPos == None:
                                currPos = x                             # Note that this X represents the distance in from the stop codon boundary
                        if codon.lower() == 'cat' and atgPos == None:
                                atgPos = x
                # Compare this position to the original based on length and codon type
                if atgPos != None:
                        if currentStart.lower() != 'cat':
                                accepted = [atgPos, True]               # We'll accept any length replacement if it's replacing a non-ATG with an ATG
                        elif atgPos >= extensionLength:                 # We'll only replace an ATG with another ATG if it increases length by a predefined amount
                                accepted = [atgPos, True]
                        else:
                                accepted = [0, False]
                else:
                        accepted = [0, False]
                if currPos != None:
                        if currPos >= accepted[0] + extensionLength:    # We always set accepted as a value above, whether it be an ATG start or 0
                                accepted = [currPos, True]              # When we get here, we assume we accepted an ATG start, so we want to extend upon it even further to accept it as a valid extension
                # Recompute the accepted start position with reference to the entire contig
                acceptedPos = startCoord + accepted[0] + 1              # +1 to reconvert this to 1-based
                # Update this in our coords value
                if accepted[1] == True:
                        coords[0] = [coords[0][0], acceptedPos]
        # Determine if we need to do a stop codon crawl
        cds = make_cds(coords, genomeRecords, contigID, orientation)
        cdsRecord = Seq(cds)
        with warnings.catch_warnings():
                warnings.simplefilter('ignore')                         # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                cdsProt = str(cdsRecord.translate(table=1))
        if not cdsProt.count('*') < 2:                                  # Make sure the CDS is correct - it should be!
                return False                                            # This can happen when we strip off a micro exon prior to the ORF and, when extended, the first codon is immediately a stop.
        if cdsProt.count('*') == 1:
                assert cdsProt[-1] == '*'                               # If we have a stop codon, it should be at the end
                return coords, cds                                      # No need for a backwards crawl
        # Begin the stop codon crawl
        if orientation == '+':
                endCoord = int(coords[-1][1])
                # Trim off excess from the CDS to make sure we're in frame
                endCoord -= len(cds) % 3
                genomeSeq = genomeRecords[contigID][endCoord:]          # endCoord is 1-based; we want just after it, so accepting it as-is is correct
                for i in range(0, len(genomeSeq), 3):
                        codon = str(genomeSeq[i:i+3].seq)
                        if codon.lower() in stopCodonsPos:
                                break
                i = endCoord + i + 2 + 1                                # +2 to go to the end of the stop codon; +1 to make it 1-based
                coords[-1] = [coords[-1][0], i]
        else:
                endCoord = int(coords[-1][0])
                genomeSeq = genomeRecords[contigID][0:endCoord-1]       # endCoord is 1-based so we -1 to counter that
                for i in range(len(genomeSeq)-1, -1, -3):
                        codon = str(genomeSeq[i-2:i+1].seq)
                        if codon.lower() in stopCodonsNeg:
                                break
                i = i - 2 + 1                                           # -2 to go to the start of the codon; +1 to make it 1-based
                acceptedPos = accepted[0] + 1
                coords[-1] = [i, coords[-1][1]]
        # Make the final CDS and return
        cds = make_cds(coords, genomeRecords, contigID, orientation)
        return coords, cds

## Basic sequence operations
def make_cds(coords, genomeRecords, contigID, orientation):
        cds = []
        for i in range(len(coords)):
                splitCoord = coords[i]
                cdsBit = str(genomeRecords[contigID].seq)[int(splitCoord[0])-1:int(splitCoord[1])]
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                cds.append(cdsBit)
        cds = ''.join(cds)
        return cds

def longest_orf(record):
        longest = ''
        for frame in range(3):
                # Get nucleotide for this frame
                nucl = str(record.seq)[frame:]
                nucl = Seq(nucl)
                # Translate to protein
                with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        frameProt = str(nucl.translate(table=1))
                # Find the longest ORF
                prots = frameProt.split('*')
                tmpLongest = ''
                for i in range(len(prots)):
                        if len(prots[i]) > len(tmpLongest):
                                tmpLongest = prots[i]
                if len(tmpLongest) > len(longest):
                        longest = tmpLongest
        return longest

def ssw_simple(querySeq, targetSeq):
        # Perform SSW with parasail implementation
        profile = parasail.profile_create_sat(querySeq, parasail.blosum62)
        alignment = parasail.sw_trace_striped_profile_sat(profile, targetSeq, 10, 1)
        queryAlign = alignment.traceback.query
        targetAlign = alignment.traceback.ref
        return [queryAlign, targetAlign]

def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Complement characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

## Overlap collapsing
def compare_novels(inputDict, genomeRecords):
        def readd_removed_models(modelID, rejectedByDict, modelSets):
                if modelID in rejectedByDict:
                        for model in rejectedByDict[modelID]:
                                modelSets.append(model)
                        del rejectedByDict[modelID] # This just helps to lower memory usage slightly
                return rejectedByDict, modelSets
        def add_to_rejectedByDict(removedModel, rejectedByID, rejectedByDict):
                if rejectedByID not in rejectedByDict:
                        rejectedByDict[rejectedByID] = [removedModel]
                else:
                        rejectedByDict[rejectedByID].append(removedModel)
                return rejectedByDict
        # Setup
        rejectedByDict = {} # This value will allow us to re-add models that were rejected by models that themselves were later rejected
        # Find our contig IDs from the genome
        contigIDs = list(genomeRecords.keys())
        # Loop through our contigs and compare gene models
        acceptedModels = []
        for contig in contigIDs:
                # Gather all models on this contig
                contigModels = []
                for key, value in inputDict.items():
                        if value[1] == contig:
                                contigModels.append([key, value])
                # Convert to sets
                modelSets = []
                for model, value in contigModels:
                        coordSet = set()
                        for coord in value[0]:
                                splitCoord = coord
                                coordSet = coordSet.union(set(range(int(splitCoord[0]), int(splitCoord[1]) + 1))) # +1 to make the range up to and including the last digit
                        modelSets.append([model, coordSet, value[0], value[1], value[2]])
                # Compare sets to find overlaps
                loopEnd = False
                while True:
                        if loopEnd == True:
                                break
                        loopEnd = True
                        for i in range(len(modelSets)-1):
                                for x in range(i+1, len(modelSets)):
                                        set1 = modelSets[i][1]
                                        set2 = modelSets[x][1]
                                        # Check for overlap
                                        ovl = set1 & set2
                                        # If no overlap, continue
                                        if ovl == set():
                                                continue                                        
                                        # Handle overlaps
                                        # Filter 1: Microexons
                                        shortestExon1 = None
                                        shortestExon2 = None
                                        microLen = 30   # Arbitrary; exons longer than 30bp aren't considered microexons (this seems to be agreed upon in literature I briefly viewed)
                                        for coord in modelSets[i][2]:
                                                splitCoord = coord
                                                exonLen = int(splitCoord[1]) - int(splitCoord[0]) + 1
                                                if shortestExon1 == None:
                                                        shortestExon1 = exonLen
                                                elif exonLen < shortestExon1:
                                                        shortestExon1 = exonLen
                                        for coord in modelSets[x][2]:
                                                splitCoord = coord
                                                exonLen = int(splitCoord[1]) - int(splitCoord[0]) + 1
                                                if shortestExon2 == None:
                                                        shortestExon2 = exonLen
                                                elif exonLen < shortestExon2:
                                                        shortestExon2 = exonLen
                                        if shortestExon1 > shortestExon2:
                                                if shortestExon2 < microLen:
                                                        rejectedByDict = add_to_rejectedByDict(modelSets[x], modelSets[i][0], rejectedByDict)
                                                        rejectedByDict, modelSets = readd_removed_models(modelSets[x][0], rejectedByDict, modelSets)
                                                        del modelSets[x]
                                                        loopEnd = False
                                                        break
                                        elif shortestExon2 > shortestExon1:
                                                if shortestExon1 < microLen:
                                                        rejectedByDict = add_to_rejectedByDict(modelSets[i], modelSets[x][0], rejectedByDict)
                                                        rejectedByDict, modelSets = readd_removed_models(modelSets[i][0], rejectedByDict, modelSets)
                                                        del modelSets[i]
                                                        loopEnd = False
                                                        break
                                        # Filter 2: Length
                                        if len(set1) > len(set2):
                                                rejectedByDict = add_to_rejectedByDict(modelSets[x], modelSets[i][0], rejectedByDict)
                                                rejectedByDict, modelSets = readd_removed_models(modelSets[x][0], rejectedByDict, modelSets)
                                                del modelSets[x]
                                                loopEnd = False
                                                break
                                        elif len(set2) > len(set1):
                                                rejectedByDict = add_to_rejectedByDict(modelSets[i], modelSets[x][0], rejectedByDict)
                                                rejectedByDict, modelSets = readd_removed_models(modelSets[i][0], rejectedByDict, modelSets)
                                                del modelSets[i]
                                                loopEnd = False
                                                break
                                        # Filter 3: Splice rules
                                        spliceTypes1 = splice_sites(modelSets[i][2], genomeRecords, modelSets[i][3], modelSets[i][4])
                                        spliceTypes2 = splice_sites(modelSets[x][2], genomeRecords, modelSets[x][3], modelSets[x][4])
                                        canonPct1 = spliceTypes1[0] / sum(spliceTypes1)
                                        canonPct2 = spliceTypes2[0] / sum(spliceTypes2)
                                        noncanonPct1 = sum(spliceTypes1[1:3]) / sum(spliceTypes1)
                                        noncanonPct2 = sum(spliceTypes2[1:3]) / sum(spliceTypes2)
                                        if canonPct1 != canonPct2:
                                                if canonPct1 > canonPct2:
                                                        rejectedByDict = add_to_rejectedByDict(modelSets[x], modelSets[i][0], rejectedByDict)
                                                        rejectedByDict, modelSets = readd_removed_models(modelSets[x][0], rejectedByDict, modelSets)
                                                        del modelSets[x]
                                                        loopEnd = False
                                                        break
                                                else:
                                                        rejectedByDict = add_to_rejectedByDict(modelSets[i], modelSets[x][0], rejectedByDict)
                                                        rejectedByDict, modelSets = readd_removed_models(modelSets[i][0], rejectedByDict, modelSets)
                                                        del modelSets[i]
                                                        loopEnd = False
                                                        break
                                        elif noncanonPct1 != noncanonPct2:
                                                if noncanonPct1 > noncanonPct2:
                                                        rejectedByDict = add_to_rejectedByDict(modelSets[x], modelSets[i][0], rejectedByDict)
                                                        rejectedByDict, modelSets = readd_removed_models(modelSets[x][0], rejectedByDict, modelSets)
                                                        del modelSets[x]
                                                        loopEnd = False
                                                        break
                                                else:
                                                        rejectedByDict = add_to_rejectedByDict(modelSets[i], modelSets[x][0], rejectedByDict)
                                                        rejectedByDict, modelSets = readd_removed_models(modelSets[i][0], rejectedByDict, modelSets)
                                                        del modelSets[i]
                                                        loopEnd = False
                                                        break
                                        # If we pass all of these filters, we need to make a decision somehow
                                        ## Final decision 1: Shortest exon length
                                        if shortestExon1 > shortestExon2:
                                                rejectedByDict = add_to_rejectedByDict(modelSets[x], modelSets[i][0], rejectedByDict)
                                                rejectedByDict, modelSets = readd_removed_models(modelSets[x][0], rejectedByDict, modelSets)
                                                del modelSets[x]
                                                loopEnd = False
                                                break
                                        elif shortestExon2 > shortestExon1:
                                                rejectedByDict = add_to_rejectedByDict(modelSets[i], modelSets[x][0], rejectedByDict)
                                                rejectedByDict, modelSets = readd_removed_models(modelSets[i][0], rejectedByDict, modelSets)
                                                del modelSets[i]
                                                loopEnd = False
                                                break
                                        ## Final decision 2: Lower path number
                                        pathNum1 = int(modelSets[i][0].rsplit('.path', maxsplit=1)[1])
                                        pathNum2 = int(modelSets[x][0].rsplit('.path', maxsplit=1)[1])
                                        if pathNum1 < pathNum2:
                                                rejectedByDict = add_to_rejectedByDict(modelSets[x], modelSets[i][0], rejectedByDict)
                                                rejectedByDict, modelSets = readd_removed_models(modelSets[x][0], rejectedByDict, modelSets)
                                                del modelSets[x]
                                                loopEnd = False
                                                break
                                        elif pathNum2 < pathNum1:
                                                rejectedByDict = add_to_rejectedByDict(modelSets[i], modelSets[x][0], rejectedByDict)
                                                rejectedByDict, modelSets = readd_removed_models(modelSets[i][0], rejectedByDict, modelSets)
                                                del modelSets[i]
                                                loopEnd = False
                                                break
                                        ## Final decision 3: How!?!? Just kill x
                                        del modelSets[x]
                                        loopEnd = False
                                        break
                # Hold onto accepted models
                for entry in modelSets:
                        acceptedModels.append(entry[0])
        # Cull models from the dictionary that don't pass curation
        dictKeys = list(inputDict.keys())
        for key in dictKeys:
                if key not in acceptedModels:
                        del inputDict[key]
        # Return modified dictionary
        return inputDict
                                                
def splice_sites(coords, genomeRecords, contigID, orientation):
        # Extract bits to left and right of exon
        splices = []
        for i in range(len(coords)):
                splitCoord = coords[i]
                start = int(splitCoord[0])
                end = int(splitCoord[1])
                leftSplice = str(genomeRecords[contigID].seq)[start-1-2:start-1]        # -1 to make 1-based;-2 to go back 2 spots... -1 at end for going back to the base before the CDS
                rightSplice = str(genomeRecords[contigID].seq)[end-1+1:end-1+1+2]       # -1 to make 1-based;+1 to go to the base after CDS...-1 to make 1-based;+1 to go to the base after CDS;+2 to go to the end of the splice site
                # Hold onto splices with respect to orientation
                if i == 0:
                        if orientation == '+':
                                splices.append(rightSplice)
                        else:
                                splices.append(leftSplice)
                elif i == len(coords) - 1:
                        if orientation == '+':
                                splices.append(leftSplice)
                        else:
                                splices.append(rightSplice)
                else:
                        if orientation == '+':
                                splices.append(leftSplice)
                                splices.append(rightSplice)
                        else:
                                splices.append(rightSplice)
                                splices.append(leftSplice)
        # Reverse comp if necessary
        if orientation == '-':
                for i in range(len(splices)):
                        splices[i] = reverse_comp(splices[i])
        # Assess the number of canonical, somewhat common non-canonical, and very rare splices
        canonical = ['GT', 'AG']
        noncanonical = ['GC', 'AG']
        rare = ['AT', 'AC']
        spliceTypes = [0, 0, 0, 0]      # Refers to [canonical, noncanonical, rare, unknown]
        for i in range(0, len(splices), 2):
                left = splices[i].upper()
                right = splices[i+1].upper()
                if left == canonical[0] and right == canonical[1]:
                        spliceTypes[0] += 1
                elif left == noncanonical[0] and right == noncanonical[1]:
                        spliceTypes[1] += 1
                elif left == rare[0] and right == rare[1]:
                        spliceTypes[2] += 1
                else:
                        spliceTypes[3] += 1
        # Return value
        return spliceTypes

def merge_dictionaries(dictList):
        mergedDict = {}
        for dictionary in dictList:
                for key, value in dictionary.items():
                        assert key not in mergedDict    # This is the reason we have this as a function
                        mergedDict[key] = value         # The input files shouldn't have redundant names
        return mergedDict

## Curation checks
def remove_bad_splices(inputDict, genomeRecords):
        outputDict = {}
        for key, value in inputDict.items():
                # Arbitrary cutoffs
                minimumSize = 200       # An ORF that is >= 200 AA long is probably real assuming it is highly similar to its origin (which it is with a 0.85 cutoff elsewhere in the code)
                unknownCutoff = 0.33    # If more than 1/3 introns have unknown splices, there's a good chance this model isn't real
                # Extract details
                spliceTypes = splice_sites(value[0], genomeRecords, value[1], value[2])
                unknownPct = spliceTypes[3] / sum(spliceTypes)
                seqLen = 0
                for coord in value[0]:
                        splitCoord = coord
                        seqLen += int(splitCoord[1]) - int(splitCoord[0]) + 1
                # Drop any models that lack canonical splices entirely
                if spliceTypes[0] == 0: # This works really well, since it forces 2-exon genes to have a canonical splice and hence removes most obviously bad models
                        continue
                # Drop any short models with unknown splices
                if spliceTypes[3] != 0 and seqLen/3 < minimumSize:
                        continue
                # Drop any models with > cutoff percent of unknown splices
                if round(unknownPct, 2) > unknownCutoff:
                        continue
                # Hold onto anything that passes cutoff
                outputDict[key] = value
        return outputDict

def remove_weird_models(inputDict, ovlDict):
        outputDict = {}
        # Set up arbitrary cutoffs
        intronLenMin = 50       # Any model which only contains super short introns is probably not a good model
        intronLenMax = 10000    # Models which overlap existing ones with large introns are often GMAP fishing for a domain alignment - not a good model
        for key, value in inputDict.items():
                # Check 1: Only micro introns (probably a sign that GMAP has done some weird stuff to get an in-frame match)
                intronCoords, intronLens = intron_detail_extract(value[0], value[2])
                if max(intronLens) < intronLenMin:
                        continue
                # Check 2: Huge introns when gene overlaps existing (probably a sign that we're incorrectly aligning against a domain region)
                if ovlDict[key] != 0 and max(intronLens) > intronLenMax:
                        continue
                # If we pass the above checks, hold onto this result
                outputDict[key] = value
        return outputDict

def intron_detail_extract(coords, orientation):
        intronCoords = []
        intronLens = []
        for i in range(0, len(coords)-1):
                # Extract details
                splitCoord1 = coords[i]
                start1 = int(splitCoord1[0])
                stop1 = int(splitCoord1[1])
                splitCoord2 = coords[i+1]
                start2 = int(splitCoord2[0])
                stop2 = int(splitCoord2[1])
                # Identify intron gap depending on orientation
                if orientation == '+':
                        intronCoord = str(stop1+1) + '-' + str(start2-1)
                        intronLen = (start2-1) - (stop1+1)
                else:
                        intronCoord = str(stop2+1) + '-' + str(start1-1)
                        intronLen = (start1-1) - (stop2+1)
                intronCoords.append(intronCoord)
                intronLens.append(intronLen)
        return intronCoords, intronLens

## NCLS RELATED
def gff3_parse_ncls_mrna(gff3File):                             # This function will make a NCLS object which can be used to find gene model overlaps; note that this is the whole gene's range, not separate exon ranges
        gff3Loc = {}
        starts = []
        ends = []
        ids = []
        ongoingCount = 0
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        line = line.replace('\r', '')   # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                        # Skip unneccessary lines
                        if line.startswith('#') or line == '\n':
                                continue
                        sl = line.split('\t')
                        if len(sl) < 8:                 # If the length is shorter than this, it's not a gene detail line
                                continue
                        # Skip non-mRNA lines
                        if sl[2] != 'mRNA':
                                continue
                        # Get details from line including start, stop, and orientation
                        contigID = sl[0]
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        orient = sl[6]
                        details = sl[8].split(';')
                        detailDict = {}
                        for i in range(len(details)):
                                if details[i] == '':
                                        continue
                                splitDetail = details[i].split('=')
                                detailDict[splitDetail[0]] = splitDetail[1]
                        # Add to our NCLS
                        starts.append(contigStart)
                        ends.append(contigStop+1)       # NCLS indexes 0-based like a range (up to but not including end), so +1 to make this more logically compliant with gff3 1-based system.
                        ids.append(ongoingCount)
                        gff3Loc[ongoingCount] = [contigStart, contigStop, orient, detailDict['ID'], contigID]
                        ongoingCount += 1
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        return ncls, gff3Loc

def ncls_finder(ncls, locDict, start, stop, featureID, featureIndex):
        overlaps = ncls.find_overlap(start, stop+1)             # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        dictEntries = []
        for result in overlaps:
                dictEntries.append(locDict[result[2]])
        dictEntries = copy.deepcopy(dictEntries)                # Any time we're deleting things from a section of a dictionary we need to build a deepcopy to keep the original dictionary intact.
        # Narrow down our dictEntries to hits to the same feature
        dictEntries = ncls_feature_narrowing(dictEntries, featureID, featureIndex)
        # Return list
        return dictEntries

def ncls_feature_narrowing(nclsEntries, featureID, featureIndex):
        for k in range(len(nclsEntries)-1, -1, -1):
                if nclsEntries[k][featureIndex] != featureID:
                        del nclsEntries[k]
        return nclsEntries

def output_func(inputDict, outFileName):
        with open(outFileName, 'w') as fileOut:
                for key, value in inputDict.items():
                        # Format base triticeae details
                        pathID = value[1] + '.' + key
                        name = 'gmap_gene_find_' + key
                        mrnaID = pathID.replace('.path', '.mrna')       # Could theoretically be a problem if the gene triticeae contains .path in its actual triticeae, but this isn't the case with my transposon and shouldn't be with others
                        # Extract details
                        firstCoord = value[0][0]
                        firstInts = [int(firstCoord[0]), int(firstCoord[1])]
                        lastCoord = value[0][-1]
                        lastInts = [int(lastCoord[0]), int(lastCoord[1])]
                        protein = value[3]
                        # Determine gene start and end coordinates with respect to orientation
                        if value[2] == '+':
                                start = min(firstInts)
                                end = max(lastInts)
                        else:
                                end = max(firstInts)
                                start = min(lastInts)
                        # Format start comment
                        startComment = '# GMAP_GENE_FIND: ' + mrnaID + ' automatic model build'
                        fileOut.write(startComment + '\n')
                        # Format gene line
                        typeCol = 'gmap_gene_find'
                        geneLine = '\t'.join([value[1], typeCol, 'gene', str(start), str(end), '.', value[2], '.', 'ID=' + pathID +';Name=' + name])
                        fileOut.write(geneLine + '\n')
                        # Format mRNA line
                        mrnaLine = '\t'.join([value[1], typeCol, 'mRNA', str(start), str(end), '.', value[2], '.', 'ID=' + mrnaID +';Name=' + name + ';Parent=' + pathID])
                        fileOut.write(mrnaLine + '\n')
                        # Derive phasing information from coordinates
                        totalLen = 0
                        phasing = []
                        for i in range(len(value[0])):
                                coord = value[0][i]
                                segmentLen = coord[1] - coord[0] + 1
                                totalLen += segmentLen
                                if i == 0:
                                        phasing.append('0')     # GGF always returns ORFs which are 0-phased on the first CDS bit
                                else:
                                        prevLen = totalLen - segmentLen
                                        leftover = prevLen % 3
                                        if leftover == 0:
                                                phase = 0
                                        else:
                                                phase = 3 - leftover
                                        phasing.append(str(phase))
                        # Iterate through coordinates and write exon/CDS lines
                        ongoingCount = 1
                        for i in range(len(value[0])):
                                start, end = value[0][i]
                                phase = phasing[i]
                                # Format exon line
                                exonLine = '\t'.join([value[1], typeCol, 'exon', str(start), str(end), '.', value[2], '.', 'ID=' + mrnaID + '.exon' + str(ongoingCount) + ';Name=' + name + ';Parent=' + mrnaID])
                                fileOut.write(exonLine + '\n')
                                # Format CDS line
                                cdsLine = '\t'.join([value[1], typeCol, 'CDS', str(start), str(end), '.', value[2], phase, 'ID=' + mrnaID + '.cds' + str(ongoingCount) + ';Name=' + name + ';Parent=' + mrnaID])
                                fileOut.write(cdsLine + '\n')
                                ongoingCount += 1
                        # Format end comment
                        endComment = '#PROT ' + mrnaID + ' ' + pathID + '\t' + protein
                        fileOut.write(endComment + '\n')

# Main block
def main():
        def validate_args(args):
                # Ensure no None arguments exist
                for key, value in vars(args).items():
                        if value == None:
                                print(key + ' argument was not specified. Fix this and try again.')
                                quit()
                # Validate input file locations
                for gmapFile in args.gmapFiles:
                        if not os.path.isfile(gmapFile):
                                print('I am unable to locate the input GMAP gff3 file (' + gmapFile + ')')
                                print('Make sure you\'ve typed the file triticeae or location correctly and try again.')
                                quit()
                for cdsFile in args.cdsFiles:
                        if not os.path.isfile(cdsFile):
                                print('I am unable to locate the input CDS FASTA file (' + cdsFile + ')')
                                print('Make sure you\'ve typed the file triticeae or location correctly and try again.')
                                quit()
                if len(args.gmapFiles) != len(args.cdsFiles):
                        print('There is a different number of arguments provided for gmapFiles and cdsFiles.')
                        print('Each gmap GFF3 file should be matched with its respective CDS FASTA file. Try again.')
                        quit()
                if not os.path.isfile(args.genomeFile):
                        print('I am unable to locate the input genome FASTA file (' + args.genomeFile + ')')
                        print('Make sure you\'ve typed the file triticeae or location correctly and try again.')
                        quit()
                if not os.path.isfile(args.annotationFile):
                        print('I am unable to locate the input genome annotation gff3 file (' + args.annotationFile + ')')
                        print('Make sure you\'ve typed the file triticeae or location correctly and try again.')
                        quit()
                # Validate numerical arguments
                if not 0 <= args.coverageCutoff <= 100.0:
                        print('Coverage cut-off must be any number >= 0.0 and <= 100.0. Try again.')
                        quit()
                if not 0 <= args.identityCutoff <= 100.0:
                        print('Identity cut-off must be any number >= 0.0 and <= 100.0. Try again.')
                        quit()
                if not 0 <= args.alignPctCutoff <= 100.0:
                        print('Identity cut-off must be any number >= 0.0 and <= 100.0. Try again.')
                        quit()
                args.alignPctCutoff = args.alignPctCutoff / 100 # I think it's more intuitive on the commandline to deal with percentages 0-100 rather than ratios 0-1
                # Handle file overwrites
                if os.path.isfile(args.outputFileName):
                        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                        quit()
                return args
        
        #### USER INPUT SECTION
        usage = """%(prog)s reads in 1 or more input GMAP GFF3 files built using -f 2 formatting and will compare
        these against the current annotation, returning a GFF3 file of curated gene models that pass a variety of checks.
        These include intron splicing checks, overlap checks between new models and existing, and checks that the new model
        is similar to the old one. The result should need minimal manual curation, and likely contains many models that are
        either inside introns or are extra copies of other models that weren't assembled due to PASA's default GMAP settings
        resulting in only 1 path.
        """
        
        # Reqs
        p = argparse.ArgumentParser(description=usage)
        p.add_argument("-gm", "-gmapFiles", dest="gmapFiles", nargs="+",
                           help="Input GMAP gene gff3 (-f 2) file triticeae(s).")
        p.add_argument("-cd", "-cdsFiles", dest="cdsFiles", nargs="+",
                           help="Input CDS fasta file triticeae(s) (these file[s] were used for GMAP alignment).")
        p.add_argument("-ge", "-genomeFile", dest="genomeFile",
                           help="Input genome FASTA file triticeae.")
        p.add_argument("-an", "-annotationFile", dest="annotationFile",
                           help="Input current genome annotation gff3 file triticeae.")
        p.add_argument("-co", "-coverage", dest="coverageCutoff", type=float,
                           help="Coverage cut-off (must have coverage >= provided value; accepted range 0.0->100.0; default == 70.0).", default=70.0)
        p.add_argument("-id", "-identity", dest="identityCutoff", type=float,
                           help="Identity cut-off (must have identity >= provided value; accepted range 0.0->100.0; default == 90.0).", default=90.0)
        p.add_argument("-al", "-alignPctCutoff", dest="alignPctCutoff", type=float,
                           help="Alignment percent cut-off (new sequence must align against original >= provided value; accepted range 0.0->100.0; default == 90.0).", default=90.0)
        p.add_argument("-o", "-outputFile", dest="outputFileName",
                           help="Output file triticeae.")
        p.add_argument("-mexon", dest="allowMicroExon", action='store_true',
                            help="Optional parameter to set if you are using transcripts (NOT CDS) for the alignment, to change how terminal micro-exons are handled.", default=False)
        
        args = p.parse_args()
        args = validate_args(args)
        
        # Load in genome fasta file as dict
        genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'r'), 'fasta'))
        
        # Parse main annotation GFF3 as NCLS and model
        annotation_ncls, annotation_ncls_loc = gff3_parse_ncls_mrna(args.annotationFile)
        annotation_gff3 = Gff3(args.annotationFile)
        
        # Main loop: find good multipath gene models
        gmapDicts = [] # Hold onto dictionaries from each iteration of gmap files
        OVERLAP_CUTOFF = 0.35 # Arbitrary value; only sharing about 1/3 of its length with known genes seems appropriate
        ovlDict = {} # Hold onto ovl percentages for bad model checking
        for gmap_file, cds_file in zip(args.gmapFiles, args.cdsFiles):
                # Load in CDS fasta file as dict
                cdsRecords = SeqIO.to_dict(SeqIO.parse(open(cds_file, 'r'), 'fasta'))
                # Parse GMAP GFF3 for gene models
                gmap_gff3 = Gff3(gmap_file)
                # Detect well-suported models
                coordDict = {}
                for key in gmap_gff3.mrna_values:
                        value = gmap_gff3.index_dict[key][key]
                        # Skip if the current path is a 1-exon gene; some 1-exon genes will be real, but there's a much higher chance of it being a pseudogene or transposon-related gene
                        if len(value['exon']['coords']) == 1:
                                continue
                        # Skip if the transcript lacks a stop codon and thus is likely to be a fragment
                        transcriptID = key.rsplit('.', maxsplit=1)[0] # This removes the '.mrna#' suffix which won't be present in the original FASTA file
                        lastCodon = str(cdsRecords[transcriptID].seq)[-3:]
                        if lastCodon.lower() not in ['tag', 'taa', 'tga']:
                                continue
                        # Cut-off checks
                        decision = check_model(value['attributes'], args.coverageCutoff, args.identityCutoff)
                        if decision == False:
                                continue
                        result = cds_build(value['exon']['coords'], value['contig_id'], value['orientation'], cdsRecords, genomeRecords, transcriptID, args.alignPctCutoff, args.allowMicroExon) # Split off the '.mrna#' suffix ## CHECK THIS OUT
                        if result == False:
                                continue
                        else:
                                coords, protein = result
                                coordDict[key] = [coords, value['contig_id'], value['orientation'], protein]
                # Remove overlaps of existing genes
                outputValues = {}
                for key, value in coordDict.items():
                        # Find overlaps for each exon
                        dictEntries = []
                        for coord in value[0]:
                                start, stop = coord
                                # Find overlaps
                                dictEntries += ncls_finder(annotation_ncls, annotation_ncls_loc, start, stop, value[1], 4) # 4 refers to the index position in the annotation_ncls_loc dictionary for the contigID
                        # Convert coordinates to set values for overlap calculation
                        valueSet = set()
                        for coord in value[0]:
                                start, stop = coord
                                valueSet = valueSet.union(set(range(start, stop+1)))
                        # Compare overlaps to see if this gene overlaps existing genes
                        checked = []
                        overlapped = set()
                        for entry in dictEntries:
                                # Handle redundancy
                                if entry in checked:
                                        continue
                                checked.append(entry)
                                # Retrieve this mRNA value from the gff3Dict object
                                mrnaHit = annotation_gff3.index_dict[entry[3]][entry[3]]
                                # Calculate the overlap of the current value against this mRNA model
                                mrnaSet = set()
                                for coord in mrnaHit['CDS']['coords']:
                                        start, stop = coord
                                        mrnaSet = mrnaSet.union(set(range(start, stop+1)))
                                        # Store how much overlap is present
                                        overlapped = overlapped.union(valueSet & mrnaSet)
                        # Drop any models which overlap existing genes enough to suggest that it's either 1) fragmented, or 2) an isoform
                        ovlPct = (len(valueSet) - (len(valueSet) - len(overlapped))) / len(valueSet) # Originally I was going to drop anything with overlap, but manual inspection shows there's justification for allowing some overlap (but not much)
                        if ovlPct > OVERLAP_CUTOFF:
                                continue
                        # Hold onto things which pass this check
                        else:
                                outputValues[gmap_gff3.index_dict[key]['attributes']['ID']] = value
                                ovlDict[gmap_gff3.index_dict[key]['attributes']['ID']] = ovlPct
                # Hold onto the result
                gmapDicts.append(outputValues)
        
        # Merge dictionaries together
        mergedDict = merge_dictionaries(gmapDicts)
        
        # Collapse overlapping paths by selecting the 'best' according to canonical splicing and other rules
        mergedDict = compare_novels(mergedDict, genomeRecords)
        
        # Cull models with poor splicing
        mergedDict = remove_bad_splices(mergedDict, genomeRecords)
        
        # Cull weird looking models
        mergedDict = remove_weird_models(mergedDict, ovlDict)
        
        # Output to file
        output_func(mergedDict, args.outputFileName)
        
        # Done!
        print('Program completed successfully!')
        print(str(len(mergedDict)) + ' models were discovered.')

if __name__ == '__main__':
        main()
