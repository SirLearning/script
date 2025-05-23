#! python3
# gff3_handling.py
# This file is intended to serve as the central storage location for GFF3-related
# classes and functions rather than having multiple versions spread across multiple files.

# Modules imported for these functions
import re, copy
import pandas as pd
from ncls import NCLS
from Bio import SeqIO

# Gff3 class definition
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
                                details = sl[8].strip("\"").split(';')
                                detail_dict = {}
                                for i in range(len(details)):
                                        if details[i] == '':
                                                continue
                                        split_details = details[i].split('=', maxsplit=1)
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
        ## Lines
        def add_lines(self):
                # Setup
                main_types = list(self.id_values['main'].keys())
                KNOWN_HEAD_COMMENTS = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND', '# GEMOMA ANNOTATION', '# APOLLO ANNOTATION') # These are the comment lines we'll handle within this code; anything not like this is ignored
                KNOWN_FOOT_COMMENTS = ('#PROT')
                assert self.file_location != None
                # Main loop
                with open(self.file_location, 'r') as file_in:
                        for line in file_in:
                                line = line.replace('\r', '') # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                                # Skip filler lines
                                if line == '\n' or set(line.rstrip('\n')) == {'#'} or set(line.rstrip('\n')) == {'#', '\t'}: # If this is true, it's a blank line or a comment line with no information in it
                                        continue
                                sl = line.rstrip('\n').split('\t')
                                # Handle known header comment lines
                                if line.startswith(KNOWN_HEAD_COMMENTS):
                                        # Extract gene ID
                                        mrna_ID = line.split(': ')[1].split(' ')[0].rstrip(',') # According to known header comments, the mRNA ID will be found inbetween ': ' and ' ' with a possible comma at the end which we can strip off
                                        gene_ID = self.index_dict[mrna_ID]['attributes']['ID'] # mrna_ID indexes back to the main gene dict object, and from here we can get the geneID from its attributes field
                                        # Add to lines dict
                                        if 'lines' not in self.index_dict[gene_ID]:
                                                self.index_dict[gene_ID]['lines'] = {0: [line], 1: [], 2: []}
                                        else:
                                                self.index_dict[gene_ID]['lines'][0].append(line)
                                # Handle known footer comment lines
                                elif line.startswith(KNOWN_FOOT_COMMENTS):
                                        # Extract gene ID
                                        gene_ID = line.split()[2] # According to known footer comments, the gene ID will be the third 1-based value (e.g., ['#PROT', 'evm.model.utg0.34', 'evm.TU.utg0.34', 'MATEDAP....'])
                                        # Add to lines dict
                                        if 'lines' not in self.index_dict[gene_ID]:
                                                self.index_dict[gene_ID]['lines'] = {0: [], 1: [], 2: [line]}
                                        else:
                                                self.index_dict[gene_ID]['lines'][2].append(line)
                                # Handle feature detail lines
                                elif not line.startswith('#'):
                                        # Extract gene ID
                                        attributes = sl[8].split(';')
                                        if sl[2] in main_types:
                                                for attr in attributes:
                                                        if attr.startswith('ID='): # For main-type lines, the ID= is our gene/feature ID
                                                                gene_ID = attr[3:].strip('\n') # This trims off the ID= bit and any new lines
                                        else:
                                                gene_or_mrna_ID = None
                                                for attr in attributes:
                                                        if attr.startswith('Parent='): # For every other type of line, the Parent= field should tell us the geneID or mrnaID
                                                                gene_or_mrna_ID = attr[7:].strip('\n') # This trims off the Parent= bit and any new lines
                                                if gene_or_mrna_ID == None: # This will handle biological_region (ctrl+f for this reference in gff3_index()) and other values which lack ID= and Parent= fields; we don't index these since they are (currently) of no interest
                                                        continue
                                                if gene_or_mrna_ID in self.index_dict:
                                                        gene_ID = self.index_dict[gene_or_mrna_ID]['attributes']['ID'] # This lets us handle the ambiguity of our geneORmrnaID and make sure we're looking at the geneID
                                                elif ',' in gene_or_mrna_ID: # This is for specific scenarios like in TAIR9 where a feature has multiple parents
                                                        gene_ID = gene_or_mrna_ID.split(',')
                                        # Add to lines dict
                                        if type(gene_ID) != list:
                                                if 'lines' not in self.index_dict[gene_ID]:
                                                        self.index_dict[gene_ID]['lines'] = {0: [], 1: [line], 2: []}
                                                else:
                                                        self.index_dict[gene_ID]['lines'][1].append(line)
                                        else:                                                   # This section relates to the immediately above comment when handling multiple parent features
                                                for parent in gene_ID:                          # In this case, gene_ID is a list of parents
                                                        parent_text = line.split('Parent=')[1]
                                                        parent_text = parent_text.split(';')[0] # This will extract just the bit of the comment from Parent= to any potential ; after
                                                        new_line = line.replace(parent_text, parent)
                                                        if 'lines' not in self.index_dict[parent]:
                                                                self.index_dict[parent]['lines'] = {0: [], 1: [new_line], 2: []} # We do all of this so we can separate multi-parent features into individual bits
                                                        else:
                                                                self.index_dict[parent]['lines'][1].append(new_line) # I think that multi-parent features shouldn't exist in GFF3 since, if they do, it's probably redundant or compressing information too much
                                # All other lines are ignored
        
        def add_comments(self): # This function is just add_lines but with the gene lines section gutted
                # Setup
                KNOWN_HEAD_COMMENTS = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND', '# GEMOMA ANNOTATION', '# APOLLO ANNOTATION') # These are the comment lines we'll handle within this code; anything not like this is ignored
                KNOWN_FOOT_COMMENTS = ('#PROT')
                assert self.file_location != None
                # Main loop
                with open(self.file_location, 'r') as file_in:
                        for line in file_in:
                                line = line.replace('\r', '') # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                                # Skip filler lines
                                if line == '\n' or set(line.rstrip('\n')) == {'#'} or set(line.rstrip('\n')) == {'#', '\t'}: # If this is true, it's a blank line or a comment line with no information in it
                                        continue
                                # Handle known header comment lines
                                if line.startswith(KNOWN_HEAD_COMMENTS):
                                        # Extract gene ID
                                        mrna_ID = line.split(': ')[1].split(' ')[0].rstrip(',') # According to known header comments, the mRNA ID will be found inbetween ': ' and ' ' with a possible comma at the end which we can strip off
                                        gene_ID = self.index_dict[mrna_ID]['attributes']['ID'] # mrna_ID indexes back to the main gene dict object, and from here we can get the geneID from its attributes field
                                        # Add to lines dict
                                        if 'lines' not in self.index_dict[gene_ID]:
                                                self.index_dict[gene_ID]['lines'] = {0: [line], 1: [], 2: []}
                                        else:
                                                self.index_dict[gene_ID]['lines'][0].append(line)
                                # Handle known footer comment lines
                                elif line.startswith(KNOWN_FOOT_COMMENTS):
                                        # Extract gene ID
                                        gene_ID = line.split()[2] # According to known footer comments, the gene ID will be the third 1-based value (e.g., ['#PROT', 'evm.model.utg0.34', 'evm.TU.utg0.34', 'MATEDAP....'])
                                        # Add to lines dict
                                        if 'lines' not in self.index_dict[gene_ID]:
                                                self.index_dict[gene_ID]['lines'] = {0: [], 1: [], 2: [line]}
                                        else:
                                                self.index_dict[gene_ID]['lines'][2].append(line)
                                # Handle all other lines
                                else:
                                        pass

        def denovo_lines(self): # see mtdna_gff2_order.py for an example of this function
                COMMENT_ORDER = ['ID', 'Parent', 'Name'] # Hard code some of the ordering for GFF3 comments; any other type might end up being randomised a bit
                gff3_lines = {}
                for key in self.primary_values:
                        # Produce lines-type structure
                        gff3_lines[key] = {0: [], 1: [], 2: []}
                        # Produce our main feature line
                        gff3_entry = self.index_dict[key]
                        feature_line = [gff3_entry['contig_id'], gff3_entry['source'], gff3_entry['feature_type'], str(gff3_entry['coords'][0]), 
                                       str(gff3_entry['coords'][1]), gff3_entry['score'], gff3_entry['orientation'], gff3_entry['frame']]
                        # Format comment with some degree of ordering
                        comments = ''
                        for k in COMMENT_ORDER:
                                if k in gff3_entry['attributes']:
                                        if comments == '':
                                                comments += k + '=' + gff3_entry['attributes'][k]
                                        else:
                                                comments += ';' + k + '=' + gff3_entry['attributes'][k]
                        for k, v in gff3_entry['attributes'].items():
                                if k in COMMENT_ORDER:
                                        continue
                                if comments == '':
                                        comments += k + '=' + gff3_entry['attributes'][k]
                                else:
                                        comments += ';' + k + '=' + gff3_entry['attributes'][k]
                        feature_line.append(comments)
                        gff3_lines[key][1].append(feature_line)
                        # Produce subfeature line
                        for subfeature in gff3_entry['feature_list']:
                                subgff3_entry = gff3_entry[subfeature]
                                feature_line = [subgff3_entry['contig_id'], subgff3_entry['source'], subgff3_entry['feature_type'], str(subgff3_entry['coords'][0]), 
                                               str(subgff3_entry['coords'][1]), subgff3_entry['score'], subgff3_entry['orientation'], subgff3_entry['frame']]
                                # Format comment with some degree of ordering
                                comments = ''
                                for k in COMMENT_ORDER:
                                        if k in subgff3_entry['attributes']:
                                                if comments == '':
                                                        comments += k + '=' + subgff3_entry['attributes'][k]
                                                else:
                                                        comments += ';' + k + '=' + subgff3_entry['attributes'][k]
                                for k, v in subgff3_entry['attributes'].items():
                                        if k in COMMENT_ORDER:
                                                continue
                                        if comments == '':
                                                comments += k + '=' + subgff3_entry['attributes'][k]
                                        else:
                                                comments += ';' + k + '=' + subgff3_entry['attributes'][k]
                                feature_line.append(comments)
                                gff3_lines[key][1].append(feature_line)
                                # Produce secondary subfeature line(s)
                                for subfeature_type in subgff3_entry['feature_list']:
                                        secondarysubgff3_entry = subgff3_entry[subfeature_type]
                                        for i in range(len(secondarysubgff3_entry['coords'])):
                                                feature_line = [subgff3_entry['contig_id'], subgff3_entry['source'], subfeature_type, str(secondarysubgff3_entry['coords'][i][0]),
                                                               str(secondarysubgff3_entry['coords'][i][1]), secondarysubgff3_entry['score'][i], subgff3_entry['orientation'], secondarysubgff3_entry['frame'][i]]
                                                # Format comment with some degree of ordering
                                                comments = ''
                                                for k in COMMENT_ORDER:
                                                        if k in secondarysubgff3_entry['attributes'][i]:
                                                                if comments == '':
                                                                        comments += k + '=' + secondarysubgff3_entry['attributes'][i][k]
                                                                else:
                                                                        comments += ';' + k + '=' + secondarysubgff3_entry['attributes'][i][k]
                                                for k, v in secondarysubgff3_entry['attributes'][i].items():
                                                        if k in COMMENT_ORDER:
                                                                continue
                                                        if comments == '':
                                                                comments += k + '=' + secondarysubgff3_entry['attributes'][i][k]
                                                        else:
                                                                comments += ';' + k + '=' + secondarysubgff3_entry['attributes'][i][k]
                                                feature_line.append(comments)
                                                gff3_lines[key][1].append(feature_line)
                        # Reformat lines into strings
                        for x in range(len(gff3_lines[key][1])):
                                gff3_lines[key][1][x] = '\t'.join(gff3_lines[key][1][x]) + '\n'
                        # Associate lines to gff3Index
                        self.index_dict[key]['lines'] = gff3_lines[key]
        
        # Extraction of details
        def pasaprots_extract(self):
                # Setup
                self.pasa_prots = {}
                # Main loop
                for key in self.gene_values:
                        if 'lines' not in self.index_dict[key]:
                                continue
                        foot_comments = self.index_dict[key]['lines'][2]
                        # Parse each foot comment to extract the protein sequence
                        for comment in foot_comments:
                                split_comment = comment.rstrip('\r\n').split('\t')
                                # Extract the mRNA ID
                                mrnaID = split_comment[0].split(' ')[1] # Format for PASA comments after ' ' split should be ['#PROT', mrnaID, geneID]
                                # Extract the sequence
                                sequence = split_comment[1]
                                # Add into output dict
                                assert mrnaID not in self.pasa_prots # If this assertion fails, GFF3 comment format is flawed - there is a duplicate mRNA ID
                                self.pasa_prots[mrnaID] = sequence

## GFF3 - gene ID manipulation
def gff3_idlist_compare(gff3Dict, idList):
        # Set up
        outList = []
        # Main loop
        for key, value in gff3Dict.items():
                # Extract parent details from comment-containing value
                mrnaParent = None
                comment = value[-1][4].split(';')               # Our last value will always contain the GFF3 comment; we only need this once to get the parent ID
                for section in comment:
                        if section.startswith('Parent='):
                                mrnaParent = section[7:]
                assert mrnaParent != None
                # Check if the user specified a gene ID for removal/retrieval
                found = False
                if mrnaParent in idList:
                        found = True
                # Check if the user specified a mRNA ID for removal/retrieval
                for mrna in value:
                        mrnaID = mrna[0]
                        if mrnaID in idList:
                                found = True
                # If we found this ID in some capacity (as a gene or mRNA ID) within our idList, put the parent and all mRNAs in this list
                if found == True:
                        outList.append(mrnaParent)
                        for mrna in value:
                                outList.append(mrna[0])
        # Remove redundancy that may have crept in
        outList = list(set(outList))
        return outList

## GFF3 - filtering and curation
def gff3_index_cutoff_candidates(gff3Object, attributeList, cutoffList, directionList, behaviour):       # In this function, we can provide three paired lists of attributes, cutoffs, and the direction of comparison so it is generalisable to different uses
        # Setup
        outputList = []
        # Ensure that behaviour is sensible                                                             # directionList should be a list like ['<', '>', '<=', '>=', '=='] to indicate how the attribute's value should be held up against the cutoff
        if behaviour.lower() not in ['main', 'feature']:
                print('gff3_index_candidates: behaviour must be specified as "main" or "feature"; fix the code for this section.')
                quit()
        # Convert to lists if provided as strings/ints/etc
        if type(attributeList) != list:
                attributeList = [attributeList]
        if type(cutoffList) != list:
                cutoffList = [cutoffList]
        # Ensure that attributes and cutoffs are sensible
        if len(attributeList) != len(cutoffList) or len(attributeList) != len(directionList):
                print('gff3_index_candidates: attributesList, cutoffList, and/or directionList lengths are not equivalent; these lists should be paired. Fix the code for this section.')
                quit()
        # Ensure that cutoffs are capable of float conversion
        for cutoff in cutoffList:
                try:
                        float(cutoff)
                except:
                        print('gff3_index_candidates: ' + str(cutoff) + ' is provided as a cutoff, but is not capable of conversion to float. This should not happen; fix the code for this section.')
                        quit()
        # Loop through all indexed features and return viable candidates which 1: have all the attributes mentioned in our list and 2: their values pass our cutoff
        for featureType in gff3Object.id_values['feature']:
                for feature in gff3Object.id_values['feature'][featureType]:
                        # Set up this gene object's attribute:cutoff check values
                        cutoffCheck = [0]*len(attributeList)  # If we don't find any attributes or none of them pass cutoff, the sum of this will be 0; if we find them all and they all pass cutoff, the sum will == len(attributesList)
                        # Check gene object for attributes
                        geneObj = gff3Object.index_dict[feature]
                        for key, value in geneObj['attributes'].items():
                                if key in attributeList:
                                        # Ensure that this attribute works
                                        try:
                                                float(value)
                                        except:
                                                print('gff3_index_candidates: ' + str(value) + ' was found paired to attribute "' + key + '" but is not capable of float conversion. This attribute is thus not suitable for cutoff criterion checking; fix your GFF3 or do not specify this attribute.')
                                                quit()
                                        # Retrieve the index of this attribute
                                        attribIndex = attributeList.index(key)
                                        # Check against cutoff
                                        if eval(str(value) + directionList[attribIndex] + str(cutoffList[attribIndex])):
                                                cutoffCheck[attribIndex] = 1
                        # Check features within gene object if relevant
                        if sum(cutoffCheck) != len(attributeList):
                                for featureID in geneObj['feature_list']:
                                        featureObj = geneObj[featureID]
                                        for key, value in featureObj['attributes'].items():
                                                if key in attributeList:
                                                        # Ensure that this attribute works
                                                        try:
                                                                float(value)
                                                        except:
                                                                print('gff3_index_candidates: ' + str(value) + ' was found paired to attribute "' + key + '" but is not capable of float conversion. This attribute is thus not suitable for cutoff criterion checking; fix your GFF3 or do not specify this attribute.')
                                                                quit()
                                                        # Retrieve the index of this attribute
                                                        attribIndex = attributeList.index(key)
                                                        # Check against cutoff
                                                        if eval(str(value) + directionList[attribIndex] + str(cutoffList[attribIndex])):
                                                                cutoffCheck[attribIndex] = 1
                        # Handle outputList behaviour
                        if sum(cutoffCheck) == len(attributeList):
                                if behaviour.lower() == 'main' and geneObj['attributes']['ID'] not in outputList:
                                        outputList.append(geneObj['attributes']['ID'])
                                elif behaviour.lower() == 'feature':
                                        outputList.append(feature)
        return outputList

def gff3_index_intron_sizedict(gff3Object, behaviour):
        # Setup
        intronSize = {}
        # Ensure that behaviour is sensible                                                             # directionList should be a list like ['<', '>', '<=', '>=', '=='] to indicate how the attribute's value should be held up against the cutoff
        if behaviour.lower() not in ['main', 'feature']:
                print('gff3_index_intron_sizedict: behaviour must be specified as "main" or "feature"; fix the code for this section.')
                quit()
        # Main function
        for mainType in gff3Object.id_values['main'].keys():
                for geneID in gff3Object.id_values['main'][mainType]:
                        geneObj = gff3Object.index_dict[geneID]
                        # Validate that all features contain exon values and hence may contain introns
                        skip = False
                        for feature in geneObj['feature_list']:
                                if 'exon' not in geneObj[feature]:
                                        skip = True
                                        break
                        if skip == True:
                                continue
                        # Depending on behaviour, loop through each isoform/feature associated with the geneObj or just look at the longest "representative"
                        if behaviour.lower() == 'main':
                                featureList = longest_iso(geneObj)
                        else:
                                featureList = geneObj['feature_list']
                        # Loop through relevant feature(s)
                        intronList = []
                        for feature in featureList:
                                # Determine intron sizes from exon coords
                                intronLens = pair_coord_introns(geneObj[feature]['exon']['coords'])
                                if intronLens == []:
                                        intronLens = [0]
                                intronList.append(intronLens)
                        # Add values to our intronSize dict depending on behaviour
                        if behaviour.lower() == 'main':
                                intronSize[geneID] = intronList[0]
                        elif behaviour.lower() == 'feature':
                                for i in range(len(featureList)):
                                        intronSize[featureList[i]] = intronList[i]
        return intronSize

## Retrieve/remove function
def gff3_retrieve_remove_tofile(gff3IndexDict, outputFileName, idList, mode, behaviour):
        # Ensure mode value makes sense
        if mode.lower() not in ['retrieve', 'remove']:
                print('gff3_retrieve_remove_tofile: Input mode value is not "retrieve" or "remove" but is instead "' + str(mode) + '".')
                print('Fix the code for this section.')
                quit()
        # Ensure behaviour value makes sense
        if behaviour.lower() not in ['main', 'feature']:
                print('gff3_retrieve_remove_tofile: Input behaviour value is not "main" or "feature" but is instead "' + str(behaviour) + '".')
                print('Fix the code for this section.')
                quit()
        # Main function
        with open(outputFileName, 'w') as fileOut:
                # Generate our list for iteration (ensuring that genes come first)
                iterList = []
                if 'gene' in gff3IndexDict['idValues']['main']:
                        iterList.append(gff3IndexDict['idValues']['main']['gene'])
                for key in gff3IndexDict['idValues']['main'].keys():
                        if key != 'gene':
                                iterList.append(gff3IndexDict['idValues']['main'][key])
                # Iterate through main features and determine if they are being written to file
                for key in gff3IndexDict['geneValues']:
                        value = gff3IndexDict[key]
                        # Check if relevant sequence details are within our idList
                        found = []
                        if key in idList:                       # Checking gene ID here
                                found = True
                        elif value['contig_id'] in idList:      # Checking contig ID here
                                found = True
                        elif value['source'] in idList:         # Checking source here
                                found = True
                        elif value['orientation'] in ['+', '-'] and value['orientation'] in idList:     # Checking orientation here
                                found = True                    # We want this extra check for orientation since it can be '.' in some GFF3s and this might conflict with removing source columns with '.'
                        elif value['feature_type'] in idList:   # Checking feature type here
                                found = True
                        else:
                                for mrna in value['feature_list']:                      # Checking subfeature ID here [they won't always be mRNAs but it's just how I wrote the variable]
                                        if mrna in idList:
                                                found.append(mrna)
                                        elif value[mrna]['feature_type'] in idList:     # Checking subfeature type here
                                                found.append(mrna)
                                        elif value[mrna]['source'] in idList:           # Checking subfeature source here
                                                found.append(mrna)
                        # If we find all subfeatures, make our found == True so we know we're looking at the whole gene obj
                        if type(found) == list:
                                if len(found) == len(value['feature_list']):    # If these lengths are equivalent, we know that we found all features
                                        found = True
                        # Write (or don't write) to file depending on mode setting
                        if mode.lower() == 'retrieve' and (found == True or (found != [] and behaviour.lower() == 'main')):
                                fileOut.write(''.join(value['lines'][0]))
                                fileOut.write(''.join(value['lines'][1]))
                                fileOut.write(''.join(value['lines'][2]))
                        elif mode.lower() == 'retrieve' and found == []:
                                continue
                        elif mode.lower() == 'remove' and found == []:
                                fileOut.write(''.join(value['lines'][0]))
                                fileOut.write(''.join(value['lines'][1]))
                                fileOut.write(''.join(value['lines'][2]))
                        elif mode.lower() == 'remove' and (found == True or (found != [] and behaviour.lower() == 'main')):
                                continue
                        else:
                                # Retrieve relevant header lines for scenarios where only some subfeatures were found
                                newHeader = []
                                for line in value['lines'][0]:
                                        mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # Since we only store header lines with known format, we know exactly what we're dealing with here
                                        if (mrnaID in found and mode.lower() == 'retrieve') or (mrnaID not in found and mode.lower() == 'remove'):
                                                newHeader.append(line)
                                # Retrieve relevant footer lines
                                newFooter = []
                                for line in value['lines'][2]:
                                        mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # Since we only store header lines with known format, we know exactly what we're dealing with here
                                        if (mrnaID in found and mode.lower() == 'retrieve') or (mrnaID not in found and mode.lower() == 'remove'):
                                                newFooter.append(line)
                                # Retrieve relevant feature lines
                                newFeature = []
                                for line in value['lines'][1]:
                                        sl = line.split('\t')
                                        details = sl[8].rstrip('\r\n').split(';')
                                        idField, parentField = None, None
                                        for i in range(len(details)):
                                                if details[i].startswith('ID='):
                                                        idField = details[i][3:]
                                                elif details[i].startswith('Parent='):
                                                         parentField = details[i][7:]
                                        if ((idField in found or parentField in found) and mode.lower() == 'retrieve') or (idField not in found and parentField not in found and mode.lower() == 'remove'):
                                                newFeature.append(line)
                                # Update our first gene line to reflect potential new start, stop coordinates
                                coords = []
                                for i in range(1, len(newFeature)):
                                        sl = newFeature[i].split('\t')
                                        start, stop = int(sl[3]), int(sl[4])
                                        coords += [start, stop]
                                newGeneLine = newFeature[0].split('\t')
                                newGeneLine[3], newGeneLine[4] = str(min(coords)), str(max(coords))
                                newFeature[0] = '\t'.join(newGeneLine)
                                # Write to file
                                fileOut.write(''.join(newHeader))
                                fileOut.write(''.join(newFeature))
                                fileOut.write(''.join(newFooter))

def gff3_retrieve_remove_tolist(gff3File, idList, identifiers, behaviour):
        # Setup
        outList = []
        # Ensure behaviour value makes sense
        if behaviour.lower() not in ['retrieve', 'remove']:
                print('gff3_cull_output: Input behaviour value is not "retrieve" or "remove" but is instead "' + str(behaviour) + '".')
                print('Fix the code for this section.')
                quit()
        # Main function
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        geneID = None   # This lets us perform a check to ensure we pulled out a gene ID
                        sl = line.split()
                        # Skip filler lines
                        if line == '\n' or line == '\r\n':
                                continue
                        # Handle comment lines
                        elif '#' in line:
                                for section in sl:
                                        for ident in identifiers:               # Identifiers should be a list that contains values that will occur in every line that contains values we want to retrieve/remove
                                                if ident in section:            # By default this should be '.model' or '.path'; .model appears in every '#' and full GFF3 line of PASA-formatted files; .path is for GMAP
                                                        geneID = section.rstrip(',')    # PASA-formatted comments are written human-like and contain commas; we want to remove these
                        # Handle gene annotation lines
                        elif sl[2] == 'gene':
                                gffComment = sl[8].split(';')
                                for section in gffComment:
                                        if section.startswith('ID='):
                                                geneID = section[3:]                    # Skip the ID= at start
                                                break
                        # Handle gmap_gene_find lines specifically
                        elif sl[1] == 'gmap_gene_find':
                                gffComment = sl[8].split(';')
                                for section in gffComment:
                                        if section.startswith('Parent='):
                                                geneID = section[7:].strip('\r\n')      # Skip the Parent= at start and remove newline and return characters
                                                break
                        else:
                                gffComment = sl[8].split(';')
                                for section in gffComment:
                                        if section.startswith('Parent='):
                                                geneID = section[7:].strip('\r\n')      # Skip the Parent= at start and remove newline and return characters
                                                break
                        # Get the contig ID if applicable
                        if '#' not in line:
                                contigID = sl[0]
                        # Decide if we're writing this non-gene line (e.g., rRNA or tRNA annotations) to file based on behaviour
                        if geneID == None:
                                if contigID == None:                    # If this is a pure comment line without gene details in it (e.g., PASA head/foot comments) then just write it to file
                                        outList.append(line)
                                elif behaviour.lower() == 'retrieve':   # If we get here, it's not a pure comment line; it is likely a rRNA or tRNA annotation line (or just something non-genic)
                                        if contigID in idList:          # In this case we want to only retain (or below, remove) if there is a contigID match
                                                outList.append(line)
                                elif behaviour.lower() == 'remove':
                                        if contigID not in idList:
                                                outList.append(line)
                        # Decide if we're holding this line in a list based on behaviour
                        elif behaviour.lower() == 'retrieve':
                                if geneID in idList:
                                        outList.append(line)
                        elif behaviour.lower() == 'remove':
                                if geneID not in idList:
                                        outList.append(line)
        return outList

## GFF3
def gff3_index_contig_reorder(gff3Index, fastaFile, locationStart):     # see mtdna_gff2_order.py for an example of this function
        # Define functions integral to this one
        def coord_overflow_invert(start, end, contigLen):
                if start < 1:
                        start = contigLen - abs(start)
                if end < 1:
                        end = contigLen - abs(end)
                return start, end
        # Parse FASTA file and ensure it is sensible
        records = SeqIO.to_dict(SeqIO.parse(open(fastaFile, 'r'), 'fasta'))
        if len(records) != 1:
                print('Incompatible input detected: the FASTA file has more than 1 contig value within it.')
                print('This program is designed to handle MtDNA annotations which should occur on a single contig.')
                print('Program will exit now.')
                quit()
        if gff3Index['contigValues'][0] not in records:
                print('Incompatible input detected: the FASTA file does not have the same contig ID present as the annotation file.')
                print('Annotation contig ID = "' + gff3Index['contigValues'][0] + '"... FASTA contig ID = "' + list(records.keys())[0] + '"')
                print('Program will exit now.')
                quit()
        contigID = gff3Index['contigValues'][0]
        contigLen = len(records[contigID])
        # Validate that locationStart value is sensible
        if locationStart in gff3Index:
                locationStart = gff3Index[locationStart]['coords'][0]
                print('Provided location start value corresponds to a feature ID; all coordinates will be rearranged so ' + str(locationStart) + ' becomes the first base.')
        else:
                try:
                        locationStart = int(locationStart)
                        if locationStart < 0:
                                print('If the start location (-l) value is not a gene ID, it is assumed to be an integer.')
                                print('This has occurred here, however, this value cannot be an integer less than 0.')
                                print('Fix your inputs and try again; program will exit now.')
                                quit()
                        elif locationStart > contigLen:
                                print('If the start location (-l) value is not a gene ID, it is assumed to be an integer.')
                                print('This has occurred here, however, this value cannot be greater than the contig length (which is == ' + str(contigLen))
                                print('Fix your inputs and try again; program will exit now.')
                                quit()
                except:
                        print('gff3_contig_reorder: The provided start location "' + locationStart + '" is not correct.')
                        print('It is neither capable of conversion to integer nor is it a gene ID or feature ID.')
                        print('Fix your inputs and try again; program will exit now.')
                        quit()
        # Update all coordinates for values in gff3Index
        '''Note that GFF3 positions start at 1, which means that if our locationStart is 100,
        we need to subtract 99 from all positions to get their new location. Features which are 
        earlier than 100 also have 99 subtracted from them which should result in negative numbers.
        However, they have to go through 0 on the way to become negative integers, and in this situation
        0 is equivalent to contigLen. This gives the formula for negative integers of:
                newLocation = contigLen - (oldLocation - locationStart)'''
        for primaryID in gff3Index['primaryValues']:
                '''Note here that we need to take a ground-up approach to the coordinate revisioning. If we go
                top-down with genes that have introns, we might not allow certain genes to be broken up correctly
                across introns. This might normally be a problem, but since MITOS2 doesn't order its genes properly
                anyway, it's something that the user will need to sort out themself manually (sorry, can't automate it...)'''
                # Ground-up loop
                for subfeature in gff3Index[primaryID]['feature_list']:
                        # Update secondary subfeature values
                        for subfeatType in gff3Index[primaryID][subfeature]['feature_list']:
                                for i in range(len(gff3Index[primaryID][subfeature][subfeatType]['coords'])):
                                        newStart = gff3Index[primaryID][subfeature][subfeatType]['coords'][i][0] - (locationStart - 1)
                                        newEnd = gff3Index[primaryID][subfeature][subfeatType]['coords'][i][1] - (locationStart - 1)
                                        if (newStart < 1 and newEnd >= 1) or newStart > contigLen or newEnd > contigLen:
                                                print('The new start location results in an exon being split into two fragments at the start and end of the MtDNA genome')
                                                print('Feature in question = "' + primaryID + '"')
                                                print('This isn\'t a good way to present the annotation, so I\'ve not been coded to handle this scenario.')
                                                print('Choose a start location that doesn\'t result in such fragmentation and try again; program will exit now.')
                                                quit()
                                        newStart, newEnd = coord_overflow_invert(newStart, newEnd, contigLen)
                                        gff3Index[primaryID][subfeature][subfeatType]['coords'][i] = [newStart, newEnd]
                        # Update primary subfeature values based on secondary subfeatures
                        secondarySubCoords = [x for coordPair in gff3Index[primaryID][subfeature]['exon']['coords'] for x in coordPair]
                        gff3Index[primaryID][subfeature]['coords'] = [min(secondarySubCoords), max(secondarySubCoords)]
                # Update primary feature values based on primary subfeatures
                subfeatCoords = [x for subfeature in gff3Index[primaryID]['feature_list'] for x in gff3Index[primaryID][subfeature]['coords']]
                gff3Index[primaryID]['coords'] = [min(subfeatCoords), max(subfeatCoords)]
        return gff3Index, locationStart

## NCLS-related functions
def gff3_parse_ncls(gff3File, featureTypes):
        gff3Loc = {}
        starts = []
        ends = []
        ids = []
        ongoingCount = 0
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        # Skip unneccessary lines
                        if line.startswith('#') or line == '\n' or line == '\r\n':
                                continue
                        sl = line.split('\t')
                        if len(sl) < 3:
                                continue
                        # Skip lines that aren't being stored
                        if sl[2] not in featureTypes:
                                continue
                        # Get details from line including start, stop, and orientation
                        contigID = sl[0]
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        orient = sl[6]
                        details = sl[8].split(';')
                        detail_dict = {}
                        for i in range(len(details)):
                                if details[i] == '' or details[i] == '\n':
                                        continue
                                split_details = details[i].split('=', maxsplit=1)
                                detail_dict[split_details[0]] = split_details[1].rstrip('\r\n')
                        if 'ID' not in detail_dict:      # Don't index things which lack IDs; these might include things like TAIR9's 'protein' features
                                continue
                        # Add to our NCLS
                        starts.append(contigStart)
                        ends.append(contigStop+1)       # NCLS indexes 0-based like a range (up to but not including end), so +1 to make this more logically compliant with gff3 1-based system.
                        ids.append(ongoingCount)
                        gff3Loc[ongoingCount] = [contigStart, contigStop, orient, detail_dict['ID'], contigID]
                        ongoingCount += 1
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        return ncls, gff3Loc

def ncls_finder(ncls, locDict, start, stop):
        overlaps = ncls.find_overlap(start, stop+1)             # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        dictEntries = []
        for result in overlaps:
                dictEntries.append(locDict[result[2]])
        dictEntries = copy.deepcopy(dictEntries)                # Any time we're deleting things from a section of a dictionary we need to build a deepcopy to keep the original dictionary intact.
        # Return list
        return dictEntries                                      # This list will consist of all overlaps within the same coordinate range; these may be across multiple contigs/features, hence the need for narrowing

def ncls_feature_narrowing(nclsEntries, featureID, featureIndex):       # This code will narrow the results of ncls_finder to objects with a specific featureID
        for k in range(len(nclsEntries)-1, -1, -1):                     # featureIndex should correspond to the index of the feature in the dictEntries sublist objects output by ncls_finder
                if nclsEntries[k][featureIndex] != featureID:           # See gff3_merge or gmap_gene_find for an example of this code
                        del nclsEntries[k]
        return nclsEntries

## GFF3 transposon structure manipulation
def longest_iso(geneDictObj):
        longestMrna = ['', 0]           # We pick out the representative gene based on length. If length is identical, we'll end up picking the entry listed first in the gff3 file since our > condition won't be met. I doubt this will happen much or at all though.
        for mrna in geneDictObj['feature_list']:
                mrnaLen = 0
                if 'CDS' in geneDictObj[mrna]:
                        coordList = geneDictObj[mrna]['CDS']['coords']
                else:
                        coordList = geneDictObj[mrna]['exon']['coords']
                for coord in coordList:
                        mrnaLen += (int(coord[1]) - int(coord[0]) + 1)
                if mrnaLen > longestMrna[1]:
                        longestMrna = [mrna, mrnaLen]
        mrnaList = [longestMrna[0]]
        return mrnaList

def pair_coord_introns(inputList):
        introns = []
        for i in range(0, len(inputList)-1):
                exon1 = inputList[i]
                exon2 = inputList[i+1]
                start = min(int(exon1[1]), int(exon2[0])) + 1   # +1 gets our first bp of intron.
                end = max(int(exon1[1]), int(exon2[0])) - 1     # -1 does the same on the opposite border of the intron.
                intLen = end - start + 1                        # +1 as above scenario in pair_coord_exons
                introns.append(intLen)
        return introns

## GFF3 isoform handling operations
def gff3_merge_and_isoclust(mainGff3Lines, newGff3Lines, isoformDict, excludeList, outFileName):        # See gff3_merge.py for example of this function
        # Set up
        processedPaths = []
        # Main function
        with open(outFileName, 'w') as fileOut:
                # Merging isoform clusters
                for key in mainGff3Lines['geneValues']:
                        if key in isoformDict:
                                # Write opening comments for main gene
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][0]))
                                # Loop into associated isoforms and write header comments (if relevant) & hold onto coordinates
                                mrnaCoords = []
                                for mrna in isoformDict[key]:
                                        # Get this mRNA's header line specifically (if it has one)
                                        mrnaHead = None
                                        for line in newGff3Lines[mrna]['lines'][0]:
                                                if mrna in line or newGff3Lines[mrna]['attributes']['ID'] in line:              # i.e., if the mRNA or gene ID is in the line
                                                        mrnaHead = line.replace(newGff3Lines[mrna]['attributes']['ID'], mrna)   # if the gene ID is in the line, we want it to become the mRNA ID
                                        if mrnaHead != None:
                                                fileOut.write(mrnaHead)
                                        # Get the mRNA coordinates
                                        mrnaCoords.append(newGff3Lines[mrna][mrna]['coords'])
                                        processedPaths.append(mrna)
                                # Get minimum/maximum coordinates for the mRNAs being clustered into this gene as isoforms
                                minMrna = None
                                maxMrna = None
                                for coord in mrnaCoords:
                                        if minMrna == None:
                                                minMrna, maxMrna = coord[0], coord[1]
                                        if coord[0] < minMrna:
                                                minMrna = coord[0]
                                        if coord[1] > maxMrna:
                                                maxMrna = coord[1]
                                # Update our gene start/stop coordinates if relevant
                                mainGff3Lines[key]['coords'] = [min(mainGff3Lines[key]['coords'][0], minMrna), max(mainGff3Lines[key]['coords'][1], maxMrna)]
                                newGeneLine = mainGff3Lines[key]['lines'][1][0].split('\t')
                                newGeneLine[3], newGeneLine[4] = list(map(str, mainGff3Lines[key]['coords']))
                                mainGff3Lines[key]['lines'][1][0] = '\t'.join(newGeneLine)
                                # Write main gene and mRNA lines
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][1]))
                                # Loop into associated isoforms and write their mRNA lines
                                for mrna in isoformDict[key]:
                                        # Retrieve the lines specifically mapping to this mRNA
                                        mrnaLines = []
                                        for line in newGff3Lines[mrna]['lines'][1][1:]:                 # Skip the first gene line
                                                if 'ID=' + mrna in line or 'Parent=' + mrna in line:    # This is a simple way to check if we have the correct value in our attributes fields when parsing the line as a string directly
                                                        mrnaLines.append(line)
                                        # Write lines to file after editing their attributes field appropriately
                                        for line in mrnaLines:
                                                sl = line.rstrip('\\n').split('\t')                      # Need to strip the newline character off so we can work with attributes at the end of the line; we'll add this back in later
                                                attributes = sl[8].split(';')
                                                for i in range(len(attributes)):
                                                        if attributes[i].startswith('Parent='):
                                                                if sl[2] == 'mRNA':
                                                                        attributes[i] = 'Parent=' + key  # For mRNA lines, the parent is the main gene ID which is represented by 'key' currently
                                                                else:
                                                                        attributes[i] = 'Parent=' + mrna # For all other feature types (e.g., exon, CDS) the parent is the mRNA ID which is represented by 'mrna' currently
                                                attributes = ';'.join(attributes)
                                                sl[8] = attributes
                                                line = '\t'.join(sl)
                                                fileOut.write(line + '\n')
                                # Write closing comments for main gene
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][2]))
                                # Loop into associated isoforms and write their closing comments
                                for mrna in isoformDict[key]:
                                        # Get this mRNA's footer line specifically (if it has one)
                                        mrnaFoot = None
                                        for line in newGff3Lines[mrna]['lines'][2]:
                                                if mrna in line or newGff3Lines[mrna]['attributes']['ID'] in line:              # i.e., if the mRNA or gene ID is in the line
                                                        mrnaFoot = line.replace(newGff3Lines[mrna]['attributes']['ID'], key)    # Similar to the header comment, we need to replace the original gene ID; this time it's with the new gene ID
                                        if mrnaFoot != None:
                                                fileOut.write(mrnaFoot)
                        # Write genes without clustered isoforms to file directly
                        else:
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][0]))
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][1]))
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][2]))
                # Drop any new values not clustered as isoforms into the file
                for geneID in newGff3Lines['geneValues']:
                        # Figure out which of this gene's mRNAs were not already clustered as isoforms
                        nonisoMrnas = []
                        for mrnaID in newGff3Lines[geneID]['feature_list']:
                                if mrnaID not in processedPaths and mrnaID not in excludeList:
                                        nonisoMrnas.append(mrnaID)
                        if nonisoMrnas == []:
                                continue
                        # If no changes are required for this gene, write it to file like normal [If these sets are equivalent we didn't grab anything from this gene for isoform clustering/exclude any mRNAs and don't need to bother with more elaborate handling]
                        if set(nonisoMrnas) == set(newGff3Lines[geneID]['feature_list']):
                                fileOut.write(''.join(newGff3Lines[geneID]['lines'][0]))
                                fileOut.write(''.join(newGff3Lines[geneID]['lines'][1]))
                                fileOut.write(''.join(newGff3Lines[geneID]['lines'][2]))
                        else:
                                # Write header lines for this gene's mRNAs
                                mrnaHeads = []
                                for mrnaID in nonisoMrnas:
                                        mrnaHead = None
                                        for line in newGff3Lines[geneID]['lines'][0]:
                                                if mrnaID in line or geneID in line:
                                                        mrnaHead = line
                                        if mrnaHead != None:
                                                # Handle gene ID duplication
                                                if geneID in mainGff3Lines['idValues'][0]:
                                                        if geneID in mrnaHead:
                                                                mrnaHead = mrnaHead.replace(geneID, geneID + '_gff3_merge_separated')
                                                        if mrnaID in mrnaHead:
                                                                mrnaHead = mrnaHead.replace(mrnaID, mrnaID + '_gff3_merge_separated')
                                                if mrnaHead not in mrnaHeads:   # We need to do this since we're looking for mRNA OR gene IDs in the header comment; this is necessary for GGF but might cause redundancy with other GFF3 formats
                                                        mrnaHeads.append(mrnaHead)
                                fileOut.write(''.join(mrnaHeads))
                                # Get minimum/maximum coordinates for the mRNAs being clustered into this gene as isoforms
                                minMrna = None
                                maxMrna = None
                                for mrnaID in nonisoMrnas:
                                        coord = newGff3Lines[geneID][mrnaID]['coords']
                                        if minMrna == None:
                                                minMrna, maxMrna = coord[0], coord[1]
                                        if coord[0] < minMrna:
                                                minMrna = coord[0]
                                        if coord[1] > maxMrna:
                                                maxMrna = coord[1]
                                # Update our gene start/stop coordinates if relevant
                                newGff3Lines[geneID]['coords'] = [minMrna, maxMrna]
                                newGeneLine = newGff3Lines[geneID]['lines'][1][0].split('\t')
                                newGeneLine[3], newGeneLine[4] = list(map(str, newGff3Lines[geneID]['coords']))
                                if geneID in mainGff3Lines['idValues'][0]:
                                        # Handle gene ID duplication
                                        newGeneLine[8] = newGeneLine[8].replace('ID=' + geneID, 'ID=' + geneID + '_gff3_merge_separated')
                                newGff3Lines[geneID]['lines'][1][0] = '\t'.join(newGeneLine)
                                # Write main gene and mRNA lines
                                fileOut.write(''.join(newGff3Lines[geneID]['lines'][1][0]))
                                for mrnaID in nonisoMrnas:
                                        for line in newGff3Lines[geneID]['lines'][1][1:]:                       # Skip the first gene line
                                                if 'ID=' + mrnaID in line or 'Parent=' + mrnaID in line:        # This is a simple way to check if we have the correct value in our attributes fields when parsing the line as a string directly
                                                        # Handle gene ID duplication
                                                        if geneID in mainGff3Lines['idValues'][0] and 'ID=' + mrnaID in line:   # We only need to change the parent ID for mRNA lines and only when we're dealing with duplicate gene ID
                                                                line = line.replace('Parent=' + geneID, 'Parent=' + geneID + '_gff3_merge_separated')
                                                        fileOut.write(line)
                                # Write footer lines for this gene's mRNAs
                                mrnaFoots = []
                                for mrnaID in nonisoMrnas:
                                        # Get this mRNA's footer line specifically (if it has one)
                                        mrnaFoot = None
                                        for line in newGff3Lines[geneID]['lines'][2]:
                                                if mrnaID in line or geneID in line:
                                                        mrnaFoot = line
                                        if mrnaFoot != None:
                                                # Handle gene ID duplication
                                                if geneID in mainGff3Lines['idValues'][0]:
                                                        if geneID in mrnaFoot:
                                                                mrnaFoot = mrnaFoot.replace(geneID, geneID + '_gff3_merge_separated')
                                                        if mrnaID in mrnaFoot:
                                                                mrnaFoot = mrnaFoot.replace(mrnaID, mrnaID + '_gff3_merge_separated')
                                                if mrnaFoot not in mrnaFoots:
                                                        mrnaFoots.append(mrnaFoot)
                                fileOut.write(''.join(mrnaFoots))
                # Write non-gene lines to file if relevant
                valueList = []
                for key in mainGff3Lines['idValues']['main'].keys():
                        if key != 'gene':
                                valueList.append(mainGff3Lines['idValues']['main'][key])
                for value in valueList:
                        for key in value:
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][0]))
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][1]))
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][2]))
                                excludeList.append(key)         # This helps with preventing redundancy with "note" entries like line_type == "chromosome"
                valueList = []
                for key in newGff3Lines['idValues']['main'].keys():
                        if key != 'gene':
                                valueList.append(newGff3Lines['idValues']['main'][key])
                for value in valueList:
                        for key in value:
                                found = False
                                for feature in newGff3Lines[key]['feature_list']:
                                        if feature in excludeList:
                                                found = True
                                if found == False and key not in excludeList:
                                        fileOut.write(''.join(newGff3Lines[key]['lines'][0]))
                                        fileOut.write(''.join(newGff3Lines[key]['lines'][1]))
                                        fileOut.write(''.join(newGff3Lines[key]['lines'][2]))

def edit_parent(gff3Line, parentID):
        # Handle unsplit values
        if type(gff3Line) != list:
                gff3Line = gff3Line.split('\t')
                gff3Line[-1] = gff3Line[-1].rstrip('\r\n')      # Need to make sure there isn't a new line at the end; we'll handle this in the main loop
        # Edit parent comment
        commentValues = gff3Line[8].split(';')
        for i in range(len(commentValues)):
                if commentValues[i].startswith('Parent='):
                        commentValues[i] = 'Parent=' + parentID
        commentValues = ';'.join(commentValues)
        gff3Line[8] = commentValues
        gff3Line = '\t'.join(gff3Line)
        return gff3Line

## GFF3 sequence extraction & handling
def gff3_object_sequence_extract(gff3Object, mrna, genomeRecords, seqType): # gff3Object should be produced by the Gff3 class; mrna is a string which should correspond to a subfeature key in the Gff3.index_dict; genomeRecords should be a Biopython SeqIO.parse() object of the genome's contigs
        # Setup
        cdsWarning = False
        # Ensure that seqType is sensible
        seqType = seqType.lower()
        if seqType not in ['transcript', 'cds', 'both']:
                print('gff3_index_sequence_extract: seqType value is not sensible; you need to fix the inputs to this function.')
                quit()
        # Obtain the indexed gene object
        value = gff3Object.index_dict[mrna]
        # Retrieve genomic sequence
        try:
                genomeSeq = str(genomeRecords[value['contig_id']].seq)
        except:
                print('Contig ID "' + value['contig_id'] + '" is not present in your FASTA file; mRNA ID "' + mrna + '" cannot be handled.')
                print('This represents a major problem with your inputs, so I\'m going to stop processing now. Make sure you are using the correct FASTA file and try again.')
                quit()
        # Sort coords lists for consistency [this can be relevant since not all GFF3's are ordered equivalently]
        ## Exon coords
        if value[mrna]['orientation'] == '+':
                value[mrna]['exon']['coords'].sort(key = lambda x: (int(x[0]), int(x[1])))
        elif value[mrna]['orientation'] == '-':
                value[mrna]['exon']['coords'].sort(key = lambda x: (-int(x[0]), -int(x[1])))
        ## CDS coords
        if 'CDS' in value[mrna]: # This check here (and below) is a way of ensuring that we only produce CDS outputs for features that are annotated as coding in the GFF3 file
                cdsSort = list(zip(value[mrna]['CDS']['coords'], value[mrna]['CDS']['frame']))
                if value[mrna]['orientation'] == '+':
                        cdsSort.sort(key = lambda x: (int(x[0][0]), int(x[0][1])))
                        value[mrna]['CDS']['coords'] = [coord for coord,frame in cdsSort]
                        value[mrna]['CDS']['frame'] = [frame for coord,frame in cdsSort]
                elif value[mrna]['orientation'] == '-':
                        cdsSort.sort(key = lambda x: (-int(x[0][0]), -int(x[0][1])))
                        value[mrna]['CDS']['coords'] = [coord for coord,frame in cdsSort]
                        value[mrna]['CDS']['frame'] = [frame for coord,frame in cdsSort]
                else:
                        print(mrna + ' lacks proper orientation specification within GFF3 (it is == "' + str(value[mrna]['orientation']) + '"; this may result in problems.')
        elif cdsWarning == False and seqType == 'both':
                print('Warning: there are \'gene\' features which lack CDS subfeatures; your .trans file will contain more entries than .aa or .nucl files.')
                cdsWarning = True
        # Reverse the coord lists if we're looking at a '-' model so we start at the 3' end of the gene model
        if value['orientation'] == '-':
                value[mrna]['exon']['coords'].reverse()
                if 'CDS' in value[mrna]:
                        value[mrna]['CDS']['coords'].reverse()
        # Join sequence segments
        if seqType == 'transcript' or seqType == 'both':
                transcript = ''
                for coord in value[mrna]['exon']['coords']:
                        segment = genomeSeq[int(coord[0])-1:int(coord[1])] # Make it 1-based by -1 to the first coordinate
                        transcript += segment
                # Reverse comp if necessary
                if value['orientation'] == '-':
                        transcript = reverse_comp(transcript)
        if seqType == 'cds' or seqType == 'both':
                cds = None      # This lets us return something when CDS doesn't exist
                if 'CDS' in value[mrna]:
                        cds = ''
                        for coord in value[mrna]['CDS']['coords']:
                                segment = genomeSeq[int(coord[0])-1:int(coord[1])] # Make it 1-based by -1 to the first coordinate
                                cds += segment
                        # Reverse comp if necessary
                        if value['orientation'] == '-':
                                cds = reverse_comp(cds)
        # Return transcript and/or CDS sequence
        if seqType == 'transcript':
                return transcript
        elif seqType == 'cds':
                return cds
        elif seqType == 'both':
                return transcript, cds
