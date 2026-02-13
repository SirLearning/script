import sys

# Static Data Strings
_POSITION_FILE_S = (
    "0\t0\t480980714\tchrUn\t0\t480980714\n"
    "1\t0\t471304005\tchr1A\t0\t471304005\n"
    "2\t0\t122798051\tchr1A\t471304005\t594102056\n"
    "3\t0\t438720154\tchr1B\t0\t438720154\n"
    "4\t0\t251131716\tchr1B\t438720154\t689851870\n"
    "5\t0\t452179604\tchr1D\t0\t452179604\n"
    "6\t0\t43273582\tchr1D\t452179604\t495453186\n"
    "7\t0\t462376173\tchr2A\t0\t462376173\n"
    "8\t0\t318422384\tchr2A\t462376173\t780798557\n"
    "9\t0\t453218924\tchr2B\t0\t453218924\n"
    "10\t0\t348037791\tchr2B\t453218924\t801256715\n"
    "11\t0\t462216879\tchr2D\t0\t462216879\n"
    "12\t0\t189635730\tchr2D\t462216879\t651852609\n"
    "13\t0\t454103970\tchr3A\t0\t454103970\n"
    "14\t0\t296739669\tchr3A\t454103970\t750843639\n"
    "15\t0\t448155269\tchr3B\t0\t448155269\n"
    "16\t0\t382674495\tchr3B\t448155269\t830829764\n"
    "17\t0\t476235359\tchr3D\t0\t476235359\n"
    "18\t0\t139317064\tchr3D\t476235359\t615552423\n"
    "19\t0\t452555092\tchr4A\t0\t452555092\n"
    "20\t0\t292033065\tchr4A\t452555092\t744588157\n"
    "21\t0\t451014251\tchr4B\t0\t451014251\n"
    "22\t0\t222603248\tchr4B\t451014251\t673617499\n"
    "23\t0\t451004620\tchr4D\t0\t451004620\n"
    "24\t0\t58852447\tchr4D\t451004620\t509857067\n"
    "25\t0\t453230519\tchr5A\t0\t453230519\n"
    "26\t0\t256543224\tchr5A\t453230519\t709773743\n"
    "27\t0\t451372872\tchr5B\t0\t451372872\n"
    "28\t0\t261776885\tchr5B\t451372872\t713149757\n"
    "29\t0\t451901030\tchr5D\t0\t451901030\n"
    "30\t0\t114179647\tchr5D\t451901030\t566080677\n"
    "31\t0\t452440856\tchr6A\t0\t452440856\n"
    "32\t0\t165638404\tchr6A\t452440856\t618079260\n"
    "33\t0\t452077197\tchr6B\t0\t452077197\n"
    "34\t0\t268911281\tchr6B\t452077197\t720988478\n"
    "35\t0\t450509124\tchr6D\t0\t450509124\n"
    "36\t0\t23083594\tchr6D\t450509124\t473592718\n"
    "37\t0\t450046986\tchr7A\t0\t450046986\n"
    "38\t0\t286659250\tchr7A\t450046986\t736706236\n"
    "39\t0\t453822637\tchr7B\t0\t453822637\n"
    "40\t0\t296797748\tchr7B\t453822637\t750620385\n"
    "41\t0\t453812268\tchr7D\t0\t453812268\n"
    "42\t0\t184873787\tchr7D\t453812268\t638686055\n"
    "43\t0\t452528\tMit\t0\t452528\n"
    "44\t0\t134545\tChl\t0\t134545"
)

_CENTROMERE_FILE_S = (
    "chr1A\t210200000\t215800000\n"
    "chr1B\t237700000\t243500000\n"
    "chr1D\t166200000\t173800000\n"
    "chr2A\t326300000\t327000000\n"
    "chr2B\t344400000\t351300000\n"
    "chr2D\t264400000\t272500000\n"
    "chr3A\t316900000\t319900000\n"
    "chr3B\t345800000\t347000000\n"
    "chr3D\t237100000\t243200000\n"
    "chr4A\t264100000\t267900000\n"
    "chr4B\t303900000\t304400000\n"
    "chr4D\t182300000\t188200000\n"
    "chr5A\t252500000\t255100000\n"
    "chr5B\t198900000\t202500000\n"
    "chr5D\t185600000\t188700000\n"
    "chr6A\t283300000\t288700000\n"
    "chr6B\t323000000\t327500000\n"
    "chr6D\t211900000\t217400000\n"
    "chr7A\t360200000\t363800000\n"
    "chr7B\t308000000\t310100000\n"
    "chr7D\t336300000\t341700000"
)

# Global Maps
_chrIDChromosomeMap = {}
_chrIDLengthMap = {}
_chromosomeHalfLengthMap = {}
_centromereStartMap = {}
_centromereEndMap = {}


def _build_maps():
    global _chrIDChromosomeMap, _chrIDLengthMap, _chromosomeHalfLengthMap
    global _centromereStartMap, _centromereEndMap

    # 1. Parse Position File
    lines = _POSITION_FILE_S.split('\n')
    for i, line in enumerate(lines):
        if not line.strip(): continue
        temp = line.split('\t')
        
        chrID = int(temp[0])
        length = int(temp[2])
        chrom = temp[3].replace("chr", "")
        
        _chrIDChromosomeMap[chrID] = chrom
        _chrIDLengthMap[chrID] = length
        
        # In Java: if (i%2 != 0) continue; 
        # Java loop was 0-indexed on the array of lines.
        # Lines 1, 2, 3... correspond to indices 0, 1, 2...
        # So even indices (0, 2...) are processed for half length.
        if i % 2 == 0:
            _chromosomeHalfLengthMap[chrom] = length

    # 2. Parse Centromere File
    lines = _CENTROMERE_FILE_S.split('\n')
    for line in lines:
        if not line.strip(): continue
        # Use split() to handle whitespace/tabs generically like PStringUtils.fastSplit often does, 
        # or explicit split('\t') if consistent. Java code used fastSplit which usually handles whitespace.
        parts = line.split() 
        if len(parts) >= 3:
            chrom = parts[0].replace("chr", "")
            _centromereStartMap[chrom] = int(parts[1])
            _centromereEndMap[chrom] = int(parts[2])

# Initialize
_build_maps()

# ==================================================================================
# Public API
# ==================================================================================

def get_chromosomes():
    """Return all sorted chromosomes (e.g. ['1A', '1B'...])"""
    chroms = list(_centromereStartMap.keys())
    chroms.sort()
    return chroms

def get_chr_ids():
    """Return all sorted chrIDs (1-42)"""
    ids = list(_chrIDChromosomeMap.keys())
    ids.sort()
    return ids

def get_chr_ids_of_subgenome(sub):
    """Return chrIDs for a subgenome ('A', 'B', or 'D')"""
    chromosomes = get_chromosomes_of_subgenome(sub)
    chr_ids = []
    for chrom in chromosomes:
        chr_ids.append(get_chr_id(chrom, 1))
        # Java used Integer.MAX_VALUE to trigger logic for second part
        chr_ids.append(get_chr_id(chrom, 2147483647)) 
    chr_ids.sort()
    return chr_ids

def get_chr_ids_of_subgenome_a():
    return get_chr_ids_of_subgenome("A")

def get_chr_ids_of_subgenome_b():
    return get_chr_ids_of_subgenome("B")

def get_chr_ids_of_subgenome_d():
    return get_chr_ids_of_subgenome("D")

def get_chromosomes_of_subgenome(sub):
    """Return chromosomes of a subgenome (e.g., all ending in 'A')"""
    all_chroms = get_chromosomes()
    sub_chroms = [c for c in all_chroms if c.endswith(sub)]
    sub_chroms.sort()
    return sub_chroms

def get_chromosomes_of_subgenome_a():
    return get_chromosomes_of_subgenome("A")

def get_chromosomes_of_subgenome_b():
    return get_chromosomes_of_subgenome("B")

def get_chromosomes_of_subgenome_d():
    return get_chromosomes_of_subgenome("D")

def get_chromosome_list():
    """Return a sorted list of chromosome names (Same as get_chromosomes but List in Java)"""
    return get_chromosomes()

def get_centromere_start(chromosome):
    return _centromereStartMap.get(chromosome)

def get_centromere_end(chromosome):
    return _centromereEndMap.get(chromosome)

def get_chr_id_length(chrID):
    return _chrIDLengthMap.get(chrID)

def get_chromosome_length(chromosome):
    """Return total length of a chromosome (sum of its two parts)"""
    # Part 1
    len1 = get_chr_id_length(get_chr_id(chromosome, 1))
    # Part 2
    len2 = get_chr_id_length(get_chr_id(chromosome, 2147483647))
    return len1 + len2

def get_chromosome(chrID, position=None):
    """Return the chromosome name for a given chrID (position is ignored)"""
    return _chrIDChromosomeMap.get(chrID)

def get_pos_on_chromosome(chrID, position):
    """Return the linear position on the full chromosome from a genome position of chrID"""
    # If odd (1, 3, 5...), it's the first part.
    if chrID % 2 == 1:
        return position
    else:
        # If even, it's the second part. Add length of the previous part.
        prev_len = _chrIDLengthMap.get(chrID - 1)
        return position + prev_len

def get_chr_id(chromosome, position):
    """Return chrID from position on chromosome"""
    halfLength = _chromosomeHalfLengthMap.get(chromosome, 0)
    
    # Base calculation based on first char (Genome number 1-7)
    # Java: (Integer.parseInt(chromosome.substring(0, 1))-1)*6
    # Python equivalent:
    genome_num = int(chromosome[0])
    chr_base = (genome_num - 1) * 6
    
    subgenome = chromosome[1] # 'A', 'B', 'D'
    
    val = chr_base
    if subgenome == 'A':
        val += 1
    elif subgenome == 'B':
        val += 3
    elif subgenome == 'D':
        val += 5
        
    # Check if position is in second half
    if position > halfLength:
        val += 1
        
    return val

def get_pos_on_chr_id(chromosome, position):
    """Return position of ChrID from a position on chromosome"""
    chr_val = get_chr_id(chromosome, position)
    
    # If it's the second part (even ID), subtract the length of the first part
    if chr_val % 2 == 0:
        return position - _chromosomeHalfLengthMap.get(chromosome, 0)
    else:
        return position

def get_subgenome_from_chr_id(chrID):
    """Return the subgenome of chrID"""
    chrom = get_chromosome(chrID, 1)
    if chrom:
        return chrom[-1]
    return None
