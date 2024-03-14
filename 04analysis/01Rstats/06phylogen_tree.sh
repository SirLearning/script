# build phylogenetic tree using IQ-TREE
iqtree -s RT.aln -bb 1000 \
    -pre Ty1_copia.iqtree 1>iqtree.log 2>&1 &
# polish the tree using iTOL: https://itol.embl.de/