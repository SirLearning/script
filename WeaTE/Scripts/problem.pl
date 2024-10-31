sub Identifier {
##Major program, align the boundaries
chomp ($date=`date`);
print "$date\tModules 2-5: Start to analyze the structure of candidates...
\t\t\t\tThe terminal motif, TSD, boundary, orientation, age, and superfamily will be identified in this step.\n\n";
my $index=$_[0];
`perl $script_path/bin/get_range.pl $index.retriever.scn $index.ltrTE.stg1 -x -strict 1`;
`cat $index.retriever.scn.extend |sort -fu > $index.retriever.scn.extend.unq; mv $index.retriever.scn.extend.unq $index.retriever.scn.extend`;
`perl $script_path/bin/call_seq_by_list.pl $index.retriever.scn.extend -C $genome > $index.retriever.scn.extend.fa`; #full TE sequence with 50bp-extended on each side

# use TE HMM to classify candidates
`perl $script_path/bin/Six-frame_translate.pl $index.retriever.scn.extend.fa > $index.retriever.scn.extend.fa.aa`; ##six-frame translate candidate sequences
`${hmmer}hmmsearch --tblout $index.retriever.scn.extend.fa.aa.tbl --notextw --cpu $threads -E 0.05 --domE 0.05 --noali $TEhmm $index.retriever.scn.extend.fa.aa > $index.retriever.scn.extend.fa.aa.scn`;
`touch $index.retriever.scn.extend.fa.aa.tbl` unless -s "$index.retriever.scn.extend.fa.aa.tbl";
`perl $script_path/bin/annotate_TE.pl $index.retriever.scn.extend.fa.aa.tbl > $index.retriever.scn.extend.fa.aa.anno`;

# use more TE HMM from TEsorter to classify candidates
`$TEsorter $index.retriever.scn.extend.fa --disable-pass2 -p $threads 2>/dev/null`;
`touch $index.retriever.scn.extend.fa.rexdb.cls.tsv` unless -s "$index.retriever.scn.extend.fa.rexdb.cls.tsv";
`awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$6"\\t"\$7}' $index.retriever.scn.extend.fa.rexdb.cls.tsv | perl -nle 's/([0-9]+\\.\\.[0-9]+)_(.*:[0-9]+\\.\\.[0-9]+)/\$1\\|\$2/; print \$_' >> $index.retriever.scn.extend.fa.aa.anno 2>/dev/null`;

# identify intact LTR retrotransposons
`perl $script_path/bin/LTR.identifier.pl $index -list $index.retriever.scn -seq $index.retriever.scn.extend.fa -anno $index.retriever.scn.extend.fa.aa.anno -flanksim $flanksim -flankmiss $flankmiss -flankaln $flankaln -minlen $minlen $tsdaln -u $miu -threads $threads -blastplus $blastplus -motif @motif > $index.defalse`;
`perl -nle 'next unless /pass/i; print \$_ unless /notLTR/i or /mixture/i' $index.defalse > $index.ltrTE.pass.list`;
my $count=0;
$count=`grep -v -c \'#\' $index.ltrTE.pass.list`;
$count=~s/\s+//g;
return $count;
}