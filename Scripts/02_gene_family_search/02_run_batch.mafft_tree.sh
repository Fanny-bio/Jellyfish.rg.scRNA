

cat pfam_CDD_name.list | while read pfam CDD name; do 

PF_pref="$( echo "$pfam" | sed 's/pfam/PF/' )" &&  hmm="$(grep "$PF_pref" ~/bin/hmmer-3.3.1/Pfam_seed/Pfam-A.hmm | sed 's/ACC.*PF/PF/')"  && seqtk subseq ${name}/$hmm.scan.dedup.protein.fa ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.txt > ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.protein.fa && mafft --genafpair --maxiterate 1000 --ep 0 --quiet  ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.protein.fa > ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.protein.align.fa && echo `date` "==> Done mafft: "$hmm" <==" && ~/bin/FastTree -gamma -quote  ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.protein.align.fa >  ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.protein.align.tre && sed "s/'//g" ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.protein.align.tre > ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.protein.align.dequote.tre && Rscript ladderize_tree.R -i ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.protein.align.dequote.tre -o ${name}/$hmm.scan.dedup.protein.rpsblast.confirm.protein.align.dequote.vis & 

done



