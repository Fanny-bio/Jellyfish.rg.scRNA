cat pfam_CDD_name.list | while read pfam CDD name; do 

##protein file contain the longest-transcripts of 15 species 
prot="15species.fa"

PF_pref="$( echo "$pfam" | sed 's/pfam/PF/' )"
hmm="$(grep "$PF_pref" ~/bin/hmmer-3.3.1/Pfam_seed/Pfam-A.hmm | sed 's/ACC.*PF/PF/')" 

mkdir $name
cd $name

hmmfetch ~/bin/hmmer-3.3.1/Pfam_seed/Pfam-A.hmm $hmm > $hmm.hmm

hmmpress $hmm.hmm

hmmscan --domtblout $hmm.scan.out.txt -E 1e-5 $hmm.hmm $prot


awk '!/#/ {print $4}' $hmm.scan.out.txt | awk '!a[$0]++' > $hmm.scan.dedup.txt 

seqtk subseq  $prot $hmm.scan.dedup.txt  > $hmm.scan.dedup.protein.fa 

sed 's/|.*//' $hmm.scan.dedup.txt | awk '{A[$1]++}END{for(i in A)print i,A[i]}' > $hmm.scan.dedup.count

echo `date` "Done hmmscan for pfam "$name" ("$hmm")"

~/bin/RpsbProc-x64-linux/rpsblast -query  $hmm.scan.dedup.protein.fa -out $hmm.scan.dedup.protein.rpsblast.out -evalue 0.01 -outfmt 6 -db /home/sean/bin/RpsbProc-x64-linux/db/Pfam -num_threads 4

awk '!a[$1$2]++' $hmm.scan.dedup.protein.rpsblast.out |  grep "CDD:"$CDD"" | cut -f1 > $hmm.scan.dedup.protein.rpsblast.confirm.txt

sed 's/|.*//' $hmm.scan.dedup.protein.rpsblast.confirm.txt | awk '{A[$1]++}END{for(i in A)print i,A[i]}' > $hmm.scan.dedup.protein.rpsblast.confirm.count

echo `date` "done rpsblast "$name" ("$hmm")"


cd ../


done
