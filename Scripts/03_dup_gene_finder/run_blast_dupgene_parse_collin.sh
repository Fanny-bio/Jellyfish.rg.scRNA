##Prepare the input files (longest-transcript protein files and gff files) according to https://github.com/qiao-xin/DupGen_finder, and compare each two species as below

num_threads="10"
id1="Species1" ##Species1
id2="Species2"  ##Species2
db_name="$(echo ""$id1"_"$id2"")"

cat $id1.fa $id2.fa > ${id1}_${id2}.fa 
cat $id1.gff $id2.gff > ${id1}_${id2}.gff
prot="$(echo ""$id1"_"$id2".fa")"

diamond makedb --in $prot -d $db_name

diamond blastp -d $db_name --query $prot --evalue  1e-10 --outfmt 6 --threads $num_threads -o $db_name.blast.out.txt


grep ""$id1"|.*"$id1"|" $db_name.blast.out.txt > ${id1}.blast
grep ""$id2"|.*"$id2"|" $db_name.blast.out.txt > ${id2}.blast

grep ""$id1"|.*"$id2"|\|"$id2"|.*"$id1"|" $db_name.blast.out.txt > $db_name.blast

ln -s ${id1}_${id2}.gff ${id2}_${id1}.gff 
ln -s $db_name.blast ${id2}_${id1}.blast


###Dupgenefinder
~/bin/DupGen_finder/DupGen_finder-unique.pl -i ./ -t $id1 -c $id2 -o DupGenFinder_${id1}_${id2}_result
~/bin/DupGen_finder/DupGen_finder-unique.pl -i ./ -t $id2 -c $id1 -o DupGenFinder_${id2}_${id1}_result

##Parse collinearity pairs
cat DupGenFinder_${id1}_${id2}_result/${id1}_${id2}.collinearity | grep -v '#' | sed 's/.*:\s\{1,\}//' | sed 's/\s\{1,\}[0-9e-]\{1,\}//' > ${id1}_${id2}.collinearity.pairs.txt

cat DupGenFinder_${id1}_${id2}_result/${id1}.collinearity | grep -v '#' | sed 's/.*:\s\{1,\}//' | sed 's/\s\{1,\}[0-9e-]\{1,\}//' > ${id1}.collinearity.pairs.txt

cat DupGenFinder_${id2}_${id1}_result/${id2}.collinearity | grep -v '#' | sed 's/.*:\s\{1,\}//' | sed 's/\s\{1,\}[0-9e-]\{1,\}//' > ${id2}.collinearity.pairs.txt

##done
