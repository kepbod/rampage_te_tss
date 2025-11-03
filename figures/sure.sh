for i in {29..36};do gzip -d -c SRR40035${i}.fastq.gz|cutadapt -g GCTAGCTAACTATAACGGTCCTAAGGTAGCGAA --info-file=SRR40035${i}_info.txt - > /dev/null;done &
for i in {29..36};do cut -f 5 SRR40035${i}_info.txt|sort -S 50% --parallel=10|uniq -c|perl -alne '$l=length($F[1]);print "$F[0]\t$F[1]" if $l==20 and $F[1]!~/N/'|gzip -c >! SRR40035${i}_trimmed_table.txt.gz;done &

data_processing/code/merge-iPCR-cDNA-plDNA.bash -s samples.txt -o sure_count.txt.gz -i iPCR/iPCR_merge.txt.gz -l merge_count.log -n 20 &

zcat sure_count.txt.gz|perl -alne 'next if $.==1;$a[$_]+=$F[$_] for 5..$#F;END{print join "\t",@a[5..$#a]}' >! total_read_count.txt &
zcat sure_count.txt.gz|pl -alne 'next if $F[3] ne "+";$n=sum(@F[7..12]);$"="\t";print "@F[0..2]" for 1..$n'|bedtools genomecov -bg -i stdin -g $zlab/genome/hg38/chrom.sizes|sort --parallel=10 -k1,1 -k2,2n >! sure_cDNA_count_plus.bg &
zcat sure_count.txt.gz|pl -alne 'next if $F[3] ne "-";$n=sum(@F[7..12]);$"="\t";print "@F[0..2]" for 1..$n'|bedtools genomecov -bg -i stdin -g $zlab/genome/hg38/chrom.sizes|sort --parallel=10 -k1,1 -k2,2n >! sure_cDNA_count_minus.bg &
zcat sure_count.txt.gz|pl -alne 'next if $F[3] ne "+";$n=sum(@F[5,6])+1;$"="\t";print "@F[0..2]" for 1..$n'|bedtools genomecov -bg -i stdin -g $zlab/genome/hg38/chrom.sizes|sort --parallel=10 -k1,1 -k2,2n >! sure_plDNA_count_plus.bg &
zcat sure_count.txt.gz|pl -alne 'next if $F[3] ne "-";$n=sum(@F[5,6])+1;$"="\t";print "@F[0..2]" for 1..$n'|bedtools genomecov -bg -i stdin -g $zlab/genome/hg38/chrom.sizes|sort --parallel=10 -k1,1 -k2,2n >! sure_plDNA_count_minus.bg &

bigwigCompare --bigwig1 sure_cDNA_count_plus.bw --bigwig2 sure_plDNA_count_plus.bw --scaleFactors 34378305:69297597 --pseudocount 0 0 --skipZeroOverZero --binSize 1 -p 20 --outFileName sure_score_plus.bw --outFileFormat bigwig --skipNonCoveredRegions &
bigwigCompare --bigwig1 sure_cDNA_count_minus.bw --bigwig2 sure_plDNA_count_minus.bw --scaleFactors 34378305:69297597 --pseudocount 0 0 --skipZeroOverZero --binSize 1 -p 20 --outFileName sure_score_minus.bw --outFileFormat bigwig --skipNonCoveredRegions &

bigwigCompare --bigwig1 sure_cDNA_count_plus.bw --bigwig2 sure_plDNA_count_plus.bw --scaleFactors 34378305:69297597 --pseudocount 0 0 --skipZeroOverZero --binSize 1 -p 20 --outFileName sure_ratio_plus.bw --outFileFormat bigwig --skipNonCoveredRegions --operation ratio &
bigwigCompare --bigwig1 sure_cDNA_count_minus.bw --bigwig2 sure_plDNA_count_minus.bw --scaleFactors 34378305:69297597 --pseudocount 0 0 --skipZeroOverZero --binSize 1 -p 20 --outFileName sure_ratio_minus.bw --outFileFormat bigwig --skipNonCoveredRegions  --operation ratio &

perl -alne '$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[-1] ne "other" and $F[5] eq "+"' K562_rampage_TE_gene_all_labeled.txt >! K562_TE_TSS_plus.bed
perl -alne '$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[-1] eq "other" and $F[5] eq "+"' K562_rampage_TE_gene_all_labeled.txt|shuf|head -500 >! K562_other_TE_TSS_plus.bed
perl -alne '$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[5] eq "+"' K562_rampage_TE_gene_high.txt >! K562_high_TE_TSS_plus.bed
perl -alne '$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[-1] ne "other" and $F[5] eq "-"' K562_rampage_TE_gene_all_labeled.txt >! K562_TE_TSS_minus.bed
perl -alne '$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[-1] eq "other" and $F[5] eq "-"' K562_rampage_TE_gene_all_labeled.txt|shuf|head -500 >! K562_other_TE_TSS_minus.bed
perl -alne '$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[5] eq "-"' K562_rampage_TE_gene_high.txt >! K562_high_TE_TSS_minus.bed
perl -alne '$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[-3] eq "known" and abs($F[-2])<=50 and $F[5] eq "+"' K562_assigned_peaks.txt|shuf|head -500 >! K562_mRNA_TSS_plus.bed
perl -alne '$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[-3] eq "known" and abs($F[-2])<=50 and $F[5] eq "-"' K562_assigned_peaks.txt|shuf|head -500 >! K562_mRNA_TSS_minus.bed

computeMatrix reference-point -S sure_ratio_plus.bw -R K562_TE_TSS_plus.bed -a 5000 -b 5000 --referencePoint center -out K562_TE_TSS_plus_sure.gz &
computeMatrix reference-point -S sure_ratio_plus.bw -R K562_high_TE_TSS_plus.bed -a 5000 -b 5000 --referencePoint center -out K562_high_TE_TSS_plus_TE_TSS_plus_sure.gz &
computeMatrix reference-point -S sure_ratio_minus.bw -R K562_TE_TSS_minus.bed -a 5000 -b 5000 --referencePoint center -out K562_TE_TSS_minus_sure.gz &
computeMatrix reference-point -S sure_ratio_minus.bw -R K562_high_TE_TSS_minus.bed -a 5000 -b 5000 --referencePoint center -out K562_high_TE_TSS_minus_TE_TSS_minus_sure.gz &
computeMatrix reference-point -S sure_ratio_plus.bw -R K562_other_TE_TSS_plus.bed -a 5000 -b 5000 --referencePoint center -out K562_other_TE_TSS_plus_TE_TSS_plus_sure.gz &
computeMatrix reference-point -S sure_ratio_minus.bw -R K562_other_TE_TSS_minus.bed -a 5000 -b 5000 --referencePoint center -out K562_other_TE_TSS_minus_TE_TSS_minus_sure.gz &
./combine_sure_score.py K562_TE_TSS_plus_sure.gz K562_TE_TSS_minus_sure.gz all >! K562_TE_TSS_sure_plot.txt
./combine_sure_score.py K562_high_TE_TSS_plus_TE_TSS_plus_sure.gz K562_high_TE_TSS_minus_TE_TSS_minus_sure.gz high >> K562_TE_TSS_sure_plot.txt
./combine_sure_score.py K562_other_TE_TSS_plus_TE_TSS_plus_sure.gz K562_other_TE_TSS_minus_TE_TSS_minus_sure.gz other >> K562_TE_TSS_sure_plot.txt

computeMatrix reference-point -S sure_ratio_plus.bw -R K562_mRNA_TSS_plus.bed -a 5000 -b 5000 --referencePoint center -out K562_mRNA_TSS_plus_TE_TSS_plus_sure.gz &
computeMatrix reference-point -S sure_ratio_minus.bw -R K562_mRNA_TSS_minus.bed -a 5000 -b 5000 --referencePoint center -out K562_mRNA_TSS_minus_TE_TSS_minus_sure.gz &
./combine_sure_score.py K562_TE_TSS_plus_sure.gz K562_TE_TSS_minus_sure.gz all >! K562_TE_TSS_all_sure_plot.txt
./combine_sure_score.py K562_high_TE_TSS_plus_TE_TSS_plus_sure.gz K562_high_TE_TSS_minus_TE_TSS_minus_sure.gz high >> K562_TE_TSS_all_sure_plot.txt
./combine_sure_score.py K562_other_TE_TSS_plus_TE_TSS_plus_sure.gz K562_other_TE_TSS_minus_TE_TSS_minus_sure.gz other >> K562_TE_TSS_all_sure_plot.txt
./combine_sure_score.py K562_mRNA_TSS_plus_TE_TSS_plus_sure.gz K562_mRNA_TSS_minus_TE_TSS_minus_sure.gz mRNA >> K562_TE_TSS_all_sure_plot.txt

perl -alne '$a=$F[3];$b=$F[3]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[4]" if $F[4] eq "+"' TE_TSS_cluster.txt >! TE_TSS.plus.bed
perl -alne '$a=$F[3];$b=$F[3]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[4]" if $F[4] eq "-"' TE_TSS_cluster.txt >! TE_TSS.minus.bed
perl -alne '$a=$F[3];$b=$F[3]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[4]" if $F[4] eq "+"' high_TE_TSS_cluster.txt >! high_TE_TSS.plus.bed
perl -alne '$a=$F[3];$b=$F[3]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[4]" if $F[4] eq "-"' high_TE_TSS_cluster.txt >! high_TE_TSS.minus.bed
perl -alne '$a=$F[3];$b=$F[3]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[4]" if $F[4] eq "+"' other_TE_TSS_cluster.txt|shuf|head -10000 >! other_TE_TSS.plus.bed
perl -alne '$a=$F[3];$b=$F[3]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[4]" if $F[4] eq "-"' other_TE_TSS_cluster.txt|shuf|head -10000 >! other_TE_TSS.minus.bed
computeMatrix reference-point -S sure_ratio_plus.bw -R TE_TSS.plus.bed -a 5000 -b 5000 --referencePoint center -out TE_TSS_plus_sure.gz &
computeMatrix reference-point -S sure_ratio_plus.bw -R high_TE_TSS.plus.bed -a 5000 -b 5000 --referencePoint center -out high_TE_TSS_plus_TE_TSS_plus_sure.gz &
computeMatrix reference-point -S sure_ratio_minus.bw -R TE_TSS.minus.bed -a 5000 -b 5000 --referencePoint center -out TE_TSS_minus_sure.gz &
computeMatrix reference-point -S sure_ratio_minus.bw -R high_TE_TSS.minus.bed -a 5000 -b 5000 --referencePoint center -out high_TE_TSS_minus_TE_TSS_minus_sure.gz &
computeMatrix reference-point -S sure_ratio_plus.bw -R other_TE_TSS.plus.bed -a 5000 -b 5000 --referencePoint center -out other_TE_TSS_plus_TE_TSS_plus_sure.gz &
computeMatrix reference-point -S sure_ratio_minus.bw -R other_TE_TSS.minus.bed -a 5000 -b 5000 --referencePoint center -out other_TE_TSS_minus_TE_TSS_minus_sure.gz &
./combine_sure_score.py TE_TSS_plus_sure.gz TE_TSS_minus_sure.gz all >! TE_TSS_sure_plot.txt
./combine_sure_score.py high_TE_TSS_plus_TE_TSS_plus_sure.gz high_TE_TSS_minus_TE_TSS_minus_sure.gz high >> TE_TSS_sure_plot.txt
./combine_sure_score.py other_TE_TSS_plus_TE_TSS_plus_sure.gz other_TE_TSS_minus_TE_TSS_minus_sure.gz other >> TE_TSS_sure_plot.txt

zcat sure_count.txt.gz|pl -alne '$l=$F[2]-$F[1];next if $.==1 or $l>1500;$sum=sum(@F[7..$#F]);$"="\t";print "@F[0..2]\tsure\t$sum\t$F[3]" if $sum'|bgzip >! sure_cDNA_count.bed.gz
tabix -p bed sure_cDNA_count.bed.gz

perl -alne '$"="\t";if($ARGV=~/high/){$a{"@F[0..5]"}++}else{$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[-1] ne "other" and not $a{"@F[0..5]"} and $F[5] eq "+"}' K562_rampage_TE_gene_high.txt K562_rampage_TE_gene_all_labeled.txt >! K562_rampage_TE_TSS_plus.bed
perl -alne '$"="\t";if($ARGV=~/high/){$a{"@F[0..5]"}++}else{$a=$F[6];$b=$F[6]+1;print "$F[0]\t$a\t$b\tTSS\t0\t$F[5]" if $F[-1] ne "other" and not $a{"@F[0..5]"} and $F[5] eq "-"}' K562_rampage_TE_gene_high.txt K562_rampage_TE_gene_all_labeled.txt >! K562_rampage_TE_TSS_minus.bed
computeMatrix reference-point -S sure_ratio_plus.bw -R K562_rampage_TE_TSS_plus.bed -a 5000 -b 5000 --referencePoint center -out K562_rampage_TE_TSS_plus_sure.gz &
computeMatrix reference-point -S sure_ratio_minus.bw -R K562_rampage_TE_TSS_minus.bed -a 5000 -b 5000 --referencePoint center -out K562_rampage_TE_TSS_minus_sure.gz &
./combine_sure_score.py K562_rampage_TE_TSS_plus_sure.gz K562_rampage_TE_TSS_minus_sure.gz rampage >! K562_TE_TSS_all_sure_plot_new.txt
./combine_sure_score.py K562_high_TE_TSS_plus_TE_TSS_plus_sure.gz K562_high_TE_TSS_minus_TE_TSS_minus_sure.gz high >> K562_TE_TSS_all_sure_plot_new.txt
./combine_sure_score.py K562_other_TE_TSS_plus_TE_TSS_plus_sure.gz K562_other_TE_TSS_minus_TE_TSS_minus_sure.gz other >> K562_TE_TSS_all_sure_plot_new.txt
./combine_sure_score.py K562_mRNA_TSS_plus_TE_TSS_plus_sure.gz K562_mRNA_TSS_minus_TE_TSS_minus_sure.gz mRNA >> K562_TE_TSS_all_sure_plot_new.txt
