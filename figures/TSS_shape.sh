perl -alne 'if($#F!=7){$ARGV=~s/shrink_region\/|_rampage_TE_gene_shrink.txt//g;$a{"$ARGV|$F[1]|$F[2]|$F[6]|$F[5]"}=$F[-1]}else{for(split /,/,$F[-1]){print "high\t$a{$_}"}}' shrink_region/*_rampage_TE_gene_shrink.txt high_TE_TSS_cluster.txt >! len_dis.txt
perl -alne 'if($#F!=7){$ARGV=~s/shrink_region\/|_rampage_TE_gene_shrink.txt//g;$a{"$ARGV|$F[1]|$F[2]|$F[6]|$F[5]"}=$F[-1]}else{for(split /,/,$F[-1]){print "all\t$a{$_}"}}' shrink_region/*_rampage_TE_gene_shrink.txt TE_TSS_cluster.txt >> len_dis.txt
perl -alne '$l=$F[7]-$F[6];print "gene\t$l"' hg38-rPeaks.bed >> len_dis.txt

## TATA
perl -alne '$"="|";$id="@F[0..3]";$a=$F[3]-50;$b=$F[3]+50;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' high_TE_TSS_cluster.txt >! high_TE_TSS_region.bed
perl -alne '$"="|";$id="@F[0..3]";$a=$F[3]-50;$b=$F[3]+50;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' TE_TSS_cluster.txt >! TE_TSS_region.bed
bedtools getfasta -fi $zlab/genome/hg38/hg38_all.fa -bed high_TE_TSS_region.bed -fo high_TE_TSS_region.fa -s -name
bedtools getfasta -fi $zlab/genome/hg38/hg38_all.fa -bed TE_TSS_region.bed -fo TE_TSS_region.fa -s -name
fimo --norc --thresh 1e-2 --oc high_TE_TSS_region_TATA MA0108.2.meme high_TE_TSS_region.fa
fimo --norc --thresh 1e-2 --oc TE_TSS_region_TATA MA0108.2.meme TE_TSS_region.fa
perl -alne '$id=(split /\(/,$F[2])[0];if($F[3]<=25 and $F[4]>=15){if($a{$id}){$a{$id}=$a{$id}<$F[7]?$a{$id}:$F[7]}else{$a{$id}=$F[7]}};END{print "$_\t$a{$_}" for keys %a}' TE_TSS_region_TATA/fimo.tsv >! TE_TSS_TATA.txt
perl -alne '$id=(split /\(/,$F[2])[0];if($F[3]<=25 and $F[4]>=15){if($a{$id}){$a{$id}=$a{$id}<$F[7]?$a{$id}:$F[7]}else{$a{$id}=$F[7]}};END{print "$_\t$a{$_}" for keys %a}' high_TE_TSS_region_TATA/fimo.tsv >! high_TE_TSS_TATA.txt

perl -alne '$a=$F[3]-50;$b=$F[3]+50; print "$F[0]\t$a\t$b\t$F[1]\t0\t$F[2]"' all_peak_summit.txt >! all_TSS_region.bed
bedtools getfasta -fi $zlab/genome/hg38/hg38_all.fa -bed all_TSS_region.bed -fo all_TSS_region.fa -s -name
fimo --norc --thresh 1e-2 --oc all_TSS_region_TATA MA0108.2.meme all_TSS_region.fa
perl -alne '$id=(split /\(/,$F[2])[0];if($F[3]<=25 and $F[4]>=15){if($a{$id}){$a{$id}=$a{$id}<$F[7]?$a{$id}:$F[7]}else{$a{$id}=$F[7]}};END{print "$_\t$a{$_}" for keys %a}' all_TSS_region_TATA/fimo.tsv >! all_TSS_TATA.txt

## CpG
perl -alne '$"="\t";print "@F[1..3]\tcpg\t0\t+\t$F[-1]"' cpgIslandExtUnmasked.txt >! cpg.bed
perl -alne '$"="|";$id="@F[0..3]";$a=$F[3];$b=$F[3]+1;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' high_TE_TSS_cluster.txt >! high_TE_TSS_summit.bed
perl -alne '$"="|";$id="@F[0..3]";$a=$F[3];$b=$F[3]+1;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' TE_TSS_cluster.txt >! TE_TSS_summit.bed
bedtools intersect -a high_TE_TSS_summit.bed -b cpg.bed -wa -wb|perl -alne 'print "@F[3]\t$F[-1]"' >! high_TE_TSS_cpg.txt
bedtools intersect -a TE_TSS_summit.bed -b cpg.bed -wa -wb|perl -alne 'print "@F[3]\t$F[-1]"' >! TE_TSS_cpg.txt

perl -alne '$a=$F[3];$b=$F[3]+1;print "$F[0]\t$a\t$b\t$F[1]\t$F[2]"' all_peak_summit.txt >! all_TSS_summit.bed
bedtools intersect -a all_TSS_summit.bed -b cpg.bed -wa -wb|perl -alne 'print "@F[3]\t$F[-1]"' >! all_TSS_cpg.txt

## summary
shuf all_TSS_region.bed|head -n 30000 >! random_TSS_region.bed
perl -alne 'if($#F==1){$ARGV=~s/TE_TSS_|.txt//g;$a{$F[0]}{$ARGV}=$F[1]}else{$ta=$a{$F[3]}{TATA}?$a{$F[3]}{TATA}:"0.01";$cg=$a{$F[3]}{cpg}?$a{$F[3]}{cpg}:"0.5";print "TE_TSS\t$F[3]\t$ta\t$cg"}' TE_TSS_TATA.txt TE_TSS_cpg.txt TE_TSS_region.bed|sort -k 4,4g >! TATA_cpg_summary.txt
perl -alne 'if($#F==1){$ARGV=~s/high_TE_TSS_|.txt//g;$a{$F[0]}{$ARGV}=$F[1]}else{$ta=$a{$F[3]}{TATA}?$a{$F[3]}{TATA}:"0.01";$cg=$a{$F[3]}{cpg}?$a{$F[3]}{cpg}:"0.5";print "high_TE_TSS\t$F[3]\t$ta\t$cg"}' high_TE_TSS_TATA.txt high_TE_TSS_cpg.txt high_TE_TSS_region.bed|sort -k 4,4g >> TATA_cpg_summary.txt
perl -alne 'if($#F==1){$ARGV=~s/all_TSS_|.txt//g;$a{$F[0]}{$ARGV}=$F[1]}else{$ta=$a{$F[3]}{TATA}?$a{$F[3]}{TATA}:"0.01";$cg=$a{$F[3]}{cpg}?$a{$F[3]}{cpg}:"0.5";print "random_TSS\t$F[3]\t$ta\t$cg"}' all_TSS_TATA.txt all_TSS_cpg.txt random_TSS_region.bed|sort -k 4,4g >> TATA_cpg_summary.txt
perl -alne '$F[2]=$F[2]==0.01?"NA":-log($F[2])/log(10);$F[3]=$F[3]==0.5?"NA":$F[3];$"="\t";print "$.\t@F"' TATA_cpg_summary.txt >! TATA_cpg_summary_labeled.txt
ggplot(a, aes(factor('x'),V1,fill=V4)) + geom_tile() + theme_classic() + facet_wrap(~V2, nrow=3) + scale_fill_continuous(high='#800026', low='#fc4e2a', na.value='#f0f0f0') + coord_flip()
ggplot(a, aes(factor('x'),V1,fill=V5)) + geom_tile() + theme_classic() + facet_wrap(~V2, nrow=3) + scale_fill_continuous(high='#045a8d', low='#74a9cf', na.value='#f0f0f0') + coord_flip()

## NP or BP
perl -alne 'if($#F!=7){$ARGV=~s/shrink_region\/|_rampage_TE_gene_shrink.txt//g;$a{"$ARGV|$F[1]|$F[2]|$F[6]|$F[5]"}=$F[-1]}else{$t=0;$n=0;for(split /,/,$F[-1]){$t++;$n++ if $a{$_}<=10};$"="|";$id="@F[0..3]";$t=$n/$t>=0.5?1:0;print "$id\t$t"}' shrink_region/*_rampage_TE_gene_shrink.txt high_TE_TSS_cluster.txt >! high_TE_TSS_type.txt
perl -alne 'if($#F!=7){$ARGV=~s/shrink_region\/|_rampage_TE_gene_shrink.txt//g;$a{"$ARGV|$F[1]|$F[2]|$F[6]|$F[5]"}=$F[-1]}else{$t=0;$n=0;for(split /,/,$F[-1]){$t++;$n++ if $a{$_}<=10};$"="|";$id="@F[0..3]";$t=$n/$t>=0.5?1:0;print "$id\t$t"}' shrink_region/*_rampage_TE_gene_shrink.txt TE_TSS_cluster.txt >! TE_TSS_type.txt
perl -alne '$l=$F[7]-$F[6];$t=$l<=10?1:0;print "$F[3]\t$t"' hg38-rPeaks.bed >! random_TSS_type.txt
perl -alne 'if($#F==1){$ARGV=~s/_type.txt//;$a{"$ARGV\t$F[0]"}=$F[1]}else{print "$_\t".$a{"$F[1]\t$F[2]"}}' high_TE_TSS_type.txt TE_TSS_type.txt random_TSS_type.txt TATA_cpg_summary_labeled.txt >! TATA_cpg_summary_TSS_type.txt
pl -alne 'if($#F!=7){$ARGV=~s/shrink_region\/|_rampage_TE_gene_shrink.txt//g;$a{"$ARGV|$F[1]|$F[2]|$F[6]|$F[5]"}=$F[-1]}else{@x=();for(split /,/,$F[-1]){push @x,$a{$_}};$n=$#x+1;next if $n<2;$"="|";$id="@F[0..3]";$m=mean(@x);$s=sd(@x);$l=$m<=10?0:1;print "$id\t$n\t$m\t$s\t$l\thigh"}' shrink_region/*_rampage_TE_gene_shrink.txt high_TE_TSS_cluster.txt >! TE_TSS_type_sd.txt
pl -alne 'if($#F!=7){$ARGV=~s/shrink_region\/|_rampage_TE_gene_shrink.txt//g;$a{"$ARGV|$F[1]|$F[2]|$F[6]|$F[5]"}=$F[-1]}else{@x=();for(split /,/,$F[-1]){push @x,$a{$_}};$n=$#x+1;next if $n<2;$"="|";$id="@F[0..3]";$m=mean(@x);$s=sd(@x);$l=$m<=10?0:1;print "$id\t$n\t$m\t$s\t$l\tall"}' shrink_region/*_rampage_TE_gene_shrink.txt TE_TSS_cluster.txt >> TE_TSS_type_sd.txt


## for poliii_alu
perl -alne 'if($#F!=7){$ARGV=~s/shrink_region\/|_rampage_TE_gene_shrink.txt//g;$a{"$ARGV|$F[1]|$F[2]|$F[6]|$F[5]"}=$F[-1]}else{$t=0;$n=0;for(split /,/,$F[-1]){$t++;$n++ if $a{$_}<=10};if($F[4] eq "+"){$a=$F[3]-34;$b=$F[3]+6}else{$a=$F[3]-5;$b=$F[3]+35};print "$F[0]\t$a\t$b\tsummit\t0\t$F[4]\t$t\t$n" if $n/$t>=0.5}' shrink_region/*_rampage_TE_gene_shrink.txt high_TE_TSS_cluster.txt >! high_TE_TSS_NP_summit.bed
perl -alne 'if($#F!=7){$ARGV=~s/shrink_region\/|_rampage_TE_gene_shrink.txt//g;$a{"$ARGV|$F[1]|$F[2]|$F[6]|$F[5]"}=$F[-1]}else{$t=0;$n=0;for(split /,/,$F[-1]){$t++;$n++ if $a{$_}<=10};if($F[4] eq "+"){$a=$F[3]-34;$b=$F[3]+6}else{$a=$F[3]-5;$b=$F[3]+35};print "$F[0]\t$a\t$b\tsummit\t0\t$F[4]\t$t\t$n" if $n/$t>=0.5}' shrink_region/*_rampage_TE_gene_shrink.txt TE_TSS_cluster.txt >! TE_TSS_NP_summit.bed

## CG%
perl -alne '$a=$F[1]-50;$b=$F[1]+50;print "$F[0]\t$a\t$b\t$F[3]\t0\t$F[5]"' high_TE_TSS_summit.bed >! high_TE_TSS_summit_region.bed
perl -alne '$a=$F[1]-50;$b=$F[1]+50;print "$F[0]\t$a\t$b\t$F[3]\t0\t$F[5]"' TE_TSS_summit.bed >! TE_TSS_summit_region.bed
perl -alne '$a=$F[1]-50;$b=$F[1]+50;print "$F[0]\t$a\t$b\t$F[3]\t0\t$F[5]"' all_TSS_summit.bed >! all_TSS_summit_region.bed
bedtools getfasta -fi $zlab/genome/hg38/hg38_all.fa -bed high_TE_TSS_summit_region.bed -fo high_TE_TSS_summit_region.fa -s -name
bedtools getfasta -fi $zlab/genome/hg38/hg38_all.fa -bed TE_TSS_summit_region.bed -fo TE_TSS_summit_region.fa -s -name
bedtools getfasta -fi $zlab/genome/hg38/hg38_all.fa -bed all_TSS_summit_region.bed -fo all_TSS_summit_region.fa -s -name
perl -alne 'next if /^>/;$n=0;while(/C|G/ig){$n++};print "high_TE_TSS\t$n"' high_TE_TSS_summit_region.fa >! GC.txt
perl -alne 'next if /^>/;$n=0;while(/C|G/ig){$n++};print "TE_TSS\t$n"' TE_TSS_summit_region.fa >> GC.txt
perl -alne 'next if /^>/;$n=0;while(/C|G/ig){$n++};print "all_TSS\t$n"' all_TSS_summit_region.fa >> GC.txt
perl -alne 'next if /^>/;$cg=0;$c=0;$g=0;while(/C/ig){$c++};while(/G/ig){$g++};while(/CG/ig){$cg++};$cg/=100;$c/=100;$g/=100;$r=$cg/(($c+$g)/2)**2;print "high_TE_TSS\t$r"' high_TE_TSS_summit_region.fa >! GC_content.txt
perl -alne 'next if /^>/;$cg=0;$c=0;$g=0;while(/C/ig){$c++};while(/G/ig){$g++};while(/CG/ig){$cg++};$cg/=100;$c/=100;$g/=100;$r=$cg/(($c+$g)/2)**2;print "TE_TSS\t$r"' TE_TSS_summit_region.fa >> GC_content.txt
perl -alne 'next if /^>/;$cg=0;$c=0;$g=0;while(/C/ig){$c++};while(/G/ig){$g++};while(/CG/ig){$cg++};$cg/=100;$c/=100;$g/=100;$r=$cg/(($c+$g)/2)**2;print "all_TSS\t$r"' all_TSS_summit_region.fa >> GC_content.txt
