### exclude biosamples with fewer than 10 million uniquely mappped reads
### exclude biosamples with fewer than 25% of non-redundant reads
perl -alne 'if($#F==1){$a{$F[0]}++}else{$l=(split /\//,$F[2])[3];s/ TE/ filtered_TE/;print if $a{$l}}' /data/zusers/moorej3/moorej.ghpcc.project/ENCODE/RAMPAGE/RAMPAGE-List-Filtered.txt ln.sh >! filter_ln.sh

### tss dis cutoff (<1kb)
perl -alne '$"="\t";print "TETSS\t$ARGV\t$F[-2]"' TE_genes/*rampage_TE_gene.txt >! tss_dis_cutoff.txt
perl -alne '$"="\t";print "other\t$ARGV\t$F[-3]" if $F[-1] eq "other" and $F[-3] ne "None"' TE_genes/*rampage_TE_gene_all_labeled.txt >> tss_dis_cutoff.txt

### peaks
perl -alne '$ARGV=~s/TE_genes\/|_rampage_TE_gene_all_labeled.txt//g;$g{$ARGV}++ if $F[-1] ne "other"; $r{$ARGV}++ if $F[-1] eq "other";END{for(keys %r){print "$_\tall\t$g{$_}\n$_\tother\t$r{$_}"}}' TE_genes/*rampage_TE_gene_all_labeled.txt >! peak_assign_count_tss_type.txt
#### Fig 1B
perl -alne '$ARGV=~s/TE_genes\/|_rampage_TE_gene.txt//g;$high{$ARGV}++ if abs($F[-2])<=1000;END{for(keys %high){print "$_\thigh\t$high{$_}"}}' TE_genes/*rampage_TE_gene.txt >> peak_assign_count_tss_type.txt
for i in TE_genes/*_rampage_TE_gene_all_labeled.txt;do perl -alne '$F[-1]=$F[-1]ne"other"?"all":$F[-1];print "$F[-1]\t$F[12]"' $i;done >! tss_type_exp.txt
for i in TE_genes/high_confidence/*_rampage_TE_gene_high.txt;do perl -alne 'print "high\t$F[12]"' $i;done >> tss_type_exp.txt
#### Fig S1E
perl -alne '$ARGV=~s/TE_genes\/high_confidence\/|_rampage_TE_gene_high.txt//g;$a{$ARGV}{$F[-3]}++;END{for(keys %a){$k=$a{$_}{known}?$a{$_}{known}:0;$n=$a{$_}{novel}?$a{$_}{novel}:0;print "$_\tknown\t$k\n$_\tnovel\t$n"}}' TE_genes/high_confidence/*rampage_TE_gene_high.txt >! peak_assign_count_RNA_tss.txt
for i in TE_genes/high_confidence/*_rampage_TE_gene_high.txt;do perl -alne '$ARGV=~s/TE_genes\/high_confidence\/|_rampage_TE_gene_high.txt//g;print "$ARGV\t$F[-3]\t$F[12]"' $i;done >! RNA_tss_type_exp.txt

### TEs
#perl -alne '$ARGV=~s/TE_genes\/|_rampage_TE_gene.txt//g;$t=(split /\|/,$F[19]) [0];$t=~s/\?//;$x=(split /\|/,$F[19])[1];$x=~s/\?//; $x=$x=~/^(Alu|MIR|L1|L2|ERV|SVA)/?$x:"other";$a{$ARGV}{"$t|$x"}++; END{for$i(keys %a){for("SINE|Alu", "SINE|MIR", "LINE|L1", "LINE|L2", "LTR|ERV1", "LTR|ERVK", "LTR|ERVL", "LTR|ERVL-MaLR", "Retroposon|SVA", "DNA|other"){$x=$a{$i}{$_}?$a{$i}{$_}:0;print "$i\t$_\t$x"}}}' TE_genes/*rampage_TE_gene.txt >! te_type.txt
perl -alne '$ARGV=~s/TE_genes\/|_rampage_TE_gene.txt//g;$t=(split /\|/,$F[19]) [0];$t=~s/\?//;$x=(split /\|/,$F[19])[1];$x=~s/\?//; $x=$x=~/^(Alu|MIR|L1|L2|ERV|SVA)/?$x:"other";$a{$ARGV}{"$t|$x"}++; END{for$i(keys %a){for("SINE|Alu", "SINE|MIR", "SINE|other", "LINE|L1", "LINE|L2", "LINE|other", "LTR|ERV1", "LTR|ERVK", "LTR|ERVL", "LTR|ERVL-MaLR", "LTR|other", "Retroposon|SVA", "DNA|other"){$x=$a{$i}{$_}?$a{$i}{$_}:0;print "$i\t$_\t$x"}}}' TE_genes/*rampage_TE_gene.txt >! te_type.txt
perl -alne '$ARGV=~s/TE_genes\/|_rampage_TE_gene.txt//g;$t=(split /\|/,$F[19])[0];$t=~s/\?//;$a{$ARGV}{$t}++;END{for$i(keys %a){for("SINE", "LINE", "LTR", "Retroposon", "DNA"){$x=$a{$i}{$_}?$a{$i}{$_}:0;print "$i\t$_\t$x"}}}' TE_genes/*rampage_TE_gene.txt >! te_class_type.txt
ggplot(a, aes(factor(V2, levels=c('SINE', 'LINE', 'LTR', 'Retroposon', 'DNA')), V3)) + geom_boxplot(color='#E99B4A', fill='#E99B4A') + theme_classic()
ggplot(a, aes(factor(V2, levels=c('SINE|Alu', 'SINE|MIR', 'SINE|other', 'LINE|L1', 'LINE|L2', 'LINE|other', 'LTR|ERVK', 'LTR|ERVL', 'LTR|ERVL-MaLR', 'LTR|ERV1', 'LTR|other')), V3)) + geom_boxplot(color='#E99B4A', fill='#E99B4A') + theme_classic()

cut -f 17-22 filtered_TE_peaks/*rampage_TE.txt|sort|uniq -c|perl -alne 'if($F[0]>=10){$a{10}++}else{$a{$F[0]}++};END{print "$_\t$a{$_}" for sort {$a<=>$b} keys %a}'
perl -alne '$ARGV=~s/filtered_TE_peaks\/|_rampage_TE.txt//g;$t=(split /\|/,$F[-4])[0];$t=~s/\?//;$t=$t eq "Retroposon"?"SINE":$t;$x=(split /\|/,$F[-4])[1];$a{$ARGV}{$x}++ if $t eq "DNA";END{for$i(keys %a){for(keys %{$a{$i}}){print "$i\t$_\t$a{$i}{$_}"}}}' filtered_TE_peaks/*rampage_TE.txt >! dna_type.txt
ggplot(a, aes(factor(V2, levels=c('RNA', 'RAMPAGE', 'other')), V3, color=factor(V2), fill=factor(V2))) + geom_boxplot(width=0.5) + theme_classic() + scale_y_log10() + scale_color_manual(values=c('#5E5E5E', '#0076BA', '#ED3938')) + scale_fill_manual(values=c('#5E5E5E', '#0076BA', '#ED3938'))

perl -alne '$"="\t";if($#F==8){$l{"@F[0..2]"}=$F[6];$r{"@F[0..2]"}=$F[7]}else{$t=(split /\|/,$F[3])[1];$t=~s/\?//;next if $t!~/^(Alu|MIR|L1|L2|ERV1|ERVK|ERVL|ERVL-MaLR)$/;$left=$l{"@F[0..2]"};$right=$r{"@F[0..2]"};if($F[5] eq "+"){$up=$F[1]-$left;$down=$right-$F[2]}else{$down=$F[1]-$left;$up=$right-$F[2]};print "$_\t$t\t$up\t$down"}' $zlab/genome/hg38/te_con.bed TE_genes/TE_TSS_cluster_te.bed >! TE_TSS_cluster_te_truncate_type.txt
perl -alne '$"="\t";if($#F==8){$l{"@F[0..2]"}=$F[6];$r{"@F[0..2]"}=$F[7]}else{$t=(split /\|/,$F[3])[1];$t=~s/\?//;next if $t!~/^(Alu|MIR|L1|L2|ERV1|ERVK|ERVL|ERVL-MaLR)$/;$left=$l{"@F[0..2]"};$right=$r{"@F[0..2]"};if($F[5] eq "+"){$up=$F[1]-$left;$down=$right-$F[2]}else{$down=$F[1]-$left;$up=$right-$F[2]};print "$_\t$t\t$up\t$down"}' $zlab/genome/hg38/te_con.bed TE_genes/high_confidence/high_TE_TSS_cluster_te.bed >! high_TE_TSS_cluster_te_truncate_type.txt
perl -alne '$"="\t";if($#F==8){$l{"@F[0..2]"}=$F[6];$r{"@F[0..2]"}=$F[7]}else{$t=(split /\|/,$F[3])[4];$t=~s/\?//;next if $t!~/^(Alu|MIR|L1|L2|ERV1|ERVK|ERVL|ERVL-MaLR)$/;$left=$l{"@F[0..2]"};$right=$r{"@F[0..2]"};if($F[5] eq "+"){$up=$F[1]-$left;$down=$right-$F[2]}else{$down=$F[1]-$left;$up=$right-$F[2]};print "$_\t$t\t$up\t$down"}' $zlab/genome/hg38/te_con.bed evolution/TSS_1k_te.bed >! TSS_1k_te_truncate_type.txt
perl -alne 'if($F[7]>10){$l=$F[8]>10?"both":"5truncated"}else{$l=$F[8]>10?"3truncated":"full"};$a{$F[6]}{$l}++;END{for$i(keys %a){for("full", "5truncated", "3truncated", "both"){$x=$a{$i}{$_}?$a{$i}{$_}:0;print "$i\t$_\t$x"}}}' TE_TSS_cluster_te_truncate_type.txt

### group TE tss
perl -alne 'if($#F==1){$a{$F[0]}=$F[1]}else{($g,$i)=split /\|/,$F[3];if($i!~/^ENST/){$g=$a{$g}?$a{$g}:$g};push @{$x{$g}},$_};END{for$i(keys %x){print "$_\t$i" for @{$x{$i}}}}' $zlab/genome/hg38/reftogencode_diff_name.txt $zlab/genome/hg38/ref.bed >! ref_group_labled.bed
# for i in TE_genes/high_confidence/*_rampage_TE_gene_high.txt;do j=${${i/TE_genes\/high_confidence\//}/_rampage_TE_gene_high.txt/}; perl -alne '$"="\t";if($ARGV=~/high/){$a{"@F[0..15]"}++;print "@F[23..34]\t$_\tTE_tss"}else{print "@F[16..27]\t@F[0..15]\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t@F[16..30]\tctrl_tss" if not $a{"@F[0..15]"}}' $i TE_genes/ctrl/${j}_rampage_TE_ctrl_tss.txt|bedtools intersect -a - -b ref_group_labled.bed -wao -s -split|perl -alne '$"="\t";$F[-2]=$F[-2]eq"."?"None":$F[-2];push @{$a{"@F[12..50]"}},$F[-2];END{for$k(keys %a){%s={};print "$k\t".join(",",grep(!$s{$_}++,@{$a{$k}}))}}' >! TE_TSS_group/${j}_TE_gene_group.txt;done
perl -alne '$l=$F[2]-$F[1];print if $l>1000' ref_group_labled.bed >!  ref_group_labled_long.bed
for i in TE_genes/high_confidence/*_rampage_TE_gene_high.txt;do j=${${i/TE_genes\/high_confidence\//}/_rampage_TE_gene_high.txt/}; perl -alne '$"="\t";if($ARGV!~/ctrl/){$a{"@F[0..15]"}++;print "@F[23..34]\t$_\tTE_tss"}else{print "@F[16..27]\t@F[0..15]\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t@F[16..30]\tctrl_tss" if not $a{"@F[0..15]"}}' $i TE_genes/ctrl/${j}_rampage_TE_ctrl_tss.txt|bedtools intersect -a - -b ref_group_labled_long.bed -wao -s -split|perl -alne '$"="\t";$F[-2]=$F[-2]eq"."?"None":$F[-2];push @{$a{"@F[12..50]"}},$F[-2];END{for$k(keys %a){%s={};print "$k\t".join(",",grep(!$s{$_}++,@{$a{$k}}))}}' >! TE_TSS_group/${j}_TE_gene_group.txt;done
for i in TE_genes/*_rampage_TE_gene.txt;do j=${${i/TE_genes\/}/_rampage_TE_gene.txt/}; perl -alne '$"="\t";if($ARGV!~/ctrl/){$a{"@F[0..15]"}++;print "@F[23..34]\t$_\tTE_tss"}else{print "@F[16..27]\t@F[0..15]\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t@F[16..30]\tctrl_tss" if not $a{"@F[0..15]"}}' $i TE_genes/ctrl/${j}_rampage_TE_ctrl_tss.txt|bedtools intersect -a - -b ref_group_labled_long.bed -wao -s -split|perl -alne '$"="\t";$F[-2]=$F[-2]eq"."?"None":$F[-2];push @{$a{"@F[12..50]"}},$F[-2];END{for$k(keys %a){%s={};print "$k\t".join(",",grep(!$s{$_}++,@{$a{$k}}))}}' >! all_TE_TSS_group/${j}_TE_gene_group.txt;done

### Fig 1B (from TE_TSS_cluster.txt in TE_genes and high_TE_TSS_cluster.txt in TE_genes/high_confidence
join TE_TSS_count.txt TE_TSS_high_count.txt >! TE_TSS_count_summary.txt
