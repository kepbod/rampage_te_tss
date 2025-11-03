perl -alne '$id=join"|",@F[0..4];$a=$F[3]-10;$b=$F[3]+10;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' high_TE_TSS_cluster_te.txt >! high_TE_TSS_cluster_10.txt
perl -alne '$id=join"|",@F[0..4];$a=$F[3]-100;$b=$F[3]+100;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' high_TE_TSS_cluster_te.txt >! high_TE_TSS_cluster_100.txt
perl -alne '$id=join"|",@F[0..4];$a=$F[3]-10;$b=$F[3]+10;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' TE_TSS_cluster_te.txt >! TE_TSS_cluster_10.txt
perl -alne '$id=join"|",@F[0..4];$a=$F[3]-100;$b=$F[3]+100;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' TE_TSS_cluster_te.txt >! TE_TSS_cluster_100.txt

## extract TE_TSS only within TE region
perl -alne '$site=$F[3];$flag=0;for(split /,/,$F[-1]){@a=split /\|/,$_;$flag++ if $site>=$a[1] and $site<=$a[2]};print if $flag' high_TE_TSS_cluster_te.txt >! high_TE_TSS_cluster_te_withinregion.txt
perl -alne '$site=$F[3];$flag=0;for(split /,/,$F[-1]){@a=split /\|/,$_;$flag++ if $site>=$a[1] and $site<=$a[2]};print if $flag' TE_TSS_cluster_te.txt >! TE_TSS_cluster_te_withinregion.txt

for i in {BosTau9,CalJac3,CanFam3,EquCab3,GorGor5,MicMur3,Mm10,NomLeu3,OtoGar3,PanTro6,PapAnu4,PonAbe3,RheMac10,Rn6,SusScr11,TarSyr2};do liftOver -minMatch=0.9 -minChainT=10000 -minChainQ=10000 -multiple high_TE_TSS_cluster_10.txt hg38To$i.over.chain.gz high_TE_TSS_cluster_10_$i.txt high_TE_TSS_cluster_10_${i}_unmap.txt;done &
for i in {BosTau9,CalJac3,CanFam3,EquCab3,GorGor5,MicMur3,Mm10,NomLeu3,OtoGar3,PanTro6,PapAnu4,PonAbe3,RheMac10,Rn6,SusScr11,TarSyr2};do liftOver -minMatch=0.9 -minChainT=10000 -minChainQ=10000 -multiple TE_TSS_cluster_10.txt hg38To$i.over.chain.gz TE_TSS_cluster_10_$i.txt TE_TSS_cluster_10_${i}_unmap.txt;done &

for i in {BosTau9,CalJac3,CanFam3,EquCab3,GorGor5,MicMur3,Mm10,NomLeu3,OtoGar3,PanTro6,PapAnu4,PonAbe3,RheMac10,Rn6,SusScr11,TarSyr2};do liftOver -minMatch=0.5 -minChainT=10000 -minChainQ=10000 -multiple high_TE_TSS_cluster_100.txt hg38To$i.over.chain.gz high_TE_TSS_cluster_100_$i.txt high_TE_TSS_cluster_100_${i}_unmap.txt;done &
for i in {BosTau9,CalJac3,CanFam3,EquCab3,GorGor5,MicMur3,Mm10,NomLeu3,OtoGar3,PanTro6,PapAnu4,PonAbe3,RheMac10,Rn6,SusScr11,TarSyr2};do liftOver -minMatch=0.5 -minChainT=10000 -minChainQ=10000 -multiple TE_TSS_cluster_100.txt hg38To$i.over.chain.gz TE_TSS_cluster_100_$i.txt TE_TSS_cluster_100_${i}_unmap.txt;done &


## mammalian
for i in {CanFam3,EquCab3,BosTau9,SusScr11,Rn6,Mm10};do cut -f 4 TE_TSS_cluster_10_$i.txt;done|sort -u >! TE_TSS_cluster_10_mammalian_id.txt
for i in {CanFam3,EquCab3,BosTau9,SusScr11,Rn6,Mm10};do cut -f 4 TE_TSS_cluster_100_$i.txt;done|sort -u >! TE_TSS_cluster_100_mammalian_id.txt
perl -alne 'print if $a{$_}++' TE_TSS_cluster_10_mammalian_id.txt TE_TSS_cluster_100_mammalian_id.txt >! TE_TSS_cluster_mammalian_id.txt
## primate
for i in {OtoGar3,MicMur3,TarSyr2,CalJac3};do cut -f 4 TE_TSS_cluster_10_$i.txt;done|sort -u >! TE_TSS_cluster_10_primate_id.txt
for i in {OtoGar3,MicMur3,TarSyr2,CalJac3};do cut -f 4 TE_TSS_cluster_100_$i.txt;done|sort -u >! TE_TSS_cluster_100_primate_id.txt
perl -alne 'print if $a{$_}++' TE_TSS_cluster_10_primate_id.txt TE_TSS_cluster_100_primate_id.txt|perl -alne 'if($ARGV=~/mammalian/){$a{$_}++}else{print if not $a{$_}}' TE_TSS_cluster_mammalian_id.txt - >! TE_TSS_cluster_primate_id.txt
## OWA
for i in {PapAnu4,RheMac10,NomLeu3};do cut -f 4 TE_TSS_cluster_10_$i.txt;done|sort -u >! TE_TSS_cluster_10_OWA_id.txt
for i in {PapAnu4,RheMac10,NomLeu3};do cut -f 4 TE_TSS_cluster_100_$i.txt;done|sort -u >! TE_TSS_cluster_100_OWA_id.txt
perl -alne 'print if $a{$_}++' TE_TSS_cluster_10_OWA_id.txt TE_TSS_cluster_100_OWA_id.txt|perl -alne 'if($ARGV=~/mammalian|primate/){$a{$_}++}else{print if not $a{$_}}' TE_TSS_cluster_mammalian_id.txt TE_TSS_cluster_primate_id.txt - >! TE_TSS_cluster_OWA_id.txt
## hominid
for i in {PonAbe3,GorGor5,PanTro6};do cut -f 4 TE_TSS_cluster_10_$i.txt;done|sort -u >! TE_TSS_cluster_10_hominid_id.txt
for i in {PonAbe3,GorGor5,PanTro6};do cut -f 4 TE_TSS_cluster_100_$i.txt;done|sort -u >! TE_TSS_cluster_100_hominid_id.txt
perl -alne 'print if $a{$_}++' TE_TSS_cluster_10_hominid_id.txt TE_TSS_cluster_100_hominid_id.txt|perl -alne 'if($ARGV=~/mammalian|primate|OWA/){$a{$_}++}else{print if not $a{$_}}' TE_TSS_cluster_mammalian_id.txt TE_TSS_cluster_primate_id.txt TE_TSS_cluster_OWA_id.txt - >! TE_TSS_cluster_hominid_id.txt
## human
perl -alne 'if($#F==0){$a{$_}++}else{print $F[3] if not $a{$F[3]}}' TE_TSS_cluster_mammalian_id.txt TE_TSS_cluster_primate_id.txt TE_TSS_cluster_OWA_id.txt TE_TSS_cluster_hominid_id.txt TE_TSS_cluster_10.txt >! TE_TSS_cluster_human_id.txt
## classes
perl -alne 'if($#F==0){$ARGV=~/TE_TSS_cluster_(.*)_id.txt/;$n=$1;$a{$_}=$n}else{$"="|";$id="@F[0..4]";print "$_\t$a{$id}"}' TE_TSS_cluster_{mammalian,primate,OWA,hominid,human}_id.txt TE_TSS_cluster_te_withinregion.txt >! TE_TSS_cluster_te_labeled.txt

## mammalian
for i in {CanFam3,EquCab3,BosTau9,SusScr11,Rn6,Mm10};do cut -f 4 high_TE_TSS_cluster_10_$i.txt;done|sort -u >! high_TE_TSS_cluster_10_mammalian_id.txt
for i in {CanFam3,EquCab3,BosTau9,SusScr11,Rn6,Mm10};do cut -f 4 high_TE_TSS_cluster_100_$i.txt;done|sort -u >! high_TE_TSS_cluster_100_mammalian_id.txt
perl -alne 'print if $a{$_}++' high_TE_TSS_cluster_10_mammalian_id.txt high_TE_TSS_cluster_100_mammalian_id.txt >! high_TE_TSS_cluster_mammalian_id.txt
## primate
for i in {OtoGar3,MicMur3,TarSyr2,CalJac3};do cut -f 4 high_TE_TSS_cluster_10_$i.txt;done|sort -u >! high_TE_TSS_cluster_10_primate_id.txt
for i in {OtoGar3,MicMur3,TarSyr2,CalJac3};do cut -f 4 high_TE_TSS_cluster_100_$i.txt;done|sort -u >! high_TE_TSS_cluster_100_primate_id.txt
perl -alne 'print if $a{$_}++' high_TE_TSS_cluster_10_primate_id.txt high_TE_TSS_cluster_100_primate_id.txt|perl -alne 'if($ARGV=~/mammalian/){$a{$_}++}else{print if not $a{$_}}' high_TE_TSS_cluster_mammalian_id.txt - >! high_TE_TSS_cluster_primate_id.txt
## OWA
for i in {PapAnu4,RheMac10,NomLeu3};do cut -f 4 high_TE_TSS_cluster_10_$i.txt;done|sort -u >! high_TE_TSS_cluster_10_OWA_id.txt
for i in {PapAnu4,RheMac10,NomLeu3};do cut -f 4 high_TE_TSS_cluster_100_$i.txt;done|sort -u >! high_TE_TSS_cluster_100_OWA_id.txt
perl -alne 'print if $a{$_}++' high_TE_TSS_cluster_10_OWA_id.txt high_TE_TSS_cluster_100_OWA_id.txt|perl -alne 'if($ARGV=~/mammalian|primate/){$a{$_}++}else{print if not $a{$_}}' high_TE_TSS_cluster_mammalian_id.txt high_TE_TSS_cluster_primate_id.txt - >! high_TE_TSS_cluster_OWA_id.txt
## hominid
for i in {PonAbe3,GorGor5,PanTro6};do cut -f 4 high_TE_TSS_cluster_10_$i.txt;done|sort -u >! high_TE_TSS_cluster_10_hominid_id.txt
for i in {PonAbe3,GorGor5,PanTro6};do cut -f 4 high_TE_TSS_cluster_100_$i.txt;done|sort -u >! high_TE_TSS_cluster_100_hominid_id.txt
perl -alne 'print if $a{$_}++' high_TE_TSS_cluster_10_hominid_id.txt high_TE_TSS_cluster_100_hominid_id.txt|perl -alne 'if($ARGV=~/mammalian|primate|OWA/){$a{$_}++}else{print if not $a{$_}}' high_TE_TSS_cluster_mammalian_id.txt high_TE_TSS_cluster_primate_id.txt high_TE_TSS_cluster_OWA_id.txt - >! high_TE_TSS_cluster_hominid_id.txt
## human
perl -alne 'if($#F==0){$a{$_}++}else{print $F[3] if not $a{$F[3]}}' high_TE_TSS_cluster_mammalian_id.txt high_TE_TSS_cluster_primate_id.txt high_TE_TSS_cluster_OWA_id.txt high_TE_TSS_cluster_hominid_id.txt high_TE_TSS_cluster_10.txt >! high_TE_TSS_cluster_human_id.txt
## classes
perl -alne 'if($#F==0){$ARGV=~/high_TE_TSS_cluster_(.*)_id.txt/;$n=$1;$a{$_}=$n}else{$"="|";$id="@F[0..4]";print "$_\t$a{$id}"}' high_TE_TSS_cluster_{mammalian,primate,OWA,hominid,human}_id.txt high_TE_TSS_cluster_te_withinregion.txt >! high_TE_TSS_cluster_te_labeled.txt

## TSS_1k TE
bedtools intersect -a $zlab/genome/hg38/te.bed -b TSS_1k.bed -wa -wb|cut -f 1-6|sort -u|perl -alne '$id=join"|",@F;$"="\t";print "@F[0..2]\t$id\t0\t$F[5]"' >! TSS_1k_te.bed
for i in {BosTau9,CalJac3,CanFam3,EquCab3,GorGor5,MicMur3,Mm10,NomLeu3,OtoGar3,PanTro6,PapAnu4,PonAbe3,RheMac10,Rn6,SusScr11,TarSyr2};do liftOver -minMatch=0.5 -minChainT=10000 -minChainQ=10000 -multiple TSS_1k_te.bed hg38To$i.over.chain.gz TSS_1k_te_$i.bed TSS_1k_te_${i}_unmap.txt &;done
for i in {CanFam3,EquCab3,BosTau9,SusScr11,Rn6,Mm10};do cut -f 4 TSS_1k_te_$i.bed;done|sort -u >! TSS_1k_te_mammalian_id.txt
for i in {OtoGar3,MicMur3,TarSyr2,CalJac3};do           cut -f 4 TSS_1k_te_$i.bed;done|sort -u >! TSS_1k_te_primate_id.txt
for i in {PapAnu4,RheMac10,NomLeu3};do                  cut -f 4 TSS_1k_te_$i.bed;done|sort -u >! TSS_1k_te_OWA_id.txt
for i in {PonAbe3,GorGor5,PanTro6};do                   cut -f 4 TSS_1k_te_$i.bed;done|sort -u >! TSS_1k_te_hominid_id.txt
perl -alne 'if($#F==0){$a{$_}="hominid" if $ARGV=~/hominid/;$a{$_}="OWA" if $ARGV=~/OWA/;$a{$_}="primate" if $ARGV=~/primate/;$a{$_}="mammalian" if $ARGV=~/mammalian/}else{$id=$F[3];$t=$a{$id}?$a{$id}:"human";print "$_\t$t"}' TSS_1k_te_hominid_id.txt TSS_1k_te_OWA_id.txt TSS_1k_te_primate_id.txt TSS_1k_te_mammalian_id.txt TSS_1k_te.bed >! TSS_1k_te_labled.bed

perl -alne '%a=map {(split /\|/,$_)[0]=>1} split /,/,$F[7];$n=%a;$n=$n>=3?3:$n;print "$n\t$F[-1]"' TE_TSS_cluster_te_labeled.txt|sort |uniq -c|sort -k 2,2 -k 1,1nr

perl -alne 'if($#F==6){$"="|";$te{"@F[0..5]"}=$F[6]}else{$x=(split /,/,$F[8])[0];print "all\t$F[-1]\t$te{$x}"}' ~zlab/genome/hg38/rmsk_div.bed TE_TSS_cluster_te_labeled.txt >! te_div.txt
perl -alne 'if($#F==6){$"="|";$te{"@F[0..5]"}=$F[6]}else{$x=(split /,/,$F[8])[0];print "high\t$F[-1]\t$te{$x}"}' ~zlab/genome/hg38/rmsk_div.bed high_TE_TSS_cluster_te_labeled.txt >> te_div.txt
perl -alne 'if($ARGV=~/div/){$"="|";$te{"@F[0..5]"}=$F[6]}else{$x=$F[3];print "genomic\t$F[-1]\t$te{$x}"}' ~zlab/genome/hg38/rmsk_div.bed TSS_1k_te_labled.bed >> te_div.txt
perl -alne 'if($F[1] eq "human"){$F[1]="hominid"};$"="\t";print "@F"' te_div.txt >! te_div_new.txt
pdf('te_div_new.pdf', width=10, height=5)
ggplot(a, aes(factor(V2, levels=c('mammalian', 'primate', 'OWA', 'hominid')), V3, fill=factor(V2, levels=c('mammalian', 'primate', 'OWA', 'hominid')), alpha=factor(V1, levels=c('genomic', 'all')))) + geom_violin(position=position_dodge(1)) + theme_classic() + coord_flip() + theme(legend.position="NA") + scale_fill_manual(values=c('black', '#FF9300', '#0076BA', '#4349AA')) + scale_alpha_manual(values=c(0.5, 1)) + geom_boxplot(outlier.shape=NA, width=0.2, fill='white', position=position_dodge(1))
dev.off()
perl -alne 'if($#F==6){$"="|";$te{"@F[0..5]"}=$F[6]}else{if(/testis/){$l="testis"}else{$l="non_testis"};$x=(split /,/,$F[8])[0];print "$l\t$F[-1]\t$te{$x}"}' ~zlab/genome/hg38/rmsk_div.bed TE_TSS_cluster_te_labeled.txt >! testis_te_div.txt
perl -alne 'if($F[1] eq "human"){$F[1]="hominid"};$"="\t";print "@F"' testis_te_div.txt >! testis_te_div_new.txt
pdf('testis_te_div_new.pdf', width=10, height=5)
ggplot(a, aes(factor(V2, levels=c('mammalian', 'primate', 'OWA', 'hominid')), V3, fill=factor(V2, levels=c('mammalian', 'primate', 'OWA', 'hominid')), alpha=factor(V1, levels=c('non_testis', 'testis')))) + geom_violin(position=position_dodge(1)) + theme_classic() + coord_flip() + theme(legend.position="NA") + scale_fill_manual(values=c('black', '#FF9300', '#0076BA', '#4349AA')) + scale_alpha_manual(values=c(0.5, 1)) + geom_boxplot(outlier.shape=NA, width=0.2, fill='white', position=position_dodge(1))
dev.off()
perl -alne '$l=(split /\|/,$F[8])[4];print "$l\t$F[-1]"' TE_TSS_cluster_te_labeled.txt|sort|uniq -c|sort -k 3,3 -k 1,1nr     |perl -alne 'print "$F[2]\t$F[1]\t$F[0]"'|perl -alne '$F[1]=~s/\?//;$F[1]=$F[1]=~/^(Alu|MIR|L1|L2|ERV)/?$F[1]:"other";$"="\t";$a{"$F[0]\t$F[1]"}+=$F[2];END{print "$_\t$a{$_}" for keys %a}' >! TE_TSS_cluster_te_class.txt
perl -alne '$l=(split /\|/,$F[8])[4];print "$l\t$F[-1]"' high_TE_TSS_cluster_te_labeled.txt|sort|uniq -c|sort -k 3,3 -k 1,1nr|perl -alne 'print "$F[2]\t$F[1]\t$F[0]"'|perl -alne '$F[1]=~s/\?//;$F[1]=$F[1]=~/^(Alu|MIR|L1|L2|ERV)/?$F[1]:"other";$"="\t";$a{"$F[0]\t$F[1]"}+=$F[2];END{print "$_\t$a{$_}" for keys %a}' >! high_TE_TSS_cluster_te_class.txt
perl -alne '$l=(split /\|/,$F[3])[4];print "$l\t$F[-1]"' TSS_1k_te_labled.bed|sort|uniq -c|sort -k 3,3 -k 1,1nr              |perl -alne 'print "$F[2]\t$F[1]\t$F[0]"'|perl -alne '$F[1]=~s/\?//;$F[1]=$F[1]=~/^(Alu|MIR|L1|L2|ERV)/?$F[1]:"other";$"="\t";$a{"$F[0]\t$F[1]"}+=$F[2];END{print "$_\t$a{$_}" for keys %a}' >! TSS_1k_te_class.txt
perl -alne 'print if /testis/' TE_TSS_cluster_te_labeled.txt >! testis_TE_TSS_cluster_te_labeled.txt
perl -alne 'print if !/testis/' TE_TSS_cluster_te_labeled.txt >! non_testis_TE_TSS_cluster_te_labeled.txt
perl -alne '$l=(split /\|/,$F[8])[4];print "$l\t$F[-1]"' testis_TE_TSS_cluster_te_labeled.txt|sort|uniq -c|sort -k 3,3 -k 1,1nr     |perl -alne 'print "$F[2]\t$F[1]\t$F[0]"'|perl -alne '$F[1]=~s/\?//;$F[1]=$F[1]=~/^(Alu|MIR|L1|L2|ERV)/?$F[1]:"other";$"="\t";$a{"$F[0]\t$F[1]"}+=$F[2];END{print "$_\t$a{$_}" for keys %a}' >! testis_TE_TSS_cluster_te_class.txt
perl -alne '$l=(split /\|/,$F[8])[4];print "$l\t$F[-1]"' non_testis_TE_TSS_cluster_te_labeled.txt|sort|uniq -c|sort -k 3,3 -k 1,1nr     |perl -alne 'print "$F[2]\t$F[1]\t$F[0]"'|perl -alne '$F[1]=~s/\?//;$F[1]=$F[1]=~/^(Alu|MIR|L1|L2|ERV)/?$F[1]:"other";$"="\t";$a{"$F[0]\t$F[1]"}+=$F[2];END{print "$_\t$a{$_}" for keys %a}' >! non_testis_TE_TSS_cluster_te_class.txt

## canonical TSS
perl -alne '$id=join"|",@F[0..2];$a=$F[1]-10;$b=$F[1]+10;  next if $a<0 or $F[0]!~/chr(\d{1,2}|X|Y)$/;print "$F[0]\t$a\t$b\t$id"' GENCODE24_basic_tss.bed|bedtools intersect -a - -b TE_TSS_cluster_10.txt -v >! GENCODE24_basic_tss_10.txt
perl -alne '$id=join"|",@F[0..2];$a=$F[1]-100;$b=$F[1]+100;next if $a<0 or $F[0]!~/chr(\d{1,2}|X|Y)$/;print "$F[0]\t$a\t$b\t$id"' GENCODE24_basic_tss.bed|bedtools intersect -a - -b TE_TSS_cluster_100.txt -v >! GENCODE24_basic_tss_100.txt
for i in {BosTau9,CalJac3,CanFam3,EquCab3,GorGor5,MicMur3,Mm10,NomLeu3,OtoGar3,PanTro6,PapAnu4,PonAbe3,RheMac10,Rn6,SusScr11,TarSyr2};do liftOver -minMatch=0.9 -minChainT=10000 -minChainQ=10000 -multiple GENCODE24_basic_tss_10.txt hg38To$i.over.chain.gz GENCODE24_basic_tss_10_$i.txt GENCODE24_basic_tss_10_${i}_unmap.txt;done &
for i in {BosTau9,CalJac3,CanFam3,EquCab3,GorGor5,MicMur3,Mm10,NomLeu3,OtoGar3,PanTro6,PapAnu4,PonAbe3,RheMac10,Rn6,SusScr11,TarSyr2};do liftOver -minMatch=0.5 -minChainT=10000 -minChainQ=10000 -multiple GENCODE24_basic_tss_100.txt hg38To$i.over.chain.gz GENCODE24_basic_tss_100_$i.txt GENCODE24_basic_tss_100_${i}_unmap.txt;done &

for i in {CanFam3,EquCab3,BosTau9,SusScr11,Rn6,Mm10};do cut -f 4 GENCODE24_basic_tss_10_$i.txt;done|sort -u >! GENCODE24_basic_tss_10_mammalian_id.txt
for i in {CanFam3,EquCab3,BosTau9,SusScr11,Rn6,Mm10};do cut -f 4 GENCODE24_basic_tss_100_$i.txt;done|sort -u >! GENCODE24_basic_tss_100_mammalian_id.txt
perl -alne 'print if $a{$_}++' GENCODE24_basic_tss_10_mammalian_id.txt GENCODE24_basic_tss_100_mammalian_id.txt >! GENCODE24_basic_tss_mammalian_id.txt
for i in {OtoGar3,MicMur3,TarSyr2,CalJac3};do cut -f 4 GENCODE24_basic_tss_10_$i.txt;done|sort -u >! GENCODE24_basic_tss_10_primate_id.txt
for i in {OtoGar3,MicMur3,TarSyr2,CalJac3};do cut -f 4 GENCODE24_basic_tss_100_$i.txt;done|sort -u >! GENCODE24_basic_tss_100_primate_id.txt
perl -alne 'print if $a{$_}++' GENCODE24_basic_tss_10_primate_id.txt GENCODE24_basic_tss_100_primate_id.txt|perl -alne 'if($ARGV=~/mammalian/){$a{$_}++}else{print if not $a{$_}}' GENCODE24_basic_tss_mammalian_id.txt - >! GENCODE24_basic_tss_primate_id.txt
for i in {PapAnu4,RheMac10,NomLeu3};do cut -f 4 GENCODE24_basic_tss_10_$i.txt;done|sort -u >! GENCODE24_basic_tss_10_OWA_id.txt
for i in {PapAnu4,RheMac10,NomLeu3};do cut -f 4 GENCODE24_basic_tss_100_$i.txt;done|sort -u >! GENCODE24_basic_tss_100_OWA_id.txt
perl -alne 'print if $a{$_}++' GENCODE24_basic_tss_10_OWA_id.txt GENCODE24_basic_tss_100_OWA_id.txt|perl -alne 'if($ARGV=~/mammalian|primate/){$a{$_}++}else{print if not $a{$_}}' GENCODE24_basic_tss_mammalian_id.txt GENCODE24_basic_tss_primate_id.txt - >! GENCODE24_basic_tss_OWA_id.txt
for i in {PonAbe3,GorGor5,PanTro6};do cut -f 4 GENCODE24_basic_tss_10_$i.txt;done|sort -u >! GENCODE24_basic_tss_10_hominid_id.txt
for i in {PonAbe3,GorGor5,PanTro6};do cut -f 4 GENCODE24_basic_tss_100_$i.txt;done|sort -u >! GENCODE24_basic_tss_100_hominid_id.txt
perl -alne 'print if $a{$_}++' GENCODE24_basic_tss_10_hominid_id.txt GENCODE24_basic_tss_100_hominid_id.txt|perl -alne 'if($ARGV=~/mammalian|primate|OWA/){$a{$_}++}else{print if not $a{$_}}' GENCODE24_basic_tss_mammalian_id.txt GENCODE24_basic_tss_primate_id.txt GENCODE24_basic_tss_OWA_id.txt - >! GENCODE24_basic_tss_hominid_id.txt
perl -alne 'if($#F==0){$a{$_}++}else{print $F[3] if not $a{$F[3]}}' GENCODE24_basic_tss_mammalian_id.txt GENCODE24_basic_tss_primate_id.txt GENCODE24_basic_tss_OWA_id.txt GENCODE24_basic_tss_hominid_id.txt GENCODE24_basic_tss_10.txt >! GENCODE24_basic_tss_human_id.txt
perl -alne 'if($#F==0){$ARGV=~/GENCODE24_basic_tss_(.*)_id.txt/;$n=$1;$a{$_}=$n}else{$"="|";$id="@F";print "$_\t$a{$id}" if $a{$id}}' GENCODE24_basic_tss_{mammalian,primate,OWA,hominid,human}_id.txt GENCODE29_basic_protein_coding_tss.bed >! GENCODE29_basic_protein_coding_tss_labeled.bed

perl -alne '@s=map{(split /\|/,$_)[0]} split /,/,$F[7];print "$F[0]\t$F[3]\t$F[4]\t".(join ",",@s)."\t$F[-1]"' TE_TSS_cluster_te_labeled.txt >! table3.txt

## for halLiftover
perl -alne '$id=join"|",@F[0..4];$a=$F[3]-50;$b=$F[3]+50;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' high_TE_TSS_cluster_te.txt >! high_TE_TSS_cluster_100_region.bed
perl -alne '$id=join"|",@F[0..4];$a=$F[3]-50;$b=$F[3]+50;print "$F[0]\t$a\t$b\t$id\t0\t$F[4]"' TE_TSS_cluster_te.txt >! TE_TSS_cluster_100_region.bed
perl -alne '$id=join"|",@F[0..2];$a=$F[1]-50;$b=$F[1]+50; next if $a<0 or $F[0]!~/chr(\d{1,2}|X|Y)$/;print "$F[0]\t$a\t$b\t$id"' GENCODE29_basic_protein_coding_tss.bed >! GENCODE29_basic_protein_coding_tss_100_region.bed
### another file is TSS_1k_te.bed
perl -alne 'if($#F==6){$"="|";$te{"@F[0..5]"}=$F[6]}else{$x=(split /,/,$F[8])[0];$t=(split /\|/,$x)[4];print "all\t$F[-1]\t$te{$x}\t$t"}' ~zlab/genome/hg38/rmsk_div.bed TE_TSS_cluster_te_labeled.txt >! te_div_family.txt
perl -alne 'if($#F==6){$"="|";$te{"@F[0..5]"}=$F[6]}else{$x=(split /,/,$F[8])[0];$t=(split /\|/,$x)[4];print "high\t$F[-1]\t$te{$x}\t$t"}' ~zlab/genome/hg38/rmsk_div.bed high_TE_TSS_cluster_te_labeled.txt >> te_div_family.txt
perl -alne 'if($ARGV=~/div/){$"="|";$te{"@F[0..5]"}=$F[6]}else{$x=$F[3];$t=(split /\|/,$F[3])[4];print "genomic\t$F[-1]\t$te{$x}\t$t"}' ~zlab/genome/hg38/rmsk_div.bed TSS_1k_te_labled.bed >> te_div_family.txt
perl -alne '$F[-1]=~s/\?//;$F[-1]=$F[-1]=~/^(Alu|MIR|L1|L2|ERV)/?$F[-1]:"other";print "$_\t$F[-1]"' te_div_family.txt >! te_div_family_final.txt
perl -alne 'if($#F==6){$"="|";$te{"@F[0..5]"}=$F[6]}else{if(/testis/){$l="testis"}else{$l="non_testis"};$x=(split /,/,$F[8])[0];$t=(split /\|/,$x)[4];print "$l\t$F[-1]\t$te{$x}\t$t"}' ~zlab/genome/hg38/rmsk_div.bed TE_TSS_cluster_te_labeled.txt|perl -alne '$F[-1]=~s/\?//;$F[-1]=$F[-1]=~/^(Alu|MIR|L1|L2|ERV)/?$F[-1]:"other";print "$_\t$F[-1]"' >! testis_te_div_family.txt
