perl -alne 'print "GM12878\tgene\t$F[-1]\t$F[-2]"' GM12878_assigned_peaks.txt >! assigned_ratio.txt
perl -alne 'print "K562\tgene\t$F[-1]\t$F[-2]"' K562_assigned_peaks.txt >> assigned_ratio.txt
bedtools intersect -a GM12878_assigned_peaks.txt -b $zlab/genome/hg38/te.bed -wa -wb|cut -f 1-31|sort -u|perl -alne 'print "GM12878\tte_gene\t$F[-1]\t$F[-2]"' >> assigned_ratio.txt
bedtools intersect -a K562_assigned_peaks.txt -b $zlab/genome/hg38/te.bed -wa -wb|cut -f 1-31|sort -u|perl -alne 'print "K562\tte_gene\t$F[-1]\t$F[-2]"' >> assigned_ratio.txt
perl -alne 'print if $F[0]=~/chr(\d{1,2}|X|Y)$/ and $F[12]>=1 and $F[13]>=2.5 and abs($F[-1])<100' GM12878_assigned_peaks.txt >! GM12878_assigned_peaks_high.txt
perl -alne 'print if $F[0]=~/chr(\d{1,2}|X|Y)$/ and $F[12]>=1 and $F[13]>=2.5 and abs($F[-1])<100' K562_assigned_peaks.txt >! K562_assigned_peaks_high.txt

perl -alne 'print if $F[0]=~/chr(\d{1,2}|X|Y)$/ and $F[12]>=1 and $F[13]>=2.5 and $F[-1]<100 and $F[-2] !=1' testis_assigned_peaks.txt >! testis_assigned_peaks_high.txt
perl -alne 'print if $F[0]=~/chr(\d{1,2}|X|Y)$/ and $F[12]>=1 and $F[13]>=2.5 and $F[-1]<100 and $F[-2] !=1' liver_assigned_peaks.txt >! liver_assigned_peaks_high.txt
perl -alne 'print if $F[0]=~/chr(\d{1,2}|X|Y)$/ and $F[12]>=1 and $F[13]>=2.5 and $F[-1]<100 and $F[-2] !=1' heart_assigned_peaks.txt >! heart_assigned_peaks_high.txt
perl -alne 'print if $F[0]=~/chr(\d{1,2}|X|Y)$/ and $F[12]>=10 and $F[13]>=2.5 and $F[-1]<50' testis_assigned_peaks.txt >! testis_assigned_peaks_new.txt
perl -alne 'print if $F[0]=~/chr(\d{1,2}|X|Y)$/ and $F[12]>=10 and $F[13]>=2.5 and $F[-1]<50' liver_assigned_peaks.txt >! liver_assigned_peaks_new.txt
perl -alne 'print if $F[0]=~/chr(\d{1,2}|X|Y)$/ and $F[12]>=10 and $F[13]>=2.5 and $F[-1]<50' heart_assigned_peaks.txt >! heart_assigned_peaks_new.txt

perl -alne 'print "GM12878\tgene\t$F[-1]\t$F[-2]"' GM12878_assigned_peaks_high.txt >! assigned_ratio_new.txt
perl -alne 'print "testis\tgene\t$F[-2]\t$F[-1]"' testis_assigned_peaks_new.txt >> assigned_ratio_new.txt
perl -alne 'print "K562\tgene\t$F[-1]\t$F[-2]"' K562_assigned_peaks_high.txt >> assigned_ratio_new.txt
perl -alne 'print "liver\tgene\t$F[-2]\t$F[-1]"' liver_assigned_peaks_new.txt >> assigned_ratio_new.txt
perl -alne 'print "heart\tgene\t$F[-2]\t$F[-1]"' heart_assigned_peaks_new.txt >> assigned_ratio_new.txt
bedtools intersect -a GM12878_assigned_peaks_high.txt -b $zlab/genome/hg38/te.bed -wa -wb|cut -f 1-31|sort -u|perl -alne 'print "GM12878\tte_gene\t$F[-1]\t$F[-2]"' >> assigned_ratio_new.txt
bedtools intersect -a K562_assigned_peaks_high.txt -b $zlab/genome/hg38/te.bed -wa -wb|cut -f 1-31|sort -u|perl -alne 'print "K562\tte_gene\t$F[-1]\t$F[-2]"' >> assigned_ratio_new.txt
bedtools intersect -a testis_assigned_peaks_high.txt -b $zlab/genome/hg38/te.bed -wa -wb|cut -f 1-23|sort -u|perl -alne 'print "testis\tte_gene\t$F[-2]\t$F[-1]"' >> assigned_ratio_new.txt
bedtools intersect -a liver_assigned_peaks_high.txt -b $zlab/genome/hg38/te.bed -wa -wb|cut -f 1-23|sort -u|perl -alne 'print "liver\tte_gene\t$F[-2]\t$F[-1]"' >> assigned_ratio_new.txt
bedtools intersect -a heart_assigned_peaks_high.txt -b $zlab/genome/hg38/te.bed -wa -wb|cut -f 1-23|sort -u|perl -alne 'print "heart\tte_gene\t$F[-2]\t$F[-1]"' >> assigned_ratio_new.txt

# perl -alne 'print "GM12878\tother\t$F[-2]\t$F[-3]" if $F[-1] eq "other" and $F[-2] ne "None"' GM12878_rampage_TE_gene_all_labeled.txt >> assigned_ratio_new.txt
# perl -alne 'print "K562\tother\t$F[-2]\t$F[-3]" if $F[-1] eq "other" and $F[-2] ne "None"' K562_rampage_TE_gene_all_labeled.txt >> assigned_ratio_new.txt
# perl -alne 'print "testis\tother\t$F[-2]\t$F[-3]" if $F[-1] eq "other" and $F[-2] ne "None"' testis_rampage_TE_gene_all_labeled.txt >> assigned_ratio_new.txt
