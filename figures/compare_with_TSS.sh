# Gro-cap
bedtools intersect -a K562_rampage_TE_TSS_high.bed -b K562_GROcap_peaks_rpm.txt -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t1GROcap"' >! K562_rampage_TE_TSS_high_upset.txt
bedtools intersect -a K562_rampage_TE_TSS_all.bed -b K562_GROcap_peaks_rpm.txt -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t1GROcap"' >! K562_rampage_TE_TSS_all_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_high.bed -b GM12878_GROcap_peaks_rpm.txt -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t1GROcap"' >! GM12878_rampage_TE_TSS_high_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_all.bed -b GM12878_GROcap_peaks_rpm.txt -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t1GROcap"' >! GM12878_rampage_TE_TSS_all_upset.txt

# PacBio
bedtools intersect -a K562_rampage_TE_TSS_high.bed -b K562_pacbio_peak.txt -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t2pacbio"' >> K562_rampage_TE_TSS_high_upset.txt
bedtools intersect -a K562_rampage_TE_TSS_all.bed -b K562_pacbio_peak.txt -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t2pacbio"' >> K562_rampage_TE_TSS_all_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_high.bed -b GM12878_pacbio_peak.txt -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t2pacbio"' >> GM12878_rampage_TE_TSS_high_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_all.bed -b GM12878_pacbio_peak.txt -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t2pacbio"' >> GM12878_rampage_TE_TSS_all_upset.txt

# CAGE
bedtools intersect -a K562_rampage_TE_TSS_high.bed -b hg38_fair+new_CAGE_peaks_phase1and2.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t3cage"' >> K562_rampage_TE_TSS_high_upset.txt
bedtools intersect -a K562_rampage_TE_TSS_all.bed -b hg38_fair+new_CAGE_peaks_phase1and2.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t3cage"' >> K562_rampage_TE_TSS_all_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_high.bed -b hg38_fair+new_CAGE_peaks_phase1and2.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t3cage"' >> GM12878_rampage_TE_TSS_high_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_all.bed -b hg38_fair+new_CAGE_peaks_phase1and2.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t3cage"' >> GM12878_rampage_TE_TSS_all_upset.txt

# rPEAK
bedtools intersect -a K562_rampage_TE_TSS_high.bed -b hg38-rPeaks.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t4rPeak"' >> K562_rampage_TE_TSS_high_upset.txt
bedtools intersect -a K562_rampage_TE_TSS_all.bed -b hg38-rPeaks.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t4rPeak"' >> K562_rampage_TE_TSS_all_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_high.bed -b hg38-rPeaks.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t4rPeak"' >> GM12878_rampage_TE_TSS_high_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_all.bed -b hg38-rPeaks.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t4rPeak"' >> GM12878_rampage_TE_TSS_all_upset.txt

# PLS
bedtools intersect -a K562_rampage_TE_TSS_high.bed -b GM12878.PLS.V3.bed -wa -wb|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t5PLS"' >> K562_rampage_TE_TSS_high_upset.txt
bedtools intersect -a K562_rampage_TE_TSS_all.bed -b GM12878.PLS.V3.bed -wa -wb|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t5PLS"' >> K562_rampage_TE_TSS_all_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_high.bed -b K562.PLS.V3.bed -wa -wb|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t5PLS"' >> GM12878_rampage_TE_TSS_high_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_all.bed -b K562.PLS.V3.bed -wa -wb|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t5PLS"' >> GM12878_rampage_TE_TSS_all_upset.txt

# TETSS
bedtools intersect -a K562_rampage_TE_TSS_high.bed -b human_TE_TSS.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t6database"' >> K562_rampage_TE_TSS_high_upset.txt
bedtools intersect -a K562_rampage_TE_TSS_all.bed -b human_TE_TSS.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t6database"' >> K562_rampage_TE_TSS_all_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_high.bed -b human_TE_TSS.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t6database"' >> GM12878_rampage_TE_TSS_high_upset.txt
bedtools intersect -a GM12878_rampage_TE_TSS_all.bed -b human_TE_TSS.bed -wa -wb -s|cut -f 1-7|sort -u|perl -alne '$"="|";print "@F\t6database"' >> GM12878_rampage_TE_TSS_all_upset.txt

perl -alne '$"="|";print "@F\t7TETSS"' K562_rampage_TE_TSS_high.bed >> K562_rampage_TE_TSS_high_upset.txt
perl -alne '$"="|";print "@F\t7TETSS"' K562_rampage_TE_TSS_all.bed >> K562_rampage_TE_TSS_all_upset.txt
perl -alne '$"="|";print "@F\t7TETSS"' GM12878_rampage_TE_TSS_high.bed >> GM12878_rampage_TE_TSS_high_upset.txt
perl -alne '$"="|";print "@F\t7TETSS"' GM12878_rampage_TE_TSS_all.bed >> GM12878_rampage_TE_TSS_all_upset.txt
