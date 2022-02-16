#Prep for set1
ls -1 *.gz | awk -F'_' '{print $2}' | sort -k1,1 -u | awk '{print $1"_as1B""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/"$1"_B_g1_R2.fa"}' | tail -n 5 > alsp_minus_SseriesBiQsheet_g1_R2.tsv
ls -1 *.gz | awk -F'_' '{print $2}' | sort -k1,1 -u | awk '{print $1"_as2B""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/"$1"_B_g2_R2.fa"}' | tail -n 5 > alsp_minus_SseriesBiQsheet_g2_R2.tsv

#Prep for set2
#minus, since Nnat is amplified in reverse strand and also BiQ manual suggest to use the strand which was originally amplified by PCR
ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"_as1B""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific/""N"$1"_as1B_R2.fa"}' > alsp_minus_BiQsheetN_g1_R2.tsv
ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"_as2B""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific/""N"$1"_as2B_R2.fa"}' > alsp_minus_BiQsheetN_g2_R2.tsv

cat /media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/alsp_minus_SseriesBiQsheet_g1_R2.tsv ./../alsp_minus_BiQsheetN_g1_R2.tsv > alsp_minus_BiQsheetNSeries_g1_R2.tsv
cat /media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/alsp_minus_SseriesBiQsheet_g2_R2.tsv ./../alsp_minus_BiQsheetN_g2_R2.tsv > alsp_minus_BiQsheetNSeries_g2_R2.tsv

#Note:***************************** VERY IMPORTANT: ADD HEADER to EACH TSV Files otherwise it will not take the first file ****************************

