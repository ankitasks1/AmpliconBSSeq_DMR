# Amplicon BSSeq_DMR


#New analysis amplicon Dec 2021
#Python script for bsseq_paired_nnatset.py
#Convert to fasta bsseq_paired_nnatset_tofasta.py (later conversion step was added to bsseq_paired_nnatset.py also to make a complete pipeline from one script)
#Prepare BiQ Analyzer sheet
#Follow these steps:

#Use this way since Nnat was amplified on reverse strand I used minus strand as forward strand
#Prep for set1
#ls -1 *.gz | awk -F'_' '{print $2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/"$1"_R1.fa"}' | tail -n 5 > minus_SseriesBiQsheetR1.tsv
#ls -1 *.gz | awk -F'_' '{print $2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/"$1"_R2.fa"}' | tail -n 5 > minus_SseriesBiQsheetR2.tsv

#Prep for set2
#minus, since Nnat is amplified in reverse strand and also BiQ manual suggest to use the strand which was originally amplified by PCR
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R1.fa"}' > minus_BiQsheetreplicate1.tsv
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R2.fa"}' > minus_BiQsheetreplicate2.tsv

#Note:***************************** VERY IMPORTANT: ADD HEADER to EACH TSV Files otherwise it will not take the first file ****************************
#Format will look like (header is must)
alsp_minus_BiQsheetNnSeries_as2B_R1
S54_as2B	/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1.fa	/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/S54_as2B_R1.fa
S55_as2B	/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1.fa	/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/allele_specific_2/S55_as2B_R1.fa

#Get names ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1" "$1"_L001_R1_001"" "$1"_L001_R2_001"}'

#Import BiQ Analyzer Samplesheet to BiQ HiMOD


#Run BiQ Analyzer (Double click icon --> Create new project --> Click directory where you want result --> Load from data structure table (.tsv) --> press |> run)

#Export -> All results in one TSV -> table_results.tsv

#Make another file which contain information about group sample_color.txt eg. N15_S7	Downs	N15_S7%Downs

#Run Python on reference fasta to get CG positions, determine_CpG_position_reverse.py (reverse because Nnat was reverse)


#For allele specific data Set1 and set2 both were created in one folder so they NnS series.tsv
