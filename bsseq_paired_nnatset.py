import os, sys

#Follow these steps for preparing BiQ analyser inputs:
#Get names ref should always be rc for Nnat
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"_rc""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/ref_subNnat_R1_rc.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R1_rc.fa"}' > BiQsheetreplicate1_rc.tsv
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"_rc""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/ref_subNnat_R2_rc.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R2_rc.fa"}' > BiQsheetreplicate2_rc.tsv
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/ref_subNnat_R1_rc.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R1.fa"}' > BiQsheetreplicate1.tsv
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/ref_subNnat_R2_rc.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R2.fa"}' > BiQsheetreplicate2.tsv
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"_rc""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/ref_subNnat_R1_rc.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R1_rc_top.fa"}' > BiQsheetreplicate1_rc_top.tsv
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"_rc""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/ref_subNnat_R2_rc.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R2_rc_top.fa"}' > BiQsheetreplicate2_rc_top.tsv

#Prep for set1
#ls -1 *.gz | awk -F'_' '{print $2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/"$1"_R1.fa"}' | tail -n 5 > minus_SseriesBiQsheetR1.tsv
#ls -1 *.gz | awk -F'_' '{print $2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/fastq_Laura/fastq_Laura/"$1"_R2.fa"}' | tail -n 5 > minus_SseriesBiQsheetR2.tsv

#Prep for set2
#minus, since Nnat is amplified in reverse strand and also BiQ manual suggest to use the strand which was originally amplified by PCR
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R1.fa"}' > minus_BiQsheetreplicate1.tsv
#ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/minus_ref_subNnat_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/newfastq_laura/fastq_laura/fastas/""N"$1"_R2.fa"}' > minus_BiQsheetreplicate2.tsv

#Note:***************************** VERY IMPORTANT: ADD HEADER to EACH TSV Files otherwise it will not take the first file ****************************

#Get names ls -1 *.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print "N"$1" "$1"_L001_R1_001"" "$1"_L001_R2_001"}'
list_data = ('''
N15_S7 15_S7_L001_R1_001 15_S7_L001_R2_001
N18_S8 18_S8_L001_R1_001 18_S8_L001_R2_001
N19_S9 19_S9_L001_R1_001 19_S9_L001_R2_001
N1_S1 1_S1_L001_R1_001 1_S1_L001_R2_001
N20_S10 20_S10_L001_R1_001 20_S10_L001_R2_001
N21_S11 21_S11_L001_R1_001 21_S11_L001_R2_001
N22_S12 22_S12_L001_R1_001 22_S12_L001_R2_001
N23_S13 23_S13_L001_R1_001 23_S13_L001_R2_001
N24_S14 24_S14_L001_R1_001 24_S14_L001_R2_001
N27_S15 27_S15_L001_R1_001 27_S15_L001_R2_001
N28_S16 28_S16_L001_R1_001 28_S16_L001_R2_001
N31_S17 31_S17_L001_R1_001 31_S17_L001_R2_001
N32_S18 32_S18_L001_R1_001 32_S18_L001_R2_001
N33_S19 33_S19_L001_R1_001 33_S19_L001_R2_001
N34_S20 34_S20_L001_R1_001 34_S20_L001_R2_001
N35_S21 35_S21_L001_R1_001 35_S21_L001_R2_001
N36_S22 36_S22_L001_R1_001 36_S22_L001_R2_001
N37_S23 37_S23_L001_R1_001 37_S23_L001_R2_001
N38_S24 38_S24_L001_R1_001 38_S24_L001_R2_001
N39_S25 39_S25_L001_R1_001 39_S25_L001_R2_001
N3_S2 3_S2_L001_R1_001 3_S2_L001_R2_001
N40_S26 40_S26_L001_R1_001 40_S26_L001_R2_001
N41_S27 41_S27_L001_R1_001 41_S27_L001_R2_001
N42_S28 42_S28_L001_R1_001 42_S28_L001_R2_001
N43_S29 43_S29_L001_R1_001 43_S29_L001_R2_001
N4_S3 4_S3_L001_R1_001 4_S3_L001_R2_001
N5_S4 5_S4_L001_R1_001 5_S4_L001_R2_001
N7_S5 7_S5_L001_R1_001 7_S5_L001_R2_001
N8_S6 8_S6_L001_R1_001 8_S6_L001_R2_001
NAN_S45 AN_S45_L001_R1_001 AN_S45_L001_R2_001
NA_S40 A_S40_L001_R1_001 A_S40_L001_R2_001
NC10_S34 C10_S34_L001_R1_001 C10_S34_L001_R2_001
NC12_S35 C12_S35_L001_R1_001 C12_S35_L001_R2_001
NC13_S36 C13_S36_L001_R1_001 C13_S36_L001_R2_001
NC14_S37 C14_S37_L001_R1_001 C14_S37_L001_R2_001
NC18_S38 C18_S38_L001_R1_001 C18_S38_L001_R2_001
NC19_S39 C19_S39_L001_R1_001 C19_S39_L001_R2_001
NC1_S30 C1_S30_L001_R1_001 C1_S30_L001_R2_001
NC3_S31 C3_S31_L001_R1_001 C3_S31_L001_R2_001
NC4_S32 C4_S32_L001_R1_001 C4_S32_L001_R2_001
NC7_S33 C7_S33_L001_R1_001 C7_S33_L001_R2_001
NC_S41 C_S41_L001_R1_001 C_S41_L001_R2_001
ND_S42 D_S42_L001_R1_001 D_S42_L001_R2_001
NF_S43 F_S43_L001_R1_001 F_S43_L001_R2_001
NH_S44 H_S44_L001_R1_001 H_S44_L001_R2_001
NLA_S46 LA_S46_L001_R1_001 LA_S46_L001_R2_001
NN15_S7 N15_S7_L001_R1_001 N15_S7_L001_R2_001
NT1_S47 T1_S47_L001_R1_001 T1_S47_L001_R2_001
NUndetermined_S0 Undetermined_S0_L001_R1_001 Undetermined_S0_L001_R2_001
''').strip().split('\n')

print(list_data)

for listx in list_data:
    listx = listx.split(' ')
    folder = listx[0]
    read1 = listx[1]
    read2 = listx[2]

    print 'Sample in analysis: ' + folder

    print ('Trimming data..: --phred33  (Sanger/Illumina 1.9+ encoding) for quality trimming. Default: ON, -q 20 (default), clipping R2 as first 60-70 bases were lowQ,  BS-Seq' + folder)
    os.system('/home/ankitv/tools_av/TrimGalore-0.6.6/trim_galore  --phred33 --length 36 -q 20 --clip_R2 60 --paired ' + read1 + '.fastq.gz '+ read2 + '.fastq.gz ')

    print ('QC for trimmed data')
    os.system('~/tools_av/FastQC/fastqc ' + read1 + '_val_1.fq.gz')
    os.system('~/tools_av/FastQC/fastqc ' + read2 + '_val_2.fq.gz')

    print ('Aligning data using Bismark: ' + folder)
    os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark --score_min L,0,-0.6 --genome /home/ankitv/ref_av/hg38_bs/ -1 ' + read1 + '_val_1.fq.gz -2 ' + read2 + '_val_2.fq.gz --bowtie2 --path_to_bowtie2 /home/ankitv/tools_av/bowtie2-2.3.4.1-linux-x86_64/')  

    print ('As the sample is amplicon duplicates are not removed')

    print ('Give a ID to BAM...')
    os.system('cp ' + read1 + '_val_1_bismark_bt2_pe.bam ' + folder + '_bismark_bt2_pe.bam')

    print ('Sorting data...')
    os.system('samtools sort -o ' + folder + '_bismark_bt2_pe.sort.bam ' + folder + '_bismark_bt2_pe.bam')

    print ('Index data...')
    os.system('samtools index ' + folder + '_bismark_bt2_pe.sort.bam')

    print ('Sorting data by readname...')
    os.system('samtools sort -n ' + folder + '_bismark_bt2_pe.bam -o ' + folder + '_bismark.sortedByReadname.bam')

    print ('BAM to fasta')
    os.system('samtools fasta -1 ' + folder + '_R1.fa -2 ' + folder + '_R2.fa ' + folder + '_bismark_bt2_pe.sort.bam')
    #For Nnat Since R1 is in reverse complimentary it need to be bring to forward (R1 contains SNP)

    print ('Reverse complement')
    os.system('/home/ankitv/tools_av/seqtk/seqtk seq -r ' + folder + '_R1.fa > ' + folder + '_R1_rc.fa')
    os.system('/home/ankitv/tools_av/seqtk/seqtk seq -r ' + folder + '_R2.fa > ' + folder + '_R2_rc.fa')

    os.system('head ' + folder + '_R1_rc.fa -n 8000 > ' + folder + '_R1_rc_top.fa')
    os.system('head ' + folder + '_R2_rc.fa -n 8000 > ' + folder + '_R2_rc_top.fa')

    print ('Methylation extraction...')
    #os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark_methylation_extractor --bedGraph --include_overlap --zero_based --cutoff 10 ' + read1 + '_val_1_bismark.sortedByReadname.bam')

    print 'Analysis finished for ' + folder 
    print '\n'
    print '\n'
