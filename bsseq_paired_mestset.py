import os, sys
#ls N*.fastq.gz -1 > sample_codes.txt
#awk -F'_' '{print $1"_"$2" "$1"_"$2"_"$3"_R1_001"" "$1"_"$2"_"$3"_R2_001"}'  sample_codes.txt | sort -k1,1 -u  > sample_codesr1nr2.txt 
#paste sample_codesr1nr2.txt
#Follow these steps for preparing BiQ analyser inputs:
#Get names ref should not be rc for Mest

#ls -1 N*.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/mest_fastq/fastq_Laura2/fastas/ref_subMest_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/mest_fastq/fastq_Laura2/fastas/"$1"_R1.fa"}' > BiQsheetreplicate1.tsv
#ls -1 N*.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/mest_fastq/fastq_Laura2/fastas/ref_subMest_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/mest_fastq/fastq_Laura2/fastas/"$1"_R2.fa"}' > BiQsheetreplicate2.tsv
#ls -1 N*.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/mest_fastq/fastq_Laura2/fastas/ref_subMest_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/mest_fastq/fastq_Laura2/fastas/"$1"_R1_top.fa"}' > BiQsheetreplicate1_top.tsv
#ls -1 N*.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/mest_fastq/fastq_Laura2/fastas/ref_subMest_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/mest_fastq/fastq_Laura2/fastas/"$1"_R2_top.fa"}' > BiQsheetreplicate2_top.tsv

#Note:***************************** VERY IMPORTANT: ADD HEADER to EACH TSV Files otherwise it will not take the first file ****************************

list_data = ('''
N10_S10 N10_S10_L001_R1_001 N10_S10_L001_R2_001
N11_S11 N11_S11_L001_R1_001 N11_S11_L001_R2_001
N12_S12 N12_S12_L001_R1_001 N12_S12_L001_R2_001
N13_S13 N13_S13_L001_R1_001 N13_S13_L001_R2_001
N14_S14 N14_S14_L001_R1_001 N14_S14_L001_R2_001
N15_S15 N15_S15_L001_R1_001 N15_S15_L001_R2_001
N1_S1 N1_S1_L001_R1_001 N1_S1_L001_R2_001
N22_S64 N22_S64_L001_R1_001 N22_S64_L001_R2_001
N23_S16 N23_S16_L001_R1_001 N23_S16_L001_R2_001
N24_S17 N24_S17_L001_R1_001 N24_S17_L001_R2_001
N25_S18 N25_S18_L001_R1_001 N25_S18_L001_R2_001
N26_S19 N26_S19_L001_R1_001 N26_S19_L001_R2_001
N27_S20 N27_S20_L001_R1_001 N27_S20_L001_R2_001
N28_S21 N28_S21_L001_R1_001 N28_S21_L001_R2_001
N29_S22 N29_S22_L001_R1_001 N29_S22_L001_R2_001
N2_S2 N2_S2_L001_R1_001 N2_S2_L001_R2_001
N30_S23 N30_S23_L001_R1_001 N30_S23_L001_R2_001
N31_S24 N31_S24_L001_R1_001 N31_S24_L001_R2_001
N32_S25 N32_S25_L001_R1_001 N32_S25_L001_R2_001
N33_S26 N33_S26_L001_R1_001 N33_S26_L001_R2_001
N34_S27 N34_S27_L001_R1_001 N34_S27_L001_R2_001
N35_S28 N35_S28_L001_R1_001 N35_S28_L001_R2_001
N36_S29 N36_S29_L001_R1_001 N36_S29_L001_R2_001
N37_S30 N37_S30_L001_R1_001 N37_S30_L001_R2_001
N38_S31 N38_S31_L001_R1_001 N38_S31_L001_R2_001
N39_S32 N39_S32_L001_R1_001 N39_S32_L001_R2_001
N3_S3 N3_S3_L001_R1_001 N3_S3_L001_R2_001
N40_S33 N40_S33_L001_R1_001 N40_S33_L001_R2_001
N41_S34 N41_S34_L001_R1_001 N41_S34_L001_R2_001
N42_S35 N42_S35_L001_R1_001 N42_S35_L001_R2_001
N43_S36 N43_S36_L001_R1_001 N43_S36_L001_R2_001
N4_S4 N4_S4_L001_R1_001 N4_S4_L001_R2_001
N5_S5 N5_S5_L001_R1_001 N5_S5_L001_R2_001
N6_S6 N6_S6_L001_R1_001 N6_S6_L001_R2_001
N7_S7 N7_S7_L001_R1_001 N7_S7_L001_R2_001
N8_S8 N8_S8_L001_R1_001 N8_S8_L001_R2_001
N9_S9 N9_S9_L001_R1_001 N9_S9_L001_R2_001
NAN_S61 NAN_S61_L001_R1_001 NAN_S61_L001_R2_001
NA_S55 NA_S55_L001_R1_001 NA_S55_L001_R2_001
NB_S56 NB_S56_L001_R1_001 NB_S56_L001_R2_001
NC10_S45 NC10_S45_L001_R1_001 NC10_S45_L001_R2_001
NC11_S46 NC11_S46_L001_R1_001 NC11_S46_L001_R2_001
NC12_S47 NC12_S47_L001_R1_001 NC12_S47_L001_R2_001
NC13_S48 NC13_S48_L001_R1_001 NC13_S48_L001_R2_001
NC14_S49 NC14_S49_L001_R1_001 NC14_S49_L001_R2_001
NC15_S50 NC15_S50_L001_R1_001 NC15_S50_L001_R2_001
NC17_S51 NC17_S51_L001_R1_001 NC17_S51_L001_R2_001
NC18_S52 NC18_S52_L001_R1_001 NC18_S52_L001_R2_001
NC19_S53 NC19_S53_L001_R1_001 NC19_S53_L001_R2_001
NC1_S37 NC1_S37_L001_R1_001 NC1_S37_L001_R2_001
NC20_S54 NC20_S54_L001_R1_001 NC20_S54_L001_R2_001
NC3_S38 NC3_S38_L001_R1_001 NC3_S38_L001_R2_001
NC4_S39 NC4_S39_L001_R1_001 NC4_S39_L001_R2_001
NC5_S40 NC5_S40_L001_R1_001 NC5_S40_L001_R2_001
NC6_S41 NC6_S41_L001_R1_001 NC6_S41_L001_R2_001
NC7_S42 NC7_S42_L001_R1_001 NC7_S42_L001_R2_001
NC8_S43 NC8_S43_L001_R1_001 NC8_S43_L001_R2_001
NC9_S44 NC9_S44_L001_R1_001 NC9_S44_L001_R2_001
NC_S57 NC_S57_L001_R1_001 NC_S57_L001_R2_001
ND_S58 ND_S58_L001_R1_001 ND_S58_L001_R2_001
NFL_S62 NFL_S62_L001_R1_001 NFL_S62_L001_R2_001
NF_S59 NF_S59_L001_R1_001 NF_S59_L001_R2_001
NH_S60 NH_S60_L001_R1_001 NH_S60_L001_R2_001
NLA_S63 NLA_S63_L001_R1_001 NLA_S63_L001_R2_001
''').strip().split('\n')

print(list_data)

for listx in list_data:
    listx = listx.split(' ')
    folder = listx[0]
    read1 = listx[1]
    read2 = listx[2]

    print ('Sample in analysis: ' + folder)

    print ('Trimming data..: --phred33  (Sanger/Illumina 1.9+ encoding) for quality trimming. Default: ON, -q 20 (default),  BS-Seq' + folder)
    os.system('/home/ankitv/tools_av/TrimGalore-0.6.6/trim_galore  --phred33 --length 36 -q 20 --paired ' + read1 + '.fastq.gz '+ read2 + '.fastq.gz ')

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

    os.system('head ' + folder + '_R1.fa -n 8000 > ' + folder + '_R1_top.fa')
    os.system('head ' + folder + '_R2.fa -n 8000 > ' + folder + '_R2_top.fa')

    print ('Methylation extraction...')
    #os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark_methylation_extractor --bedGraph --include_overlap --zero_based --cutoff 10 ' + read1 + '_val_1_bismark.sortedByReadname.bam')

    print ('Analysis finished for ' + folder) 
    print ('\n')
    print ('\n')
