import os, sys
#ls -1 * | awk -F'_fastq' '{print $1}' | awk -F'.fastq' '{print $1}'  |  sort -k1,1 -u | awk '{print "cp"" "$1".fastq.gz"" ""N"$1".fastq.gz"}'
#ls -1 N* | awk -F'_fastq' '{print $1}' | awk -F'.fastq' '{print $1}'  |  sort -k1,1 -u | awk -F'_L001' '{print $1}'  | sort -k1,1 -u | awk '{print $1" "$1"_L001_R1_001"" "$1"_L001_R2_001"}'

#Follow these steps for preparing BiQ analyser inputs:
#Get names ref should not be rc for Snurf

#ls -1 N*fastq.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/snurf_fastq/fastqs/fastas/ref_subSnurf_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/snurf_fastq/fastqs/fastas/"$1"_R1.fa"}' > BiQsheetreplicate1.tsv
#ls -1 N*fastq.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/snurf_fastq/fastqs/fastas/ref_subSnurf_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/snurf_fastq/fastqs/fastas/"$1"_R2.fa"}' > BiQsheetreplicate2.tsv
#ls -1 N*fastq.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/snurf_fastq/fastqs/fastas/minus_ref_subSnurf_R1.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/snurf_fastq/fastqs/fastas/"$1"_R1.fa"}' > minus_BiQsheetreplicate1.tsv
#ls -1 N*fastq.gz | awk -F'_' '{print $1"_"$2}' | sort -k1,1 -u | awk '{print $1"\t""/media/ankitv/Archivio2/ankit/bs-seq/human/snurf_fastq/fastqs/fastas/minus_ref_subSnurf_R2.fa""\t""/media/ankitv/Archivio2/ankit/bs-seq/human/snurf_fastq/fastqs/fastas/"$1"_R2.fa"}' > minus_BiQsheetreplicate2.tsv



#Note:***************************** VERY IMPORTANT: ADD HEADER to EACH TSV Files otherwise it will not take the first file ****************************

list_data = ('''
N10_S2 N10_S2_L001_R1_001 N10_S2_L001_R2_001
N11_S3 N11_S3_L001_R1_001 N11_S3_L001_R2_001
N12_S4 N12_S4_L001_R1_001 N12_S4_L001_R2_001
N14_S5 N14_S5_L001_R1_001 N14_S5_L001_R2_001
N15_S6 N15_S6_L001_R1_001 N15_S6_L001_R2_001
N23_S7 N23_S7_L001_R1_001 N23_S7_L001_R2_001
N24_S8 N24_S8_L001_R1_001 N24_S8_L001_R2_001
N25_S9 N25_S9_L001_R1_001 N25_S9_L001_R2_001
N26_S10 N26_S10_L001_R1_001 N26_S10_L001_R2_001
N27_S11 N27_S11_L001_R1_001 N27_S11_L001_R2_001
N28_S12 N28_S12_L001_R1_001 N28_S12_L001_R2_001
N29_S13 N29_S13_L001_R1_001 N29_S13_L001_R2_001
N30_S14 N30_S14_L001_R1_001 N30_S14_L001_R2_001
N31_S15 N31_S15_L001_R1_001 N31_S15_L001_R2_001
N32_S16 N32_S16_L001_R1_001 N32_S16_L001_R2_001
N33_S17 N33_S17_L001_R1_001 N33_S17_L001_R2_001
N34_S18 N34_S18_L001_R1_001 N34_S18_L001_R2_001
N35_S19 N35_S19_L001_R1_001 N35_S19_L001_R2_001
N37_S20 N37_S20_L001_R1_001 N37_S20_L001_R2_001
N38_S21 N38_S21_L001_R1_001 N38_S21_L001_R2_001
N39_S22 N39_S22_L001_R1_001 N39_S22_L001_R2_001
N41_S23 N41_S23_L001_R1_001 N41_S23_L001_R2_001
N42_S24 N42_S24_L001_R1_001 N42_S24_L001_R2_001
N43_S25 N43_S25_L001_R1_001 N43_S25_L001_R2_001
N9_S1 N9_S1_L001_R1_001 N9_S1_L001_R2_001
NAN_S45 NAN_S45_L001_R1_001 NAN_S45_L001_R2_001
NA_S40 NA_S40_L001_R1_001 NA_S40_L001_R2_001
NB_S41 NB_S41_L001_R1_001 NB_S41_L001_R2_001
NC10_S31 NC10_S31_L001_R1_001 NC10_S31_L001_R2_001
NC12_S32 NC12_S32_L001_R1_001 NC12_S32_L001_R2_001
NC13_S33 NC13_S33_L001_R1_001 NC13_S33_L001_R2_001
NC14_S34 NC14_S34_L001_R1_001 NC14_S34_L001_R2_001
NC15_S35 NC15_S35_L001_R1_001 NC15_S35_L001_R2_001
NC17_S36 NC17_S36_L001_R1_001 NC17_S36_L001_R2_001
NC18_S37 NC18_S37_L001_R1_001 NC18_S37_L001_R2_001
NC19_S38 NC19_S38_L001_R1_001 NC19_S38_L001_R2_001
NC20_S39 NC20_S39_L001_R1_001 NC20_S39_L001_R2_001
NC3_S26 NC3_S26_L001_R1_001 NC3_S26_L001_R2_001
NC6_S27 NC6_S27_L001_R1_001 NC6_S27_L001_R2_001
NC7_S28 NC7_S28_L001_R1_001 NC7_S28_L001_R2_001
NC8_S29 NC8_S29_L001_R1_001 NC8_S29_L001_R2_001
NC9_S30 NC9_S30_L001_R1_001 NC9_S30_L001_R2_001
NC_S42 NC_S42_L001_R1_001 NC_S42_L001_R2_001
NF_S43 NF_S43_L001_R1_001 NF_S43_L001_R2_001
NH_S44 NH_S44_L001_R1_001 NH_S44_L001_R2_001
NLA_S46 NLA_S46_L001_R1_001 NLA_S46_L001_R2_001
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
