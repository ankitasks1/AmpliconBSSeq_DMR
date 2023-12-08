import os, sys

#fgrep -f NNAT_snp.txt /home/ankitv/ref_av/N-masked_hg38/commonSNPrefalt_allalt_chr_GRCh38_set_Acomb_Snpfile.txt -w | awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5}' > NNAT_commonSNPrefalt_GRCh38_Acomb_Snpfile.txt
#fgrep -f NNAT_snp.txt /home/ankitv/ref_av/N-masked_hg38/commonSNPrefalt_allalt_chr_GRCh38_set_Bcomb_Snpfile.txt -w | awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5}' > NNAT_commonSNPrefalt_GRCh38_Bcomb_Snpfile.txt


list_data = ('''
S54 54_S54_L001_R1_001 54_S54_L001_R2_001
S55 55_S55_L001_R1_001 55_S55_L001_R2_001
S56 56_S56_L001_R1_001 56_S56_L001_R2_001
S57 57_S57_L001_R1_001 57_S57_L001_R2_001
S58 58_S58_L001_R1_001 58_S58_L001_R2_001
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
    os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark --score_min L,0,-0.6 --genome /home/ankitv/ref_av/N-masked_hg38/Nhg38/ -1 ' + read1 + '_val_1.fq.gz -2 ' + read2 + '_val_2.fq.gz --bowtie2 --path_to_bowtie2 /home/ankitv/tools_av/bowtie2-2.3.4.1-linux-x86_64/')  

    print ('As the sample is amplicon duplicates are not removed')

    print ('Give a ID to BAM...')
    os.system('cp ' + read1 + '_val_1_bismark_bt2_pe.bam ' + folder + '_bismark_bt2_pe.bam')

    print ('Sorting data by readname...')
    os.system('samtools sort -n ' + folder + '_bismark_bt2_pe.bam -o ' + folder + '_bismark.sortedByReadname.bam')

    print ('Copy bam in new bam as per SNPfile... ')
    os.system('cp ' + folder + '_bismark.sortedByReadname.bam ' + folder + '_bismark.sortedByReadname_Acomb.bam')
    os.system('cp ' + folder + '_bismark.sortedByReadname.bam ' + folder + '_bismark.sortedByReadname_Bcomb.bam')

    print ('SNPsplit...')
    os.system('~/tools_av/SNPsplit_v0.3.2/SNPsplit --paired --bisulfite  --snp_file /home/ankitv/ref_av/Snpfiles/NNAT_commonSNPrefalt_GRCh38_Acomb_Snpfile.txt ' + folder + '_bismark.sortedByReadname_Acomb.bam')
    os.system('~/tools_av/SNPsplit_v0.3.2/SNPsplit --paired --bisulfite  --snp_file /home/ankitv/ref_av/Snpfiles/NNAT_commonSNPrefalt_GRCh38_Bcomb_Snpfile.txt ' + folder + '_bismark.sortedByReadname_Bcomb.bam')


    print ('BAM to fasta')
    os.system('samtools fasta -1 ' + folder + '_as1A_R1.fa -2 ' + folder + '_as1A_R2.fa ' + folder + '_bismark.sortedByReadname_Acomb.genome1.bam')
    os.system('samtools fasta -1 ' + folder + '_as2A_R1.fa -2 ' + folder + '_as2A_R2.fa ' + folder + '_bismark.sortedByReadname_Acomb.genome2.bam')
    os.system('samtools fasta -1 ' + folder + '_as1B_R1.fa -2 ' + folder + '_as1B_R2.fa ' + folder + '_bismark.sortedByReadname_Bcomb.genome1.bam')
    os.system('samtools fasta -1 ' + folder + '_as2B_R1.fa -2 ' + folder + '_as2B_R2.fa ' + folder + '_bismark.sortedByReadname_Bcomb.genome2.bam')

    print ('Sorting spliited bams by Coordinates...')
    os.system('samtools sort -o ' + folder + '_bismark.sortedByCoord_Acomb.genome1.bam ' + folder + '_bismark.sortedByReadname_Acomb.genome1.bam')
    os.system('samtools sort -o ' + folder + '_bismark.sortedByCoord_Acomb.genome2.bam ' + folder + '_bismark.sortedByReadname_Acomb.genome2.bam')
    os.system('samtools sort -o ' + folder + '_bismark.sortedByCoord_Bcomb.genome1.bam ' + folder + '_bismark.sortedByReadname_Bcomb.genome1.bam')
    os.system('samtools sort -o ' + folder + '_bismark.sortedByCoord_Bcomb.genome2.bam ' + folder + '_bismark.sortedByReadname_Bcomb.genome2.bam')

    print ('Index BAM')
    os.system('samtools index ' + folder +  '_bismark.sortedByCoord_Acomb.genome1.bam')
    os.system('samtools index ' + folder +  '_bismark.sortedByCoord_Acomb.genome2.bam')
    os.system('samtools index ' + folder +  '_bismark.sortedByCoord_Bcomb.genome1.bam')
    os.system('samtools index ' + folder +  '_bismark.sortedByCoord_Bcomb.genome2.bam')

    os.system('head ' + folder + '_as1A_R1.fa -n 8000 > ' + folder + '_as1A_R1_top.fa')
    os.system('head ' + folder + '_as2A_R1.fa -n 8000 > ' + folder + '_as2A_R1_top.fa')
    os.system('head ' + folder + '_as1B_R1.fa -n 8000 > ' + folder + '_as1B_R1_top.fa')
    os.system('head ' + folder + '_as2B_R1.fa -n 8000 > ' + folder + '_as2B_R1_top.fa')
    os.system('head ' + folder + '_as1A_R2.fa -n 8000 > ' + folder + '_as1A_R2_top.fa')
    os.system('head ' + folder + '_as2A_R2.fa -n 8000 > ' + folder + '_as2A_R2_top.fa')
    os.system('head ' + folder + '_as1B_R2.fa -n 8000 > ' + folder + '_as1B_R2_top.fa')
    os.system('head ' + folder + '_as2B_R2.fa -n 8000 > ' + folder + '_as2B_R2_top.fa')

    os.system('rm *.fq.gz')

    print ('Methylation extraction...')
    #os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark_methylation_extractor --bedGraph --include_overlap --zero_based --cutoff 10 ' + read1 + '_val_1_bismark.sortedByReadname.bam')

    print 'Analysis finished for ' + folder 
    print '\n'
    print '\n'
