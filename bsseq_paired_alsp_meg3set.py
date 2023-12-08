import os, sys

list_data = ('''
S49 49_S49_L001_R1_001 49_S49_L001_R2_001
S50 50_S50_L001_R1_001 50_S50_L001_R2_001
S51 51_S51_L001_R1_001 51_S51_L001_R2_001
S52 52_S52_L001_R1_001 52_S52_L001_R2_001
S53 53_S53_L001_R1_001 53_S53_L001_R2_001
''').strip().split('\n')

print(list_data)

for listx in list_data:
    listx = listx.split(' ')
    folder = listx[0]
    read1 = listx[1]
    read2 = listx[2]
    print 'Sample in analysis: ' + folder
    print ('Trimming data: --phred33  (Sanger/Illumina 1.9+ encoding) for quality trimming. Default: ON, -q 20 (default),  BS-Seq' + folder)
    os.system('/home/ankitv/tools_av/TrimGalore-0.6.6/trim_galore  --phred33 --length 36 -q 20 --three_prime_clip_R2 60 --paired ' + read1 + '.fastq '+ read2 + '.fastq ')
    print ('Aligning data to N-masked using Bismark: ' + folder)
    os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark --score_min L,0,-0.6 --genome /home/ankitv/ref_av/N-masked_hg38/Nhg38/ -1 ' + read1 + '_val_1.fq -2 ' + read2 + '_val_2.fq --bowtie2 --path_to_bowtie2 /home/ankitv/tools_av/bowtie2-2.3.4.1-linux-x86_64/')  
    print ('As the sample is amplicon duplicates are not removed..')
    print ('Sorting data by readname...')
    os.system('samtools sort -n ' + read1 + '_val_1_bismark_bt2_pe.bam -o ' + read1 + '_val_1_bismark.sortedByReadname.bam')

    print ('Copy bam in new bam as per SNPfile... ')
    os.system('cp ' + read1 + '_val_1_bismark.sortedByReadname.bam ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.bam')
    os.system('cp ' + read1 + '_val_1_bismark.sortedByReadname.bam ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.bam')

    print ('SNPsplit...')
    os.system('~/tools_av/SNPsplit_v0.3.2/SNPsplit --paired --bisulfite  --snp_file /home/ankitv/ref_av/Snpfiles/custom_commonSNPrefalt_GRCh38_Acomb_Snpfile.txt ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.bam')
    os.system('~/tools_av/SNPsplit_v0.3.2/SNPsplit --paired --bisulfite  --snp_file /home/ankitv/ref_av/Snpfiles/custom_commonSNPrefalt_GRCh38_Bcomb_Snpfile.txt ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.bam')

    print ('Methylation extraction...')
    os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark_methylation_extractor --bedGraph --include_overlap --zero_based --cutoff 10 ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome1.bam')
    os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark_methylation_extractor --bedGraph --include_overlap --zero_based --cutoff 10 ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome2.bam')
    os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark_methylation_extractor --bedGraph --include_overlap --zero_based --cutoff 10 ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome1.bam')
    os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark_methylation_extractor --bedGraph --include_overlap --zero_based --cutoff 10 ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome2.bam')


    print ('Sorting data...')

    os.system('samtools sort -o ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome1.sort.bam ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome1.bam')
    os.system('samtools sort -o ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome2.sort.bam ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome2.bam')
    os.system('samtools sort -o ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome1.sort.bam ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome1.bam')
    os.system('samtools sort -o ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome2.sort.bam ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome2.bam')

    print ('Index data...')
    os.system('samtools index ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome1.sort.bam')
    os.system('samtools index ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome1.sort.bam')
    os.system('samtools index ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome2.sort.bam')
    os.system('samtools index ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome2.sort.bam')

    os.system('gunzip *.bismark.cov.gz')

    os.system('paste  ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome1.bismark.cov  ' + read1 + '_val_1_bismark.sortedByReadname_Acomb.genome2.bismark.cov  > '  + read1 + '_Acomb.coverage.txt')
    os.system('paste  ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome1.bismark.cov  ' + read1 + '_val_1_bismark.sortedByReadname_Bcomb.genome2.bismark.cov  > '  + read1 + '_Bcomb.coverage.txt')


    os.system('rm C*')

    print 'Analysis finished for ' + folder 
    print '\n'
    print '\n'

