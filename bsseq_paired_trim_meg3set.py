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

    print ('Trimming data..: --phred33  (Sanger/Illumina 1.9+ encoding) for quality trimming. Default: ON, -q 20 (default),  BS-Seq' + folder)
    os.system('/home/ankitv/tools_av/TrimGalore-0.6.6/trim_galore  --phred33 --length 36 -q 20 --three_prime_clip_R2 60 --paired ' + read1 + '.fastq '+ read2 + '.fastq ')

    print ('Aligning data using Bismark: ' + folder)
    os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark --score_min L,0,-0.6 --genome /home/ankitv/ref_av/hg38_bs/ -1 ' + read1 + '_val_1.fq -2 ' + read2 + '_val_2.fq --bowtie2 --path_to_bowtie2 /home/ankitv/tools_av/bowtie2-2.3.4.1-linux-x86_64/')  

    print ('As the sample is amplicon duplicates are not removed')

    print ('Sorting data...')
    os.system('samtools sort -o ' + read1 + '_val_1_bismark_bt2_pe.sort.bam ' + read1 + '_val_1_bismark_bt2_pe.bam')

    print ('Index data...')
    os.system('samtools index ' + read1 + '*_val_1_bismark_bt2_pe.sort.bam')

    print ('Sorting data by readname...')
    os.system('samtools sort -n ' + read1 + '_val_1_bismark_bt2_pe.bam -o ' + read1 + '_val_1_bismark.sortedByReadname.bam')

    print ('Methylation extraction...')
    #os.system('/home/ankitv/tools_av/Bismark-0.22.3/bismark_methylation_extractor --bedGraph --include_overlap --zero_based --cutoff 10 ' + read1 + '_val_1_bismark.sortedByReadname.bam')

    print 'Analysis finished for ' + folder 
    print '\n'
    print '\n'
