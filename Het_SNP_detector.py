import os,sys,glob

#Call variants
#unzip BS-Snper-master.zip
#cd BS-Snper-master
#sh BS-Snper.sh
#run call_variants.py


#print(os.listdir('./'))
#print(glob.glob('*_snp.list'))

for i in glob.glob('*_snp.list'):
    sampleid = i.strip().split('_')
    #print(sampleid)
    myfile = open(i , 'r')
    myfile = myfile.read().strip().split('\n')
    del myfile[0]
    #print(myfile)
    for j in myfile:
        j = j.strip().split('\t')
        if int(j[1]) == 37520583:
           #print(j)
           readcoverage = j[4]
           readcoverage = readcoverage.split(',')
           proportion = int(readcoverage[1]) / (int(readcoverage[1]) + int(readcoverage[2]))
           if float(proportion) >= 0.3 and float(proportion) <= 0.7:
              print(i,"\t" , int(readcoverage[1]),"\t" , int(readcoverage[2]),"\t" , proportion,"\t" , "Heterozygous")
           else:
              print(i,"\t" , int(readcoverage[1]),"\t"  , int(readcoverage[2]),"\t" , proportion,"\t" , "Homozygous")

