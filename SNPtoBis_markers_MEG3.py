import os,sys

list_data = ('''
S49_sam_Reads.txt
S50_sam_Reads.txt
S51_sam_Reads.txt
S52_sam_Reads.txt
S53_sam_Reads.txt
''').strip().split('\n')

#print(''.join(str(list_data)+ '_out.txt'))
#print(type('test_SNP.txt'))
#print("\n")

for myfile in list_data:
  myfile1= open(myfile ,'r')
  #myfile1.close
  file_name = myfile
  print('Analysis started for' + file_name)
  #print(myfile1.read())
  for row in myfile1:
    row = row.split(' ')
    sequence = list(row[3])
    #print(row)
    SNPfile = open("mySNPlist.txt" ,'r')
    SNPfile.close
    #print(SNPfile.read())
    Methfile = open("myImpLoci_CGs.txt" ,'r')
    Methfile.close
    #print(Methfile.read())
    print('Reading input files: mySNPlist.txt, myImpLoci_CGs.txt, '+file_name)
    for snp in SNPfile:
        snp = snp.split(' ')
        if row[0] == snp[0]:
           Refsite = snp[3]
           Altsite = snp[4]
           RefsiteID = snp[5].strip()
           AltsiteID = snp[5].strip()
           #print(snp)
           SNPstart = snp[1]
           readstart = row[1]
           readend = row[2]
           #print(SNPstart + "\n" +readstart)
           if int(SNPstart) < int(readend) and int(SNPstart) > int(readstart):
              #print(SNPstart+ "\t" + readstart + "\t" + readend)
              #print(row[4])
              #print("Refsite position =  "+SNPstart)
              #print("Refsite is =  "+Refsite)
              #print("Altsite is =  "+Altsite)

              #print(int(SNPstart))
              #print(int(row[1]))
              readsnpdifference = int(SNPstart) - int(readstart)
              #print(readsnpdifference)
              #print("Base in sequence at site is =  "+ sequence[readsnpdifference -3]+ sequence[readsnpdifference -2]+ "[" +sequence[readsnpdifference -1]+ "]" +sequence[readsnpdifference -0]+ sequence[readsnpdifference +1])

              if sequence[readsnpdifference-1] == Refsite:
                #print(Refsite, sequence[readsnpdifference-1])
                #print(row)
                for meth in Methfile:
                   meth = meth.split(' ')
                   if row[0] == meth[0]:
                      Methstart = meth[1]
                      Methsite = meth[3]
                      MethsiteID = meth[5].strip()
                      if int(Methstart) > int(readstart) and int(Methstart) < int(readend):
                         readmethdifference = int(Methstart) - int(readstart)
                         #print(''.join(row[4] + "\t"+ MethsiteID + "\t"+ RefsiteID + "\t"+ Refsite+ "\t" + Methsite + "\t" + sequence[readmethdifference-1]))
                         inforeadref = ''.join(row[4] + "\t"+ MethsiteID + "\t"+ RefsiteID + "\t"+ Refsite + "\t" + Methsite + "\t" + sequence[readmethdifference-1])
                         #print(inforeadref)
                         file1 = open(file_name+'_'+RefsiteID+'_out.txt','a')
                         file1.write(inforeadref +'\n')
                         file1.close()

              elif sequence[readsnpdifference-1] == Altsite:
                 #print(Altsite, sequence[readsnpdifference])
                 for meth in Methfile:
                     meth = meth.split(' ')
                     if row[0] == meth[0]:
                        Methstart = meth[1]
                        Methsite = meth[3]
                        MethsiteID = meth[5].strip()
                        if int(Methstart) > int(readstart) and int(Methstart) < int(readend):
                           readmethdifference = int(Methstart) - int(readstart)
                           #print(''.join(row[4] + "\t"+ MethsiteID + "\t"+ AltsiteID + "\t"+ Altsite + "\t" + Methsite+ "\t" + sequence[readmethdifference]))
                           inforeadalt = ''.join(row[4] + "\t"+ MethsiteID + "\t"+ AltsiteID + "\t"+ Altsite + "\t" + Methsite + "\t" + sequence[readmethdifference-1])
                           #print(inforeadalt)
                           file1 = open(file_name+'_'+RefsiteID+'_out.txt','a')
                           file1.write(inforeadalt +'\n')
                           file1.close()
   
              else:
               inforeadnone = ''.join(row[4] + "\t"+ "NA"  + "\t"+ "NA"+ "\t"+ "NA" + "\t" + "NA"+ "\t" + "NA")
               #print(''.join(row[4] + "\t"+ "NA"  + "\t"+ "NA" + "\t"+ "NA" + "\t" + "NA"+ "\t" + "NA"))
               #print(inforeadnone)

               file1 = open(file_name+'_'+RefsiteID+'_out.txt','a')
               file1.write(inforeadnone +'\n')
               file1.close()
print('Analysis finished for '+file_name+'\n')

