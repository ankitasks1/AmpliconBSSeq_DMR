import os,sys
if os.path.exists('CGSites.txt'):
   os.remove('CGSites.txt')
#R2 is in reverse so coordinates will be in reverse direction
minus_ref_subNnat_R2 = open(input('Enter fasta file: '), 'r')
readtype = input('Enter readtype: ')
minus_ref_subNnat_R2 = minus_ref_subNnat_R2.read()
minus_ref_subNnat_R2 = minus_ref_subNnat_R2.strip().split('\n')
#print(minus_ref_subNnat_R2)

header = minus_ref_subNnat_R2[0]
header = header.split(' ')
#print(header)
coordinates = header[1]


start  = coordinates.strip().split(':')
start  = start[1].strip().split('-')
#print(start)
sequence = minus_ref_subNnat_R2[1]
#print(sequence)
sequence = list(sequence)
SeqStart = 0
count = 0
for i in range(0,len(sequence)):
    SeqStart += 1
    #print(SeqStart)
    twobase = sequence[SeqStart: SeqStart + 2]
    twobase = ''.join(twobase)
    #print(twobase)

    if twobase == 'CG':
        count += 1
        CGsites = "".join(twobase + str(count)+ readtype)
        CGStart = "".join(str(int(start[1]) - int(SeqStart)))
        CGsiteStart = "".join(CGsites + '\t' + CGStart)
        print(CGsiteStart)
        with open('CGSites.txt', 'a') as myfile:
            myfile.write(CGsiteStart +'\n')
            myfile.close()
    else:
        pass

