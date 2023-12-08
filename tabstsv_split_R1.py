import os,sys

tabtsv = open('table_results.tsv', 'r')
tabtsv = tabtsv.read()
tabtsv = tabtsv.strip().split('\n')
#print(tabtsv)
for i in tabtsv:
    if len(i) > 80:
        #print(len(i))
        i = i.strip().split('\t')
        #print(i)
        methlist = list(i[3])
        methlist = "\t".join(methlist)
        print(i[0],i[8],i[3], methlist)
