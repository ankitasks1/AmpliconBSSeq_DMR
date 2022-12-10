import os,sys,glob

for i in glob.glob('*_as2B'):
    #print(i)
    os.system('convert -resize 1024x1024  ' + i + '/minus_ref_subNnat_R2/patternmap_Bisulfite.svg ' + i + '_minus_ref_subNnat_R2.png')


