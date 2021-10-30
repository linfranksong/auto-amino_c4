import os

path=os.getcwd()+'/'
input_file=open("%s/input.txt"%(path),"w")
dir1=path+'firstpass/'
os.chdir(dir1)
lis=os.listdir('.')
for direct in lis:
  if '_' in direct:
     metal,hfe,watermodel=direct.split('_')
     testout=open("%s%s/firstpass/test.out"%(dir1,direct),"r").readlines()
     pol1=round(float(testout[-1].split()[-1]) + 0.2,2)
     pol2=pol1+0.50
     input_file.write("%s_%s_%s_%s\n"%(metal,hfe,pol1,watermodel))
     input_file.write("%s_%s_%s_%s\n"%(metal,hfe,pol2,watermodel))
input_file.close()
os.chdir(path)
os.system("vim input.txt -c ':set nowrap'  -c ':%sort' -c 'x'")
