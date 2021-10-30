import os


import numpy as np
dir = os.getcwd()+'/'

lamb = [0,0.0025,0.005,0.0075,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,1.0,1.0]
for n in range(0,len(lamb)-1):
  i=lamb[n]
  f=float(AAA)
  m=round(f*i,2)
  h=lamb[n-1]
  j=lamb[n+1]
  os.system("rm -r %s"%(i))
  os.system("mkdir %s"%(i))
  wdir=dir+"%s/"%(i)
  os.chdir(wdir)
  os.system("cp ../temp/dis.RST .")
  os.system("cp ../temp/*1264* .")
  os.system("cp ../temp/*in .")
  os.system("cp ../temp/*pbs .")
  os.system("sed -i 's/NNN/%s/g' *.RST"%(m))
  os.system("sed -i 's/XXX/%s/g' *.pbs"%(i))
  os.system("sed -i 's/PPP/%s/g' *.pbs"%(h))
  os.system("sed -i 's/NNN/%s/g' *.pbs"%(j))
  os.chdir(dir)
os.system("cp temp/1_eq.pbs 1.0/eq.pbs")
dir2=dir+'0/'
os.chdir(dir2)
os.system("mv 0_eq.pbs eq.pbs")
os.system("sbatch mh.pbs")
os.chdir(dir)
