import os

dir = os.getcwd()+'/'

lamb=[0,0.0025,0.005,0.0075,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,1.0]
for n in range(0,len(lamb)):
  i=lamb[n]
  f=float(AAA)
  m=round(f*i,2)
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
  os.system("sbatch mh.pbs")
  os.chdir(dir)

