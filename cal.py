import os
import subprocess

path=os.getcwd()+'/'
result=open("result.txt","a")
input_file=open("input.txt","r")
for line in input_file:
  line=line.split('_')
  metal=line[0]
  hfe=line[1]
  p=line[2]
  watermodel=line[3].strip('\n')
  workdir=path+'%s_%s_%s_%s'%(metal,hfe,p,watermodel)
  os.chdir(workdir)
  #print (os.getcwd())
  #os.system("sh check.sh")
  output = subprocess.getoutput("sh cal.sh")
  os.chdir(path)
  result.write("%s_%s_%s_%s	%s\n"%(metal,hfe,p,watermodel,output))
input_file.close()
result.close()
