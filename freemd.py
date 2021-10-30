#!/mnt/home/songlin3/anaconda3/bin/python
from __future__ import print_function
import os
import numpy as np
import math

def freemd(metal,gaff,protein_ff,water_model,water_box,amber_frcmod_files,additional_frcmod_files,additional_frcmod_num,amber_frcmod_num, amber_mol2_files,amber_mol2_num,additional_mol2_files,additional_mol2_num,Na_ion,Cl_ion,cyx_num,cyx_bond,box_size,HZ_resid,coord_pol,coord_rmin,coord_epsilon,coord_atomid1,coord_atomid2,coord_atomid3,dis,ang,dihe,m_dihe,p_dihe,disrst_strength,min_steps,heat_steps,equi_steps,sampl_steps,slm_header,cutoff,time_step):
  dir1=os.getcwd()+'/'
  os.chdir(dir1)
  os.system("rm -r freemd")
  os.system("mkdir freemd")
  dir2=dir1+'freemd/'
  os.chdir(dir2)

  additional_frcmod=[]
  for num in range(0,additional_frcmod_num):
    os.system("cp %s ."%(additional_frcmod_files[num]))
    additional_frcmod.append(additional_frcmod_files[num].split('/')[-1])
  additional_mol2=[]
  for num in range(0,additional_mol2_num):
    os.system("cp %s ."%(additional_mol2_files[num]))
    additional_mol2.append(additional_mol2_files[num].split('/')[-1])
  os.system("cp ~/auto-amino_c4/source/temp.dat .")
  os.system("sed -i 's/PPPPP/%s/g' temp.dat"%(coord_pol))
  os.system("cp ../input.pdb .")



################### write tleap file ###################
  tleap=open("tleap.in","w")

  print("source leaprc.protein.%s"%(protein_ff),file=tleap)
  print("source leaprc.%s"%(gaff),file=tleap)
  print("source leaprc.water.%s"%(water_model.lower()),file=tleap)
  for num in range(0,amber_frcmod_num):
    print("loadamberparams %s"%(amber_frcmod_files[num]),file=tleap)
  for num in range(0,amber_mol2_num):
    ligand,mol2=amber_mol2_files[num].split(".mol2")
    print("%s=loadmol2 %s"%(ligand,amber_mol2_files[num]),file=tleap)
  for num in range(0,len(additional_frcmod)):
    print("loadamberparams %s"%(additional_frcmod[num]),file=tleap)
  for num in range(0,len(additional_mol2)):
    ligand,mol2=additional_mol2[num].split(".mol2")
    print("%s=loadmol2 %s"%(ligand,additional_mol2[num]),file=tleap)

  print("mol = loadpdb input.pdb",file=tleap)
  for num in range(0,cyx_num):
    num=num*2
    print("bond mol.%s.SG mol.%s.SG"%(cyx_bond[num],cyx_bond[num+1]),file=tleap)

  print("solvateoct mol %sBOX %s"%(water_box, box_size),file=tleap)
  #print("addions2 mol Na+ %s"%(Na_ion),file=tleap)
  #print("addions2 mol Cl- %s"%(Cl_ion),file=tleap)
  print("saveamberparm mol 126.prmtop 126.inpcrd",file=tleap)
  print("quit",file=tleap)

  tleap.close()

  os.system("tleap -s -f tleap.in > tleap.out")
  os.system("ambpdb -p 126.prmtop < 126.inpcrd > 126.pdb")

  addc4=open("addc4.in","w")
  print("parm 126.prmtop\nloadRestrt 126.inpcrd\nsetOverwrite True\n",file=addc4) 
  print("change AMBER_ATOM_TYPE @%s MUTE\naddLJType @%%MUTE\n"%(coord_atomid1),file=addc4)

  print("add12_6_4 :2 watermodel %s polfile MMM.txt \noutparm 1264.prmtop 1264.inpcrd\nprintLJMatrix :2"%(water_model.lower()),file=addc4)
  addc4.close()

  os.system("cp ~/auto-amino_c4/freemd/* .")

####################################### write run.sh cal.sh check.sh #################################

  runsh=open("run.sh","w")
  print('''sed -i "s/XXX/%s_%s/g" *pbs\nsed -i "s/MMM/%s/g" temp.dat\nsed -i "s/MMM/%s/g" addc4.in\nmv temp.dat %s.txt\nparmed -i addc4.in\nsbatch mh.pbs\n'''%(coord_pol,metal,coord_pol,coord_pol,coord_pol),file=runsh)
  runsh.close()

  ptrajin=open("ptraj.in","w")
  print("trajin us.netcdf\nautoimage\ndistance distance :%s @%s out bond.txt\nangle angle :%s @%s @%s out angle.txt\ndihedral dihedral :%s @%s @%s @%s out dihedral.txt\n"%(HZ_resid,coord_atomid1,HZ_resid,coord_atomid1,coord_atomid2,HZ_resid,coord_atomid1,coord_atomid2,coord_atomid3),file=ptrajin)
  ptrajin.close()

  calsh=open("cal.sh","w")
  print("awk 'NR>1 {sum+=$2} END {print sum/(FNR-1)}' bond.txt\nawk 'NR>1 {su+=$2} END {print su/(FNR-1)}' angle.txt\nawk 'NR>1 {s+=$2} END {print s/(FNR-1)}' dihedral.txt\n",file=calsh)
  calsh.close()

  checksh=open("check.sh","w")
  print("echo 'Attention: blew up'\nfgrep -l '*******' *out\n",file=checksh)
  checksh.close()
  
  os.system("sh run.sh")
