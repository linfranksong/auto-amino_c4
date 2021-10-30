#!/mnt/home/songlin3/anaconda3/bin/python
from __future__ import print_function
import os
import numpy as np
import math

def pmf(try_pol,metal,gaff,protein_ff,water_model,water_box,amber_frcmod_files,additional_frcmod_files,additional_frcmod_num,amber_frcmod_num, amber_mol2_files,amber_mol2_num,additional_mol2_files,additional_mol2_num,Na_ion,Cl_ion,cyx_num,cyx_bond,box_size,HZ_resid,coord_pol,coord_rmin,coord_epsilon,coord_atomid1,coord_atomid2,coord_atomid3,dis,ang,dihe,m_dihe,p_dihe,disrst_strength,min_steps,heat_steps,equi_steps,sampl_steps,slm_header,cutoff,time_step):
  dir1=os.getcwd()+'/'
  os.chdir(dir1)
  if os.path.isdir('pmf/'):
    dir2=dir1+'pmf/'
    os.chdir(dir2)
    os.system("cp ../freemd/126.* .")
    os.system("cp ../freemd/addc4.in .")
    os.system("cp ../freemd/%s.txt ."%(coord_pol))
    #os.system("cp ~/auto-amino_c4/sample.pbs .")
    os.system("sed -i 's/MUTE  %s/MUTE  %s/g' %s.txt"%(coord_pol,try_pol, coord_pol))
    os.system("mv %s.txt %s.txt"%(coord_pol,try_pol))
    os.system("sed -i 's/polfile %s.txt/polfile %s.txt/g' addc4.in"%(coord_pol,try_pol))
    os.system("parmed -i addc4.in")
    os.system("sed -i 's/us_md_steps  4000000/us_md_steps  1000000/g' pmf.in")
    with open("indicator.txt", "r") as indicator:
      line=indicator.readlines()
      if "%s us1 submitted"%(try_pol) not in line[-1]:
        if "%s us2 submitted"%(try_pol) not in line[-1]:
          #os.system("sbatch sample.pbs")
          print ("running python check_prepare_md.py")
          os.system("cp ~/auto-amino_c4/pmf/check_prepare_md.py .")
          os.system("python check_prepare_md.py")
          arg1=input("!!!! if no error reported by the check_prepare_md.py, please type 'yes' to run us1, otherwise type'no' to stop\n")
          if arg1=='yes':
            os.system("Auto-PMF-NVT-di_equi-cpu.py -i pmf.in -s sampling_NVT")
            indicator=open("indicator.txt","a")
            indicator.seek(0,2)
            indicator.write("%s us1 submitted\n"%(try_pol))
            print("!!! running us1 for %s !!!"%(try_pol))
          else:
            pass
        elif "%s us2 submitted"%(try_pol) in line[-1]:
          print ("!!! %s us2 submitted most recently, after jobs finished, please run check_force2, adaptive_force2, etc. for %s !!!"%(try_pol,try_pol))
      elif "%s us1 submitted"%(try_pol) in line[-1]:
        print ("!!!! Attention: %s us1 submitted most recently, before processing to running us2 for %s, please make sure 'check_force, adaptive_force, write_wham, check_wham, result' steps are completed for us1, and the generated 2ns/window pmf gives a binding affinity close to experiment. !!!!"%(try_pol,try_pol))
        arg2=input("Type 'yes' to run us2, 'no' to stop\n")
        if arg2=='yes':
          os.system("Auto-PMF-NVT-di_equi-cpu.py -i pmf.in -s us2_NVT")
          indicator=open("indicator.txt","a")
          indicator.seek(0,2)
          indicator.write("%s us2 submitted\n"%(try_pol))
        else:
          pass
  else:
    os.system("mkdir pmf")
    dir2=dir1+'pmf/'
    os.chdir(dir2)
    os.system("cp ~/auto-amino_c4/pmf/* .")
    os.system("cp ../freemd/126.* .")
    os.system("cp ../freemd/addc4.in .")
    os.system("cp ../freemd/%s.txt ."%(coord_pol))
    os.system("sed -i 's/MUTE  %s/MUTE  %s/g' %s.txt"%(coord_pol,try_pol, coord_pol))
    os.system("mv %s.txt %s.txt"%(coord_pol,try_pol))
    os.system("sed -i 's/polfile %s.txt/polfile %s.txt/g' addc4.in"%(coord_pol,try_pol))
    os.system("parmed -i addc4.in")
    pdb5=open("126.pdb","r")
    for lines in pdb5:
      line=lines.split()
      if line[0]!= "TER" and line[0]!= "END" and line[3]=='%s'%(metal) and line[4]=='%s'%(HZ_resid):
        HZ_atomid1=int(line[1])
    pdb5.close()
    os.system("sed -i 's/HZ_ATOMID1/%s/g' pmf.in"%(HZ_atomid1))
    os.system("sed -i 's/COORD_ATOMID1/%s/g' pmf.in"%(coord_atomid1))
    print ("!!! running prepare for %s !!!"%(try_pol))
    os.system("sbatch run.pbs")
    indicator=open("indicator.txt","a")
    indicator.seek(0,2)
    indicator.write("%s prepare submitted\n"%(try_pol))
