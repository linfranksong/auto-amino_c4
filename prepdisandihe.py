#!/mnt/home/songlin3/anaconda3/bin/python
from __future__ import print_function
import os
import subprocess as sbs
from optparse import OptionParser
import numpy as np
import math

def prepdisandihe(metal,hfe,gaff,protein_ff,water_model,water_box,amber_frcmod_files,additional_frcmod_files,additional_frcmod_num,amber_frcmod_num, amber_mol2_files,amber_mol2_num,additional_mol2_files,additional_mol2_num,Na_ion,Cl_ion,cyx_num,cyx_bond,box_size,HZ_resid,coord_pol,coord_rmin,coord_epsilon,coord_atomid1,coord_atomid2,coord_atomid3,dis,ang,dihe,m_dihe,p_dihe,disrst_strength,min_steps,heat_steps,equi_steps,sampl_steps,slm_header,cutoff,time_step):
  dir1=os.getcwd()+'/'
  os.system("rm -r temp")
  os.system("mkdir temp")
  dir2=dir1+'temp/'
  os.chdir(dir2)
  os.system("mkdir step0_turnoff step1 step2")
  dir0=dir2+'step0_turnoff/'
  os.chdir(dir0)
  os.system("mkdir temp/")
  dir00=dir0+'temp/'
  dir3=dir2+'step1/'
  os.chdir(dir3)
  os.system("mkdir temp/")
  dir4=dir3+'temp/'
  dir5=dir2+'step2/'
  os.chdir(dir5)
  os.system("mkdir temp/")
  dir6=dir5+'temp/'
  os.chdir(dir00)


  os.system("cp ../../../input.pdb .")
  additional_frcmod=[]
  for num in range(0,additional_frcmod_num):
    os.system("cp %s ."%(additional_frcmod_files[num]))
    additional_frcmod.append(additional_frcmod_files[num].split('/')[-1])
  additional_mol2=[]
  for num in range(0,additional_mol2_num):
    os.system("cp %s ."%(additional_mol2_files[num]))
    additional_mol2.append(additional_mol2_files[num].split('/')[-1])
  os.system("cp ~/auto-amino_c4/source/* .")
  os.system("mv c.py ..")
  os.system("sed -i 's/PPPPP/%s/g' temp.dat"%(coord_pol))




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

  HZ_resid1=HZ_resid
  HZ_resid2=HZ_resid1 + 1
  os.system("cp input.pdb input.pdb-f")
  pdb4=open("input.pdb","r")
  pdb3=open("input.pdb-f","a")
  for lines in pdb4:
    line = lines.split()
    if line[0] != "TER"  and line[0] != "END" and line[3] == '%s'%(metal) and line[4] == '%s'%(HZ_resid1):
      new_lines=lines.replace("%s     %s"%(metal,HZ_resid1),"%s     %s"%(metal,HZ_resid2))
      pdb3.write(new_lines)
  pdb4.close()
  pdb3.close()

  os.system("tleap -s -f tleap.in > tleap.out")
  os.system("ambpdb -p 126.prmtop < 126.inpcrd > 126.pdb")

  addc4=open("addc4.in","w")
  print("parm 126.prmtop\nloadRestrt 126.inpcrd\nsetOverwrite True\n",file=addc4) 
  print("change AMBER_ATOM_TYPE @%s MUTE\naddLJType @%%MUTE\n"%(coord_atomid1),file=addc4)

  print("add12_6_4 :2 watermodel %s polfile MMM.txt\noutparm 1264.prmtop 1264.inpcrd\nprintLJMatrix :2"%(water_model.lower()),file=addc4)
  addc4.close()

  #print (metal)

  #print (HZ_resid1)

  pdb5=open("126.pdb","r")
  for lines in pdb5:
    line=lines.split()
    if line[0]!= "TER" and line[0]!= "END" and line[3]=='%s'%(metal) and line[4]=='%s'%(HZ_resid1):
      HZ_atomid1=int(line[1])
      HZ_atomid2=HZ_atomid1+1
  pdb5.close()


  disrst=open("dis.RST","w")
  print("# restraint",file=disrst)
  print(" &rst  iat=%s, %s, r1=0.0, r2=%s, r3=%s, r4=100., rk2=NNN, rk3=NNN,/"%(HZ_atomid1,coord_atomid1,dis,dis),file=disrst)
  print(" &rst  iat=%s, %s, %s, r1=0.0, r2=%s, r3=%s, r4=180., rk2=NNN, rk3=NNN,/"%(HZ_atomid1,coord_atomid1,coord_atomid2,ang,ang),file=disrst)
  print(" &rst  iat=%s, %s, %s, %s, r1=%s, r2=%s, r3=%s, r4=%s, rk2=NNN, rk3=NNN,/"%(HZ_atomid1,coord_atomid1,coord_atomid2,coord_atomid3,m_dihe,dihe,dihe,p_dihe),file=disrst)
  disrst.close()
####################################### copy ti input files and pbs #################################

  os.system("cp ~/auto-amino_c4/step0_turnoff/* .")

  os.chdir(dir0)
  os.system("mv temp/set.py .")

  os.system("sed -i 's/DDDDD/%s/g' c*.py"%(dis))
  os.system("sed -i 's/AAAAA/%s/g' c*.py"%(ang))
  os.system("sed -i 's/HHHHH/%s/g' c*.py"%(dihe))

  os.chdir(dir2)


  os.system("cp -r step0_turnoff step0_para")
  os.system("cp ~/auto-amino_c4/step0_para/* step0_para/temp/")
  os.system("mv step0_para/temp/set.py step0_para/")

####################################### build folder of step1 ########################################
  os.chdir(dir4)

  os.system("cp ~/auto-amino_c4/source/* .")
  os.system("cp ../../step0_turnoff/temp/*mol2 .")
  os.system("cp ../../step0_turnoff/temp/*frcmod .")
  os.system("sed -i 's/PPPPP/%s/g' temp.dat"%(coord_pol))

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

  print("mol = loadpdb input.pdb-f",file=tleap)
  for num in range(0,cyx_num):
    num=num*2
    print("bond mol.%s.SG mol.%s.SG"%(cyx_bond[num],cyx_bond[num+1]),file=tleap)

  print("solvateoct mol %sBOX %s"%(water_box, box_size),file=tleap)
  #print("addions2 mol Na+ %s"%(Na_ion),file=tleap)
  #print("addions2 mol Cl- %s"%(Cl_ion),file=tleap)
  print("saveamberparm mol 126.prmtop 126.inpcrd",file=tleap)
  print("quit",file=tleap)

  tleap.close()

  os.system("mv ../../step0_turnoff/temp/input.pdb-f .")
  os.system("tleap -s -f tleap.in > tleap.out")

  addc4=open("addc4.in","w")
  print("parm 126.prmtop\nloadRestrt 126.inpcrd\nsetOverwrite True\n",file=addc4)
  print("change AMBER_ATOM_TYPE @%s MUTE\naddLJType @%%MUTE\n"%(coord_atomid1),file=addc4)
  print("add12_6_4 :2-3 watermodel %s polfile MMM.txt\noutparm 1264.prmtop 1264.inpcrd\nprintLJMatrix :2-3"%(water_model.lower()),file=addc4)
  addc4.close()

  os.system("cp ../../step0_turnoff/temp/dis.RST .")
  os.system("sed -i 's/NNN/MMM/g' dis.RST")
  disrst=open("dis.RST","a")
  disrst.seek(0,2)
  print(" &rst  iat=%s, %s, r1=0.0, r2=%s, r3=%s, r4=100., rk2=MMM, rk3=MMM,/"%(HZ_atomid2,coord_atomid1,dis,dis),file=disrst)
  print(" &rst  iat=%s, %s, %s, r1=0.0, r2=%s, r3=%s, r4=180., rk2=MMM, rk3=MMM,/"%(HZ_atomid2,coord_atomid1,coord_atomid2,ang,ang),file=disrst)
  print(" &rst  iat=%s, %s, %s, %s, r1=%s, r2=%s, r3=%s, r4=%s, rk2=MMM, rk3=MMM,/"%(HZ_atomid2,coord_atomid1,coord_atomid2,coord_atomid3,m_dihe,dihe,dihe,p_dihe),file=disrst)
  disrst.close()

  os.system("cp ~/auto-amino_c4/step1/* .")
  os.system("sed -i 's/:AAA/:%s/g' *in"%(HZ_resid1))
  os.system("sed -i 's/:BBB/:%s/g' *in"%(HZ_resid2))

  os.chdir(dir3)
  setpy=open("set.py","w")
  print('''import os\n\ndir = os.getcwd()+'/'\n\nlamb = [0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738, 0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078]\nfor n in range(0,len(lamb)):\n  i=lamb[n]\n  os.system("rm -r %s"%(i))\n  os.system("mkdir %s"%(i))\n  wdir=dir+"%s/"%(i)\n  os.chdir(wdir)\n  os.system("cp ../temp/dis.RST .")\n  os.system("cp ../temp/*1264* .")\n  os.system("cp ../temp/*in .")\n  os.system("cp ../temp/*pbs .")\n  os.system("sed -i 's/XXX/%s/g' *.in"%(i))\n  os.system("sed -i 's/XXX/%s/g' *.pbs"%(i))\n  os.system("sbatch mh.pbs")\n  os.chdir(dir)\n''',file=setpy)
  setpy.close()



####################################### build folder step2 #################################
  os.system("cp set.py %s"%(dir5))
  os.system("cp temp/126.* %s"%(dir6))
  os.system("cp temp/dis.RST %s"%(dir6))
  os.system("cp temp/temp.dat %s"%(dir6))
  os.chdir(dir6)

  addc4=open("addc4.in","w")
  print("parm 126.prmtop\nloadRestrt 126.inpcrd\nsetOverwrite True\n",file=addc4) 
  print("change AMBER_ATOM_TYPE @%s MUTE\naddLJType @%%MUTE\nchange AMBER_ATOM_TYPE :3 hcg\naddLJType @%%hcg\nchangeLJSingleType :3@%%hcg 0 0\nchange charge :3 0"%(coord_atomid1),file=addc4)
  print("add12_6_4 :2 watermodel %s polfile MMM.txt\noutparm 1264_step2.prmtop 1264_step2.inpcrd\nprintLJMatrix :2\nprintLJMatrix :3\nnetcharge :2\nnetcharge :3\n"%(water_model.lower()),file=addc4)
  addc4.close()

  os.system("cp ~/auto-amino_c4/step2/* .")
  os.system("sed -i 's/:AAA/:%s/g' *in"%(HZ_resid1))
  os.system("sed -i 's/:BBB/:%s/g' *in"%(HZ_resid2))
  
  os.system("cp -r ~/auto-amino_c4/calcu/ %s"%(dir3))
  os.system("cp -r ~/auto-amino_c4/calcu/ %s"%(dir5))

####################################### write run.sh cal.sh check.sh #################################

  os.chdir(dir1)  
  runsh=open("run.sh","w")
  runsh.write("for i in ")
  for num in range(len(disrst_strength)):
    disrst=disrst_strength[num]
    runsh.write("%s "%(disrst))
  print('''\ndo\nrm -r %s_$i\ncp -r temp %s_$i\ncd %s_$i/\nsed -i "s/MMM/%s/g" step*/temp/temp.dat\nsed -i "s/MMM/%s/g" step*/temp/addc4.in\nsed -i "s/MMM/$i/g" step*/temp/dis.RST\nsed -i "s/AAA/$i/g" step*/*.py\ncd step1/temp\nmv temp.dat %s.txt\nparmed -i addc4.in\ncd ../../step0_turnoff/temp/\nmv temp.dat %s.txt\nparmed -i addc4.in\ncd ../../step2/temp/\nmv temp.dat %s.txt\nparmed -i addc4.in\ncd ../../step0_para/\ncp ../step0_turnoff/temp/1264* temp/\npython set.py\ncd ../step1/\npython set.py\ncd ../step2/\npython set.py\ncd ../../\ndone'''%(coord_pol,coord_pol,coord_pol,coord_pol,coord_pol,coord_pol,coord_pol,coord_pol),file=runsh)
  runsh.close()

  hfe=abs(float(hfe))
  calsh=open("cal.sh","w")
  calsh.write("for i in ")
  for num in range(len(disrst_strength)):
    disrst=disrst_strength[num]
    calsh.write("%s "%(disrst))
  print('''\ndo\ncd %s_$i/step0_para\npython c.py > deltaG.txt\ng0=$(head -n 1 deltaG.txt)\ncd ../step1/calcu/\nsh c.sh\ng1=$(head -n 1 deltaG.txt)\ncd ../../step2/calcu/\nsh c.sh\ng2=$(head -n 1 deltaG.txt)\ncd ../../../\ng3=$(echo "-0.596161*l(e(1.5*l(%s*2))*1661/(e(1.5*l(0.0019872041*300*2*3.141592653589))*(%s*%s)*s(3.141592653589*%s/180)))" | bc -l)\ngg0=$(echo "e(0.5*l($g0*$g0))"| bc -l)\ngg1=$(echo "e(0.5*l($g1*$g1))"| bc -l)\ngg2=$(echo "e(0.5*l($g2*$g2))"| bc -l)\ngg3=$(echo "e(0.5*l($g3*$g3))"| bc -l)\nghfe=$(echo "e(0.5*l(%.2f*%.2f))"| bc -l)\necho "$ghfe-$gg0-$gg1-$gg2+$gg3" | bc -l\ndone'''%(coord_pol,disrst,dis,dis,ang,float(hfe),float(hfe)),file=calsh)
  calsh.close()

  dis_bound=float(float(dis)+0.3)
  checksh=open("check.sh","w")
  checksh.write("for i in ")
  for num in range(len(disrst_strength)):
    disrst=disrst_strength[num]
    checksh.write("%s "%(disrst))
  print('''\ndo\ncd %s_$i/step0_para/\necho 'step0_para'\nfgrep -l '*******' */*out\ncd ../step1/\necho 'step1'\nfgrep -l '*******' 0.*/*out\ncd ../step2/\necho 'step2'\nfgrep -l '*******' 0.*/*out\ncd ../../\ndone\n\n\n'''%(coord_pol),file=checksh)
  checksh.write("for i in ")
  for num in range(len(disrst_strength)):
    disrst=disrst_strength[num]
    checksh.write("%s "%(disrst))
  print('''\ndo\ncd %s_$i/step0_para/\necho 'step0_para'\nfor window in 0 0.0025 0.005 0.0075 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.4 0.6 0.8 1.0\ndo\ncd $window\navg=$(awk '{sum=sum+$2} END {print sum/FNR}' us.log)\nif (( $(echo "$avg > %.2f" |bc -l) ));then\necho $window\nfi\ncd ../\ndone\ncd ../../\ndone'''%(coord_pol,dis_bound),file=checksh)
  checksh.close()
