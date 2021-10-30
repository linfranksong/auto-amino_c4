#!/mnt/home/songlin3/anaconda3/bin/python
from __future__ import print_function
import os
from optparse import OptionParser
import numpy as np
import math
from freemd import freemd
from pmf import pmf
from firstpass import firstpass
from prepdisandihe import prepdisandihe

#----------------------------------------------------------------------------#
#                            Main Program                                    #
#----------------------------------------------------------------------------#
parser = OptionParser("Usage: auto-amino_c4.py -i inputfile -s/--step step_number \n")

parser.add_option("-i", dest="inputf", type='string',
                  help="Input file name")
parser.add_option("-s", "--step", dest="step", type='string',
                  help="Step number")
(options, args) = parser.parse_args()


# Print the title of the program
version = '1.0'

#---------------------------Default values------------------------------------
disrst_strength=["600"]

# About the program 
#if options.step not in ['prep1','prep2', 'run', 'freemd1','freemd2','qmmmmd1','qmmmmd2','prep1disandihe','prep2disandihe']:
#    raise ('Invalid step number chosen. please choose among the following values: prep1,prep2, run, freemd1, freemd2,qmmmmd1,qmmmmd2,prep1disandihe,prep1disandihe')
#----------------------------Read the input files------------------------------
rinput = open(options.inputf, 'r')
for line in rinput:
    line = line.split()
    if line[0] == "hfe":
        hfe = line[1]
    elif line[0] == "gaff":
        gaff = line[1]
    elif line[0] == "protein_ff":
        protein_ff = line[1]
    elif line[0] == "water_model":
        water_model = line[1]
        if water_model=="SPCE":
           water_box="SPC"
        else:
           water_box=water_model
    elif line[0] == "amber_frcmod_files":
        amber_frcmod_num = (len(line)-1)
        amber_frcmod_files = line[1:]
    elif line[0] == "additional_frcmod_files":
        additional_frcmod_num = (len(line)-1)
        additional_frcmod_files = line[1:]
    elif line[0] == "amber_mol2_files":
        amber_mol2_num = (len(line)-1)
        amber_mol2_files = line[1:]
    elif line[0] == "additional_mol2_files":
        additional_mol2_num = (len(line)-1)
        additional_mol2_files = line[1:]
    elif line[0] == "Na+_ion":
        Na_ion = line[1]
    elif line[0] == "CYX-CYX_bond":
        cyx_num = int((len(line)-1)/2)
        cyx_bond = line[1:]
    elif line[0] == "Cl-_ion":
        Cl_ion = line[1]
    elif line[0] == "box_size":
        box_size = line[1]
    elif line[0] == "metal":
        metal = line[1]
    elif line[0] == "HZ_resid":
        HZ_resid = int(line[1])
    elif line[0] == "try_pol":
        try_pol=round(float(line[1]),3)
    elif line[0] == "coord_Rmin_Epsilon_pol":
        coord_rmin=line[1]
        coord_epsilon=line[2]
        coord_pol=round(float(line[3]),3)
    elif line[0] == "restraint_atomid":
        coord_atomid1=line[1]
        coord_atomid2=line[2]
        coord_atomid3=line[3]
    elif line[0] == "dis_ang_dihe_value":
        dis=line[1]
        ang=line[2]
        dihe=line[3]
    elif line[0] == "dihe_mini":
        m_dihe = line[1]
    elif line[0] == "dihe_maxi":
        p_dihe = line[1]
    elif line[0] == "restraint_strength":
        disrst_num = (len(line)-1)
        disrst_strength = line[1:]
    elif line[0] == "min_steps":
        min_steps = int(line[1])
    elif line[0] == "heat_steps":
        heat_steps = int(line[1])
    elif line[0] == "equi_steps":
        equi_steps = int(line[1])
    elif line[0] == "sampl_steps":
        sampl_steps = int(line[1])
    elif line[0] == "slm_header":
        slm_header = line[1]
    elif line[0] == "time_step":
        time_step = line[1]
    elif line[0] == "cutoff":
        cutoff = line[1]
rinput.close()

#--------------------------------Print the variables---------------------------
print("The input file you are using is : %s" %options.inputf)
print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

# Necessary inputs
try:
    print('The force fileds are : ', gaff, protein_ff, water_model)
except:
    raise auto_pka_error('force fields need to be provided.')


######## step1 ########
print ()
if (options.step == 'freemd'):
  print("Building topology....")
###################### modify input-noH.pdb, temp.dat, write addC4.in, c4.txt#############
  print('Load frcmod and mol2 files : ', amber_frcmod_files,additional_frcmod_files, amber_mol2_files, additional_mol2_files)
  print('Add ions: %s Na+ and %s Cl-'%(Na_ion,Cl_ion))
  print('Water box size is solvateoct %s A'%(box_size))
  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

  freemd(metal,gaff,protein_ff,water_model,water_box,amber_frcmod_files,additional_frcmod_files,additional_frcmod_num,amber_frcmod_num, amber_mol2_files,amber_mol2_num,additional_mol2_files,additional_mol2_num,Na_ion,Cl_ion,cyx_num,cyx_bond,box_size,HZ_resid,coord_pol,coord_rmin,coord_epsilon,coord_atomid1,coord_atomid2,coord_atomid3,dis,ang,dihe,m_dihe,p_dihe,disrst_strength,min_steps,heat_steps,equi_steps,sampl_steps,slm_header,cutoff,time_step)

if (options.step == 'prepdisandihe'):
  print("Building topology....")
###################### modify input-noH.pdb, temp.dat, write addC4.in, c4.txt#############
  print('Load frcmod and mol2 files : ', amber_frcmod_files,additional_frcmod_files, amber_mol2_files, additional_mol2_files)
  print('Add ions: %s Na+ and %s Cl-'%(Na_ion,Cl_ion))
  print('Water box size is solvateoct %s A'%(box_size))
  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

  prepdisandihe(metal,hfe,gaff,protein_ff,water_model,water_box,amber_frcmod_files,additional_frcmod_files,additional_frcmod_num,amber_frcmod_num, amber_mol2_files,amber_mol2_num,additional_mol2_files,additional_mol2_num,Na_ion,Cl_ion,cyx_num,cyx_bond,box_size,HZ_resid,coord_pol,coord_rmin,coord_epsilon,coord_atomid1,coord_atomid2,coord_atomid3,dis,ang,dihe,m_dihe,p_dihe,disrst_strength,min_steps,heat_steps,equi_steps,sampl_steps,slm_header,cutoff,time_step)

if (options.step == 'pmf'):
  pmf(try_pol,metal,gaff,protein_ff,water_model,water_box,amber_frcmod_files,additional_frcmod_files,additional_frcmod_num,amber_frcmod_num, amber_mol2_files,amber_mol2_num,additional_mol2_files,additional_mol2_num,Na_ion,Cl_ion,cyx_num,cyx_bond,box_size,HZ_resid,coord_pol,coord_rmin,coord_epsilon,coord_atomid1,coord_atomid2,coord_atomid3,dis,ang,dihe,m_dihe,p_dihe,disrst_strength,min_steps,heat_steps,equi_steps,sampl_steps,slm_header,cutoff,time_step)

if (options.step == 'firstpass'):
  firstpass(metal,gaff,protein_ff,water_model,water_box,amber_frcmod_files,additional_frcmod_files,additional_frcmod_num,amber_frcmod_num, amber_mol2_files,amber_mol2_num,additional_mol2_files,additional_mol2_num,Na_ion,Cl_ion,cyx_num,cyx_bond,box_size,HZ_resid,coord_pol,coord_rmin,coord_epsilon,coord_atomid1,coord_atomid2,coord_atomid3,dis,ang,dihe,m_dihe,p_dihe,disrst_strength,min_steps,heat_steps,equi_steps,sampl_steps,slm_header,cutoff,time_step)
quit()
