import os

from math import pi,exp,log
from scipy.constants import R

lis1=[0,0.0025,0.005,0.0075,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,0.8,1.0]
#lis1=[0,0.0002,0.0004,0.0006,0.0008,0.00125,0.0025,0.00375,0.005,0.00625,0.0075,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,0.8,1.0]

len1=len(lis1)
for j in range(1,len1-1):
  deltaG=0
  lis1.remove(lis1[-2])
  #print lis1
  for i in range(0,len(lis1)-1):
    #print lis1[i]
    sum_exp=0
    line_num=0
    force=AAA*lis1[i+1]-AAA*lis1[i]
    fil=open("%s/us.log"%(lis1[i]),"r")
    line_num=0
    if '100000' not in open('%s/us.log'%(lis1[i])).read():
      print (lis1[i],'not done')
    else:
      for line in fil:
        line_num += 1 
        step,bond,angle,dihe = line.split()
        sum_exp += exp(-(force*(float(bond)-DDDDD)**2+force*(pi/180)**2*(float(angle)-AAAAA)**2+force*(pi/180)**2*(float(dihe)-HHHHH)**2)/(R*300*0.00023900573614))
      #print sum_exp,line_num
      ave=sum_exp/line_num
      #print ave
      deltaG += -R*300*0.00023900573614*log(ave)
  print (deltaG)
