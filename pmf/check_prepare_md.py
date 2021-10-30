import os

for k in ['md_min1.out','md_min2.out','md_min3.out','md_min4.out','md_min5.out','md_nvt.out','md_npt.out','md_md.out','smd_reverse_md.out','smd_forward_md.out','md_snpt.out']:
  if os.path.isfile('%s'%(k)):
    if 'Total wall time' not in open('%s'%(k)).read():
      print (os.getcwd())
      print ('XXXXXXXXXXXXXXXXXXXXXX %s not done'%(k))
  else:
    print (os.getcwd())
    print ('YYYYYYYYYYYYYYYYYYY %s file not exist'%(k))
