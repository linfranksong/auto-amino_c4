import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score


os.system("rm *csv *png")

donelist=[]
done=open("done.txt","r")
for line in done:
  line=line.split('\t')
  metal,hfe,p,watermodel=line[0].split('_')
  donelist.append("%s_%s"%(metal,watermodel))
done.close()



dic1={}
result=open("result.txt","r")
for line in result:
  line=line.split('\t')
  metal,hfe,p,watermodel=line[0].split('_')
  affinity=line[1]
  fil=open("%s_%s.csv"%(metal,watermodel),"a")
  fil.seek(0,2)
  fil.write("%s	%s"%(p,affinity))
  fil.close()
  dic1[metal]=hfe
result.close()

dic={}
expr=open("expr.txt","r")
for line in expr:
  line=line.split('\t')
  metal=line[0]
  expr_affinity=line[1].strip('\n')
  dic[metal]=expr_affinity
expr.close()

namlist=[]
for name in os.listdir("."):
  namlist.append(name)
namlist.sort()

newpu=open("new-input.txt","w")
for name in namlist:
  if ".csv" in name:
    data_pol=[]
    data_affinity=[]
    name=name.strip('.csv')
    metal,watermodel=name.split('_')
    if "%s_%s"%(metal,watermodel) in donelist:
      print ("%s_%s done"%(metal,watermodel))
    else:
      data=pd.read_csv('%s.csv'%(name), sep='\t',header=None,engine='python')
      data.columns = ["pol", "affinity"]

      data_pol=data.iloc[:,0].values.reshape(-1,1)
      data_affinity=data.iloc[:,1].values.reshape(-1,1)

      #print(metal,watermodel)
      if len(data.index) >= 2:
        regr = linear_model.LinearRegression()
        #print (data_affinity)
        #print (data_pol)
        regr.fit(data_affinity, data_pol)

        new_aray=[]
        aray=[]
        n_aray=[]
        f_aray=[]

        expr_affinity=dic[metal]
        new_aray.append(expr_affinity)
        n_aray=np.array(new_aray)
        aray=n_aray.astype(np.float64)

        #aray.reshape(-1, 1)
        #aray.reshape(1,-1)
        f_aray.append(aray)
        #print (f_aray)

        pol_predict = regr.predict(f_aray)

        #print('Coefficients: \n', regr.coef_)

        newpu.write("%s_%s_%s_%s\n"%(metal,round(float(dic1[metal]),1),round(float(pol_predict[0][0]),2),watermodel))
        #print (pol_predict[0][0])

        #plt.scatter(data_affinity,data_pol,color='black')
        #plt.plot(data_affinity, data_pol, color='red')
        #plt.savefig('%s_%s.png'%(metal,watermodel))
      else:
        print ("%s_%s does not have two pol-affinity points, please add"%(metal,watermodel))
newpu.close()
