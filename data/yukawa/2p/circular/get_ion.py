import numpy as np
import os
from os import path
import glob
import matplotlib.pyplot as plt

rootdir = os.getcwd()

Lnm = []
ion = []
pop_l1_m1 = []
pop_l1_mn1 =_l0 = []
pop_l0 = []

for folder in glob.glob(f'{rootdir}/*/'):
  if path.exists(folder+"populations.txt"):
    Lnm.append(float(folder.split("wv_")[-1][:-1]))
    with open(folder+"populations.txt","rt") as myfile:
      flag = False
      data = 0.0
      dat_l1_m1 = 0.0
      dat_l1_mn1 = 0.0
      dat_l0 = 0.0
      for myline in myfile:
        if myline[0:9] == "(2, 1, 1)":
          data += float(myline[10:].split()[0])
          dat_l1_m1 += float(myline[10:].split()[0])
          flag = True
        if myline[0:10] == "(2, 1, -1)":
          data += float(myline[11:].split()[0])
          dat_l1_mn1 += float(myline[11:].split()[0])
          flag = True
        if myline[0:9] == "(2, 1, 0)":
          data += float(myline[10:].split()[0])  
          flag = True
        if myline[0:9] == "(2, 0, 0)":
          data += float(myline[10:].split()[0])
          dat_l0 += float(myline[10:].split()[0])
          flag = True
      if flag == False:
        Lnm.pop()
      else:
        ion.append(1.0-data)
        pop_l1_m1.append(dat_l1_m1)
        pop_l1_mn1.append(dat_l1_mn1)
        pop_l0.append(dat_l0)
         
  

ion = [x for _, x in sorted(zip(Lnm,ion))]
pop_l1_m1 = [x for _, x in sorted(zip(Lnm,pop_l1_m1))]
pop_l1_mn1 = [x for _, x in sorted(zip(Lnm,pop_l1_mn1))]
pop_l0 = [x for _, x in sorted(zip(Lnm,pop_l0))]
Lnm.sort()


plt.figure()
plt.loglog(Lnm,ion,Lnm,pop_l1_m1,Lnm,pop_l1_mn1,Lnm,pop_l0)
plt.ylim([1e-7,2])
plt.grid(b=True, which='major')
plt.grid(b=True, which='minor')
plt.legend(['ion','2p-1','2p1','2s'])
plt.savefig("ion.png") 

np.savetxt("circular_Lnm_s.csv",Lnm)
np.savetxt("circular_ion_s.csv",ion)
