import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
from matplotlib.backends.backend_pdf import PdfPages
# some style:
import seaborn as sns
plt.rcParams['lines.linewidth'] = 1
os.chdir('/media/glaser/data/Nijmegen/Research/WillCausetCode/evolution_lisa/2017_code_evolution_lisa/')

files='dat/CSwSpinstate.dat'
actionf='dat/act/action.evo.act.dat'
intervalsf='dat/int/intervals.evo.act.dat'
data=np.loadtxt(files)
action=np.loadtxt(actionf)
intervals=np.loadtxt(intervalsf)

s=[]
u=[]
v=[]

for i in range(len(data)):
    if i%3==0:
        s.append(data[i])
    elif i%3==1:
        u.append(data[i])
    else:
        v.append(data[i])

for i in range(len(s)):
    plt.axes().set_aspect('equal')
    plt.scatter(u[i],v[i],c=s[i],s=100)
    plt.show()


plt.plot(action[100:,0])
plt.plot(action[100:,1])
