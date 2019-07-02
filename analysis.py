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


def atoi(text):
    try:
        return float(text)
    except ValueError:
        return text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    '''
    return [ atoi(c) for c in re.split('_', os.path.splitext(text)[0]) ]

def autocorr(x):
    y=x-x.mean()
    result = np.correlate(y, y, mode = 'full')
    maxcorr = np.argmax(result)
    #print 'maximum = ', result[maxcorr]
    result = result / result[maxcorr]     # <=== normalization

    return result[result.size/2:]

def bootstrap(data,auco,nsamp):
    means=[]
    varis=[]
    for x in range(nsamp):
        sample=data[np.random.randint(len(data),size=len(data)/auco)]
        sample=np.array(sample)
        means.append([s.mean()  for s in sample.transpose()])
        varis.append([s.var()  for s in sample.transpose()])
    #print(means)
    means=np.array(means)
    varis=np.array(varis)
    mean=[m.mean() for m in means.transpose()]

    var=[m.mean() for m in varis.transpose()]
    errmean=[np.sqrt(m.var()) for m in means.transpose()]
    errvar=[np.sqrt(m.var()) for m in varis.transpose()]
    return [mean, errmean, var, errvar]
## definitions for this geometry


#palcol=sns.diverging_palette(255, 0, l=60, n=25, center="dark")
palcol=sns.color_palette("cubehelix", 4)


os.chdir('/media/glaser/data/Nijmegen/Research/WillCausetCode/evolution_lisa/2017_code_evolution_lisa/dat/act')
os.path.abspath(os.path.curdir)

datafiles=glob.glob('*.dat')
datafiles
#datafiles=glob.glob('action_j_N_'+str(N)+'*_j_-0.1.dat')
#olddata=np.loadtxt(str(N)+'-'+str(eps)+'-MEAN')
datafiles.sort(key=natural_keys)
data=[np.loadtxt(f) for f in datafiles]

data

palcol=sns.color_palette("cubehelix",len(jval)+1)


data=np.array(data)
cutoff=5000
## let's look at autocorrelations
aucodat=[autocorr(f[cutoff:10*cutoff,0]) for f in data]

for p,b,j in zip(aucodat,beta,jising):
    plt.plot(p[0:500])
    plt.plot([1/np.exp(1) for x in range(500)])
    plt.title("b="+str(b)+" j="+j)
    plt.ylim([0,1])
    plt.show()

# this justifies taking every 1000th measurement and calling the data uncorrelated
jump=500
#%%

mean=[]
errm=[]
var=[]
errv=[]
absM=[]

for d in data:
    mm, em, v, ev = bootstrap(d[cutoff:],jump,1000)
    mean.append(mm)
    var.append(v)
    errm.append(em)
    errv.append(ev)
    mm,em,v,ev=bootstrap(np.array([np.abs(d[cutoff:,1])]).transpose(),jump,1000)
    absM.append([*mm,*em,*v,*ev])



mean=np.array(mean)
var=np.array(var)
errm=np.array(errm)
errv=np.array(errv)
absM=np.array(absM)


jsortedData=[]


for jj in bval:
    temp=[]
    for m,v,j,b,em,ev,aM in zip(mean,var,jising,beta,errm,errv,absM):
        if(b==jj and float(j)!=0):
            temp.append([float(b),float(j),*m,aM[0],*v,aM[2],*em,aM[1],*ev,aM[3]])
    jsortedData.append(np.array(temp))

#%%

for d,bf,j in zip(data,beta,jising):
    plt.plot(d[5000:]/float(j),label=" Ising/j")
    plt.title("beta="+bf+"j="+j)
    plt.legend()
    plt.show()

#%%
off=2*(len(mean[0])+1)

#%%
plt.title("Action")
for j in jsortedData:
    plt.errorbar(j[:,1],j[:,2],yerr=j[:,2+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.xlabel("$j$")
plt.ylabel("$S$")
#%%
plt.title("Magnetisation")
for j in jsortedData:
    plt.errorbar(j[:,1],np.abs(j[:,3]),yerr=j[:,3+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.xlabel("$j$")
plt.ylabel("$M$")
#%%
plt.title("Spin correlation")
for j in jsortedData:
    plt.errorbar(j[:,1],j[:,4],yerr=j[:,4+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.xlabel("$j$")
plt.ylabel("$S_c$")
#%%
plt.title("absolute magnetisation")
for j in jsortedData:
    plt.errorbar(j[:,1],j[:,5],yerr=j[:,5+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.xlabel("$j$")
plt.ylabel("$|M|$")

#%%
major_ticks = np.arange(-10, 10, 1)
minor_ticks = np.arange(-10, 10, 0.1)
#plt.title("Action by j")
for j in jsortedData:
    plt.errorbar(j[:,1],j[:,2]/j[:,1],yerr=j[:,2+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.axes().set_xticks(major_ticks)
plt.axes().set_xticks(minor_ticks, minor=True)
plt.xlabel("$j$")
plt.ylabel("$S$")
plt.savefig("../../../../IsingCauset/writing/figures/Action_by_j_random_50.pdf")
#%%
#plt.title("Action Var")
for j in jsortedData:
    plt.errorbar(j[:,1],j[:,6],yerr=j[:,6+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.axes().set_xticks(major_ticks)
plt.axes().set_xticks(minor_ticks, minor=True)
plt.xlabel("$\\beta$")
plt.ylabel("$Var$")
plt.savefig("../../../../IsingCauset/writing/figures/VarAction_random_50.pdf")
#%%
plt.title("Action Var")
for j in jsortedData:
    plt.errorbar(j[:,1],j[:,6],yerr=j[:,6+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.axes().set_xticks(major_ticks)
plt.axes().set_xticks(minor_ticks)
plt.xlabel("$\\beta$")
plt.ylabel("$Var$")
plt.xlim([3,5])
#%%
plt.title("Magnetisation Var")
for j in jsortedData:
    plt.errorbar(j[:,1],j[:,7],yerr=j[:,7+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.axes().set_xticks(major_ticks)
plt.axes().set_xticks(minor_ticks, minor=True)
plt.xlabel("$\\beta$")
plt.ylabel("$Var_M$")
#%%
plt.title("Spin Correlation Var")
for j in jsortedData:
    plt.errorbar(j[:,1],j[:,8],yerr=j[:,8+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.axes().set_xticks(major_ticks)
plt.axes().set_xticks(minor_ticks, minor=True)
plt.xlabel("$\\beta$")
plt.ylabel("$Var_s$")
plt.savefig("../../../../IsingCauset/writing/figures/VarCorrelation_random_50.pdf")
#%%
#plt.title("AbsMagnetisation var")
for j in jsortedData:
    plt.errorbar(j[:,1],j[:,9],yerr=j[:,9+off],c=palcol[jval.index(str(j[0,1]))])
plt.legend(loc=3)
plt.axes().set_xticks(major_ticks)
plt.axes().set_xticks(minor_ticks, minor=True)
plt.xlabel("$\\beta$")
plt.ylabel("$Var_{|M|}$")
plt.savefig("../../../../IsingCauset/writing/figures/VarAbsMagnetisation_random_50.pdf")
