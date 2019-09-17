import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from glob import glob

datadir = '/home/laojin/shiyan/data/'

filename = 'glg_cspec_n4_bn190114873_v02.rsp2'
savedir = '/home/laojin/shiyan/xspec_all/'

rsp_n = 7 #响应矩阵的第几个。
#ttefile = glob(datadir+'glg_cspec_n4_bn190114873_v*.pha')
#print(ttefile)

rsp2 = fits.open(datadir + filename)
num_of_rsp2 = len(rsp2)

data = rsp2[rsp_n].data

ebounds = rsp2[1].data #在列表1中，内容是各个探测器的能道的能量范围。
e_min = ebounds['E_MIN']
e_max = ebounds['E_MAX']

detector_e = np.sqrt(e_min*e_max)#为方便取对数，用这种开方的形式计算能道中心的数值



rsp_martrix = data['MATRIX']  #响应矩阵能量量化的边界
rsp_martrix[-1] = np.zeros(128)
rsp_martrix = np.vstack(rsp_martrix)
#rsp_martrix = np.array(rsp_martrix)
print(rsp_martrix[0])
print('matrix:\n',rsp_martrix)
energ_lo = data['ENERG_LO']
energ_hi = data['ENERG_HI']

rsp_e = np.sqrt(energ_lo*energ_hi)

X,Y = np.meshgrid(rsp_e,detector_e)
plt.contourf(X,Y,rsp_martrix.T)
plt.contour(X,Y,rsp_martrix.T,100)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('rsp e')
plt.ylabel('detector_e')
plt.savefig(savedir+'rsp_matrix_n_'+str(rsp_n)+'.png')
plt.close()

print(len(data['MATRIX'][0]))

#print('response matrix:\n',data['MATRIX'][0])


all_sum = 0
for i in range(len(data['MATRIX'])):
	SUM = data['MATRIX'][i].sum()
	all_sum = all_sum + SUM
print('all_sum :',all_sum)








