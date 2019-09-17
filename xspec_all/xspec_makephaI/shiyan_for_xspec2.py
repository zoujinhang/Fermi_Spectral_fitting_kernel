import numpy as np
import matplotlib.pyplot as plt
from xspec import *
import os

savedir = '/home/laojin/shiyan/xspec_all/xspec_makephaI/'
if(os.path.exists(savedir) == False):
	os.makedirs(savedir)

#brightdet = [,'A_slice_bn190114873_n7.pha','A_slice_bn190114873_n4.pha']
brightdet = ['A_slice_bn190114873_b1.pha','A_slice_bn190114873_n7.pha']

#rsplist = ['glg_cspec_n7_bn190114873_v02.rsp','glg_cspec_b1_bn190114873_v02.rsp','glg_cspec_n4_bn190114873_v02.rsp']
#baselist = ['A_slice_bn190114873_n7.bkg','A_slice_bn190114873_b1.bkg','A_slice_bn190114873_n4.bkg']

alldatastr = ' '.join(brightdet)
print(alldatastr)
print(brightdet)

AllData(alldatastr)
print('$$$$')
'''
for index,files in enumerate(brightdet):
	obj = Spectrum(files)
	obj.background = baselist[index]
	obj.response = rsplist[index]
'''
AllData.show()
AllData.ignore('1:**-200.0,40000.0-** 2-3:**-8.0,800.0-**')
print('AllData.notice')
print(AllData.notice)
print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')

#Model('cutoffpl + bbody')#模型设置
Model('grbm')

Fit.statMethod='pgstat'#拟合
Fit.nIterations=1000
Fit.query = "yes"
Fit.perform()#运行拟合


Fit.error('3.0 3')
Fit.perform()

par3=AllModels(1)(3)#第一个模型的第三个参数
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
print(par3.values)
print(par3.error)
value = par3.values[0]
value_arr1,value_arr2,ffff = par3.error
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
'''
Plot.device='/xs'
Plot.xAxis='keV'
Plot.yLog=True
'''
Plot('eeufspec')

for i in range(len(brightdet)):
	print(i)
	energies=Plot.x(i+1)
	rates=Plot.y(i+1)
	folded=Plot.model(i+1)
	xErrs=Plot.xErr(i+1)
	yErrs=Plot.yErr(i+1)
	plt.errorbar(energies,rates,xerr=xErrs,yerr=yErrs,zorder=1,ls='None')
	plt.plot(energies,folded,color='black',zorder=2)
plt.axvline(x = value,color = 'r')
plt.axvline(x = value_arr1,color = 'g')
plt.axvline(x = value_arr2,color = 'g')
plt.xlabel('Energy KeV')
plt.ylabel('rate')
plt.xscale('log')
plt.yscale('log')
plt.savefig(savedir + 'foldedspec.png')
plt.close()

Plot('eeufspec')

for i in range(len(brightdet)):
	energies=Plot.x(i+1)
	ufspec=Plot.y(i+1)
	folded=Plot.model(i+1)
	xErrs=Plot.xErr(i+1)
	yErrs=Plot.yErr(i+1)
	plt.errorbar(energies,ufspec,xerr=xErrs,yerr=yErrs,zorder=1,ls='None')
	plt.plot(energies,folded,color='black',zorder=2)
plt.axvline(x = value,color = 'r')
plt.axvline(x = value_arr1,color = 'g')
plt.axvline(x = value_arr2,color = 'g')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy KeV')
plt.ylabel('rate')
plt.savefig(savedir +'eeufspec.png')
plt.close()






