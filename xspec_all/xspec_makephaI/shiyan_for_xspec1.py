from astropy.io import fits
import numpy as np
import os
from zjh_data_analysis import r_baseline

import matplotlib.pyplot as plt


def write_phaI(spectrum_rate,bnname,detector,t1,t2,outfile):
	'''
	det:探头
	t1:开始时间
	t2:结束时间
	'''
	header0=fits.Header()#头文件基本信息设置
	header0.append(('creator', 'Zou', 'The name who created this PHA file'))
	header0.append(('telescop', 'Fermi', 'Name of mission/satellite'))
	header0.append(('bnname', bnname, 'Burst Name'))
	header0.append(('t1', t1, 'Start time of the PHA slice'))
	header0.append(('t2', t2, 'End time of the PHA slice'))

	hdu0=fits.PrimaryHDU(header=header0) #创建头

	a1 = np.arange(128)
	col1 = fits.Column(name='CHANNEL', format='1I', array=a1)                              #创建列
	col2 = fits.Column(name='COUNTS', format='1D', unit='COUNTS', array=spectrum_rate)     #创建列
	hdu1 = fits.BinTableHDU.from_columns([col1, col2])#创建一个bin列表
	header=hdu1.header
	header.append(('extname', 'SPECTRUM', 'Name of this binary table extension'))
	header.append(('telescop', 'GLAST', 'Name of mission/satellite'))
	header.append(('instrume', 'GBM', 'Specific instrument used for observation'))
	header.append(('filter', 'None', 'The instrument filter in use (if any)'))
	header.append(('exposure', 1., 'Integration time in seconds'))
	header.append(('areascal', 1., 'Area scaling factor'))
	header.append(('backscal', 1., 'Background scaling factor'))
	if outfile[-3:]=='pha':
		header.append(('backfile', 'A_slice_'+bnname+'_'+detector+'.bkg', 'Name of corresponding background file (if any)'))#这里有背景文件的名字
		header.append(('respfile', 'glg_cspec_'+detector+'_'+bnname+'_v02.rsp', 'Name of corresponding RMF file (if any)'))#这里有响应文件的名字
	else:
		header.append(('backfile', 'none', 'Name of corresponding background file (if any)'))
		header.append(('respfile', 'none', 'Name of corresponding RMF file (if any)'))
	header.append(('corrfile', 'none', 'Name of corresponding correction file (if any)'))
	header.append(('corrscal', 1., 'Correction scaling file'))
	header.append(('ancrfile', 'none', 'Name of corresponding ARF file (if any)'))
	header.append(('hduclass', 'OGIP', 'Format conforms to OGIP standard'))
	header.append(('hduclas1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'))
	header.append(('hduclas2', 'TOTAL', 'Indicates gross data (source + background)'))
	header.append(('hduclas3', 'COUNT', 'Indicates data stored as counts'))
	header.append(('hduvers', '1.2.1', 'Version of HDUCLAS1 format in use'))
	header.append(('poisserr', True, 'Use Poisson Errors (T) or use STAT_ERR (F)'))
	header.append(('chantype', 'PHA', 'No corrections have been applied'))
	header.append(('detchans', 128, 'Total number of channels in each rate'))
	header.append(('hduclas4', 'TYPEI', 'PHA Type I (single) or II (mulitple spectra)'))

	header.comments['TTYPE1']='Label for column 1'
	header.comments['TFORM1']='2-byte INTERGER'
	header.comments['TTYPE2']='Label for column 2'
	header.comments['TFORM2']='8-byte DOUBLE'
	header.comments['TUNIT2']='Unit for colum 2'

	hdul = fits.HDUList([hdu0, hdu1])#保存
	if(os.path.exists(outfile)):
		os.remove(outfile)#删除旧版本
	hdul.writeto(outfile)
	hdul.close()


#数据路径

topdir = '/home/laojin/trigdata/2019/'

bn = 'bn190114873'
ni = 'b1'

datalink = topdir+bn + '/'+'glg_tte_'+ni+'_' + bn + '_v00.fit'

savedir = '/home/laojin/shiyan/xspec_all/xspec_makephaI/'

if(os.path.exists(savedir) == False):
	os.makedirs(savedir)

hdu = fits.open(datalink)
trigtime = hdu['Primary'].header['TRIGTIME']
data_ = hdu['EVENTS'].data

t = data_.field(0) - trigtime
ch = data_.field(1)

ebound = hdu['EBOUNDS'].data

ch_n = ebound.field(0)
emin = ebound.field(1)
emax = ebound.field(2)

e_diff = emax-emin

'''
首先需要分能道扣除背景，然后设置能谱切片时间
在这里，设置切片时间是通过设置开始时间点和结束时间点来控制切片。

先获得各个能道的光变曲线，然后通过r_baseline获得背景曲线，最后总计切片中的内容

'''
time_start = -100  #光变曲线开始时间
time_stop = 300	   #光变曲线结束时间

slice_start = 1    #切片时间，单位s
slice_stop = 5

binsize = 1        #统计时间片大小

edges = np.arange(time_start,time_stop+binsize,binsize)    #生成统计时间片边界数组

total_rate = np.zeros(128)
bkg_rate = np.zeros(128)

total_uncertainty = np.zeros(128)
bkg_uncertainty = np.zeros(128)

for i in range(128):
	t_ch = t[ch == i]
	bin_n,bin_edges = np.histogram(t_ch,bins = edges)

	#bin_t = (bin_edges[1:]+bin_edges[:-1])*0.5          #获得单能道光变曲线
	bin_t = bin_edges[:-1]                               #只要前边界
	bin_rate = bin_n/binsize

	t_r,cs,bs = r_baseline(bin_t,bin_n)                 #获得单能道背景

	slice_index = np.where((bin_t>=slice_start) & (bin_t<=slice_stop))[0]
	slice_index = slice_index[:-1]                      #排除掉最后一个可能不完整的bin
	#print('slice index:\n',slice_index)
	#print('bs:\n',bs[slice_index])
	total_rate[i] = (bin_rate[slice_index]).mean()
	bkg_rate[i] = (bs[slice_index]).mean()

	exposure = len(slice_index)*binsize
	bkg_uncertainty[i] = np.sqrt(bkg_rate[i]/exposure)
	total_uncertainty[i] = np.sqrt(total_rate[i]/exposure)

write_phaI(total_rate,bn,ni,slice_start,slice_stop,savedir+'A_slice_'+bn+'_'+ni+'.pha')
write_phaI(bkg_rate,bn,ni,slice_start,slice_stop,savedir+'A_slice_'+bn+'_'+ni+'.bkg')

x = np.sqrt(emin*emax)

plt.figure(figsize=(10,10))
plt.subplot(1,1,1)
plt.errorbar(x,bkg_rate/e_diff,yerr = bkg_uncertainty/e_diff,color = 'blue')
plt.errorbar(x,total_rate/e_diff,yerr = total_uncertainty/e_diff,color = 'r')
plt.xlabel('energy KeV')
plt.ylabel('counts /N')
plt.xscale('log')
plt.yscale('log')
plt.savefig(savedir+'Z_slic_'+bn+'_'+ni+'.png')
plt.close()






