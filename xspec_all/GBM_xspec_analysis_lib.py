# General analyses of GBM catalog bursts 
# last modified: Apr. 29, 2019

from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
from glob import glob
import pandas as pd
import numpy as np
import h5py
from scipy import stats
import os
import sys
from multiprocessing import Pool
import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)
import rpy2.robjects as robjects
from rpy2.robjects import r
import rpy2.robjects.numpy2ri
robjects.numpy2ri.activate()
robjects.r("library(baseline)")
from xspec import *


databasedir='/home/yao/bn'
#databasedir='/home/yujie/downburstdata/data'

NaI=['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']
BGO=['b0','b1']
Det=['b0','b1','n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb']

ch1=3
ch2=124
ncore=10


def norm_pvalue(sigma=2.0):
	p = stats.norm.cdf(sigma)-stats.norm.cdf(-sigma)
	return p


def write_phaI(spectrum_rate,bnname,det,t1,t2,outfile):
	'''
	det:探头
	t1:开始时间
	t2:结束时间
	'''
	header0=fits.Header()#头文件基本信息设置
	header0.append(('creator', 'Shao', 'The name who created this PHA file'))
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
		header.append(('backfile', det+'.bkg', 'Name of corresponding background file (if any)'))#这些关键字会引导程序加载背景文件
		header.append(('respfile', det+'.rsp', 'Name of corresponding RMF file (if any)'))#这些关键字会引导程序加载响应文件
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
	hdul.writeto(outfile)


def copy_rspI(bnname,det,outfile):
	shortyear=bnname[2:4]
	fullyear='20'+shortyear
	datadir=databasedir+'/'+fullyear+'/'+bnname+'/'
	rspfile=glob(datadir+'/'+'glg_cspec_'+det+'_'+bnname+'_v*.rsp')
	assert len(rspfile)==1, 'response file is missing for '+'glg_cspec_'+det+'_'+bnname+'_v*.rsp'
	rspfile=rspfile[0]
	os.system('cp '+rspfile+' '+outfile)
	

class GRB:
	def __init__(self,bnname):
		self.bnname=bnname
		resultdir=os.getcwd()+'/results/'
		self.resultdir=resultdir+'/'+bnname+'/'

		shortyear=self.bnname[2:4]
		fullyear='20'+shortyear
		self.datadir=databasedir+'/'+fullyear+'/'+self.bnname+'/'
		self.dataready=True
		for i in range(14):
			ttefile=glob(self.datadir+'glg_tte_'+Det[i]+'_'+self.bnname+'_v*.fit')
			if not len(ttefile)==1:
				self.dataready=False
			else:
				hdu=fits.open(ttefile[0])
				event=hdu['EVENTS'].data.field(0)
				if len(event)<10:
					self.dataready=False
		if self.dataready:
			if not os.path.exists(resultdir):
				os.makedirs(resultdir)
			if not os.path.exists(self.resultdir):
				os.makedirs(self.resultdir)
			self.baseresultdir=self.resultdir+'/base/'
			self.phaIresultdir=self.resultdir+'/phaI/'

			# determine GTI1 and GTI2
			GTI_t1=np.zeros(14)
			GTI_t2=np.zeros(14)
			for i in range(14):
				ttefile=glob(self.datadir+'glg_tte_'+Det[i]+'_'+self.bnname+'_v*.fit')
				hdu=fits.open(ttefile[0])
				trigtime=hdu['Primary'].header['TRIGTIME']
				data=hdu['EVENTS'].data
				time=data.field(0)-trigtime
				GTI0_t1=time[0]
				GTI0_t2=time[-1]
				timeseq1=time[:-1]
				timeseq2=time[1:]
				deltime=timeseq2-timeseq1
				delindex=deltime>5 
				if len(timeseq1[delindex])>=1:
					GTItmp_t1=np.array(np.append([GTI0_t1],timeseq2[delindex]))
					GTItmp_t2=np.array(np.append(timeseq1[delindex],[GTI0_t2]))
					for kk in np.arange(len(GTItmp_t1)):
						if GTItmp_t1[kk]<=0.0 and GTItmp_t2[kk]>=0.0:
							GTI_t1[i]=GTItmp_t1[kk]
							GTI_t2[i]=GTItmp_t2[kk]
				else:
					GTI_t1[i]=GTI0_t1
					GTI_t2[i]=GTI0_t2
			self.GTI1=np.max(GTI_t1)
			self.GTI2=np.min(GTI_t2)



	def base(self,baset1=-50,baset2=300,binwidth=0.1):
		self.baset1=np.max([self.GTI1,baset1])
		self.baset2=np.min([self.GTI2,baset2])
		self.binwidth=binwidth
		self.tbins=np.arange(self.baset1,self.baset2+self.binwidth,self.binwidth)
		assert self.baset1<self.baset2, self.bnname+': Inappropriate base times!'
		if not os.path.exists(self.baseresultdir):
			os.makedirs(self.baseresultdir)
			expected_pvalue = norm_pvalue()
			f=h5py.File(self.baseresultdir+'/base.h5',mode='w')
			for i in range(14):#14个探头
				grp=f.create_group(Det[i])#创建一个探头的组
				ttefile=glob(self.datadir+'/'+'glg_tte_'+Det[i]+'_'+\
                     							self.bnname+'_v*.fit')#找到相应探头的文件名
				hdu=fits.open(ttefile[0])	#打开相应探头文件
				trigtime=hdu['Primary'].header['TRIGTIME']
				data=hdu['EVENTS'].data
				timedata=data.field(0)-trigtime
				chdata=data.field(1)
				for ch in range(128):     #这里应该就是大名鼎鼎的分能道扣除背景
					time_selected=timedata[chdata==ch]
					histvalue, histbin=np.histogram(time_selected,bins=self.tbins)
					rate=histvalue/binwidth
					r.assign('rrate',rate) 
					r("y=matrix(rrate,nrow=1)")
					fillPeak_hwi=str(int(5/binwidth))
					fillPeak_int=str(int(len(rate)/10))
					r("rbase=baseline(y,lam = 6, hwi="+fillPeak_hwi+", it=10,\
								 int ="+fillPeak_int+", method='fillPeaks')")
					r("bs=getBaseline(rbase)")
					r("cs=getCorrected(rbase)")
					bs=r('bs')[0]
					cs=r('cs')[0]
					corrections_index= (bs<0)
					bs[corrections_index]=0#baseline小于0的部分强制为0
					cs[corrections_index]=rate[corrections_index]#扣除背景的部分等于原来的部分。
					f['/'+Det[i]+'/ch'+str(ch)]=np.array([rate,bs,cs])#将每一个探头中每一个能道中的背景分别保存。
			f.flush()
			f.close()
	



													
	def phaI(self,slicet1=0,slicet2=5):
		if not os.path.exists(self.phaIresultdir):
			os.makedirs(self.phaIresultdir)
		nslice=len(os.listdir(self.phaIresultdir))
		sliceresultdir=self.phaIresultdir+'/slice'+str(nslice)+'/'
		os.makedirs(sliceresultdir)
		fig, axes= plt.subplots(7,2,figsize=(32, 30),sharex=False,sharey=False)#这里创建了一个画布，用于之后的图像输出
		sliceindex= (self.tbins >=slicet1) & (self.tbins <=slicet2)#这里进行了时间切片
		valid_bins=np.sum(sliceindex)-1
		assert valid_bins>=1, self.bnname+': Inappropriate phaI slice time!'
		'''
		assert 作用是测试一个条件(condition)是否成立，如果不成立，则抛出异常。

		'''
		f=h5py.File(self.baseresultdir+'/base.h5',mode='r')
		for i in range(14):#14个探头分别做了文件
			total_rate=np.zeros(128)
			bkg_rate=np.zeros(128)
			total_uncertainty=np.zeros(128)
			bkg_uncertainty=np.zeros(128)
			ttefile=glob(self.datadir+'/glg_tte_'+Det[i]+'_'+self.bnname+'_v*.fit')
			hdu=fits.open(ttefile[0])
			ebound=hdu['EBOUNDS'].data
			emin=ebound.field(1)
			emax=ebound.field(2)
			energy_diff=emax-emin
			energy_bins=np.concatenate((emin,[emax[-1]]))
			for ch in np.arange(128):
				base=f['/'+Det[i]+'/ch'+str(ch)][()][1]
				rate=f['/'+Det[i]+'/ch'+str(ch)][()][0]
				bkg=base[sliceindex[:-1]][:-1]
				total=rate[sliceindex[:-1]][:-1]
				bkg_rate[ch]=bkg.mean()
				total_rate[ch]=total.mean()
				#-----------------------------------------------------------------------
				exposure=len(bkg)*self.binwidth#切片长度乘bin宽度
				bkg_uncertainty[ch]=np.sqrt(bkg_rate[ch]/exposure)#这个当做背景误差
				total_uncertainty[ch]=np.sqrt(total_rate[ch]/exposure)#总的误差
				#--------------------------------------------------------------------------
			#plot both rate and bkg as count/s/keV
			write_phaI(bkg_rate,self.bnname,Det[i],slicet1,slicet2,sliceresultdir+'/'+Det[i]+'.bkg')#对每个探头进行制作
			write_phaI(total_rate,self.bnname,Det[i],slicet1,slicet2,sliceresultdir+'/'+Det[i]+'.pha')#有个小小的问题，这里面的能量是怎么换算过去的？
			copy_rspI(self.bnname,Det[i],sliceresultdir+'/'+Det[i]+'.rsp')
			bkg_diff=bkg_rate/energy_diff
			total_diff=total_rate/energy_diff
			x=np.sqrt(emax*emin)
			axes[i//2,i%2].errorbar(x,bkg_diff,yerr=bkg_uncertainty/energy_diff,linestyle='None',color='blue')
			axes[i//2,i%2].errorbar(x,total_diff,yerr=total_uncertainty/energy_diff,linestyle='None',color='red')
			bkg_diff=np.concatenate(([bkg_diff[0]],bkg_diff))
			total_diff=np.concatenate(([total_diff[0]],total_diff))
			axes[i//2,i%2].plot(energy_bins,bkg_diff,linestyle='steps',color='blue')
			axes[i//2,i%2].plot(energy_bins,total_diff,linestyle='steps',color='red')
			axes[i//2,i%2].set_xscale('log')
			axes[i//2,i%2].set_yscale('log')
			axes[i//2,i%2].tick_params(labelsize=25)
			axes[i//2,i%2].text(0.85,0.85,Det[i],transform=\
										axes[i//2,i%2].transAxes,fontsize=25)
		fig.text(0.07, 0.5, 'Rate (count s$^{-1}$ keV$^{-1}$)', ha='center',\
							va='center', rotation='vertical',fontsize=30)
		fig.text(0.5, 0.05, 'Energy (keV)', ha='center', va='center',\
															fontsize=30)	
		plt.savefig(sliceresultdir+'/PHA_rate_bkg.png')
		plt.close()
		f.close()


	def specanalyze(self,slicename):
		slicedir=self.phaIresultdir+'/'+slicename+'/'
		os.chdir(slicedir)#到数据所在的目录下
		# select the most bright two NaIs (in channels 6-118) 
		# and more bright one BGO (in channels 4-124):
		BGOtotal=np.zeros(2)
		NaItotal=np.zeros(12)
		for i in range(2):
			phahdu=fits.open(slicedir+'/'+BGO[i]+'.pha')#打开光谱文件
			bkghdu=fits.open(slicedir+'/'+BGO[i]+'.bkg')#打开背景文件
			pha=phahdu['SPECTRUM'].data.field(1)#光谱文件
			bkg=bkghdu['SPECTRUM'].data.field(1)#背景光谱文件
			src=pha-bkg
			plt.plot(src[4:125])#
			plt.savefig(BGO[i]+'.png')
			plt.close()
			BGOtotal[i]=src[4:125].sum()
		for i in range(12):
			phahdu=fits.open(slicedir+'/'+NaI[i]+'.pha')
			bkghdu=fits.open(slicedir+'/'+NaI[i]+'.bkg')
			pha=phahdu['SPECTRUM'].data.field(1)
			bkg=bkghdu['SPECTRUM'].data.field(1)
			src=pha-bkg
			plt.plot(src[6:118])
			plt.savefig(NaI[i]+'.png')
			plt.close()
			NaItotal[i]=src[6:118].sum()
		BGOindex=np.argsort(BGOtotal)#各个探头
		NaIindex=np.argsort(NaItotal)

		brightdet=[BGO[BGOindex[-1]],NaI[NaIindex[-1]],NaI[NaIindex[-2]]]#亮探头，1个bgo 两个nai
		
		# use xspec

		#alldatastr='b0.pha n4.pha n3.pha'	
		alldatastr=' '.join([det+'.pha' for det in brightdet])#整理光谱文件。
		print(alldatastr)
		#input('--wait--')
		AllData(alldatastr)#添加光谱文件
		AllData.show()
		AllData.ignore('1:**-200.0,40000.0-** 2-3:**-8.0,800.0-**')#这里前一个是bgo，后面是nai，选择的是能量范围
		print(AllData.notice)
		Model('grbm')#模型设置
		Fit.statMethod='pgstat'#拟合
		Fit.nIterations=1000
		Fit.query = "yes"
		Fit.perform()#运行拟合

		Fit.error('3.0 3')
		Fit.perform()

		par3=AllModels(1)(3)#模型参数
		print(par3.error)

		#-------
		#模型打印
		Plot.device='/xs'
		Plot.xAxis='keV'
		Plot.yLog=True
		Plot('eeufspec')

		for i in (1,2,3):
			energies=Plot.x(i)
			rates=Plot.y(i)
			folded=Plot.model(i)
			xErrs=Plot.xErr(i)
			yErrs=Plot.yErr(i)
			plt.errorbar(energies,rates,xerr=xErrs,yerr=yErrs,zorder=1,ls='None')
			plt.plot(energies,folded,color='black',zorder=2)
		plt.xscale('log')
		plt.yscale('log')
		plt.savefig('foldedspec.png')
		plt.close()
		Plot('eeufspec')

		for i in (1,2,3):
			energies=Plot.x(i)
			ufspec=Plot.y(i)
			folded=Plot.model(i)
			xErrs=Plot.xErr(i)
			yErrs=Plot.yErr(i)
			plt.errorbar(energies,ufspec,xerr=xErrs,yerr=yErrs,zorder=1,ls='None')
			plt.plot(energies,folded,color='black',zorder=2)
		plt.xscale('log')
		plt.yscale('log')
		plt.savefig('eeufspec.png')
		plt.close()

	def removebase(self):
		os.system('rm -rf '+self.baseresultdir)


grb=GRB('bn110920546')
grb.base(baset1=-50,baset2=200,binwidth=0.064)
grb.phaI(slicet1=51,slicet2=55)#这里是生成光谱文件的。两个参数为切片的起始和结束时间。
grb.specanalyze('slice'+str(0))
#grb.removebase()
	
