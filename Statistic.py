#-*-coding: UTF-8-*-
import numpy as np
import matplotlib.pyplot as plt
import random
import os
import itertools


class Regression():
	def __init__(self,data):
		self.data=data

	def __corrCoeff (self,yestimate,yrecord):
		#correlation coeffient computer
		#input:yestimate--回归得到的预测值（Nx1）列向量
			   #yrecord--实验得到的值（Nx1）列向量
		#output:correlation coeffient
		corr=np.corrcoef(yestimate.T,yrecord.T)
		coeff=corr[0,1]
		return coeff
	def __std (self,yestimate,yrecord):
		SSres=sum(np.array(np.array(yrecord)-np.array(yestimate))**2)
		std=SSres/float(len(yrecord))
		return std**0.5
		
#######################################
	def __deterCoeff (self,yestimate,Y):
		#compute coeffient of determination R2
		#input:yestimate--回归得到的预测值（Nx1）列向量
				##Y--实验得到的值（Nx1）列向量
		##output:R2
		ymean=float(sum(Y))/len(Y)
		SStot=sum(np.array(np.array(Y)-ymean)**2)
		SSres=sum(np.array(np.array(Y)-np.array(yestimate))**2)
		R2=1-np.array(SSres)/SStot
		return R2
#######################################		
	def linearRegression (self):
		#一般线性回归，data为数据，其中最后一列为y值，x的排序按xo,x1,x2....
		#y=xT*w,其中x为变量列矩阵，w为回归系数列矩阵
		#error=sum((yi-xi.trans*w)**2)
		#w=inv(xTx)*xT*y
		#return:w-regression coeffient coeff-correlation coeffient R2-coeffient of determination
		matdata=np.mat(self.data)
		X=matdata[:,:-1]
		Y=matdata[:,-1]
		xTx=X.T*X
		if np.linalg.det(xTx)==0.0:
			print ("matrix singularity")
			return
		else:
			w=xTx.I*(X.T*Y)
		yestimate=X*w
		#correlation coeffient compute
		coeff=self.__corrCoeff(yestimate,Y)
		std=self.__std(yestimate,Y)

		#compute coeffient of determination R2
		R2=self.__deterCoeff(yestimate,Y)

		return (w,coeff,R2,yestimate,std)
#######################################
	def __localRegression (self,k,xarrloc):
		#Locally weighted linear regression, k-weight parameter xarr-（x0,x1,...）单个点集
		#what=(X.T*W*X).inv*X.T*W*y
		#W is a matrix 
		#采用核来对附近的点赋予更高的权重，这里采用高斯核 w(i,i)=exp(abs(xi-x)**2/(-2*k**2))
		#点x与xi越近，w(i,i)越大
		#return local regression coeffient
		matdata=np.mat(self.data)
		X=matdata[:,0:2]
		Y=matdata[:,2] 
		m=np.shape(X)[0]
		weights=np.mat(np.eye(m))
		for i in range(m):
			diff=np.mat(xarrloc)-X[i,:]
			weights[i,i]=np.exp((diff*diff.T)/(-2.0*k**2))
		xtx=X.T*(weights*X)
		if np.linalg.det(xtx)==0.0:
			print ("matrix singularity")
			return
		else:
			wloc=xtx.I*(X.T*(weights*Y))
		yloc=xarrloc*wloc
		return yloc
#######################################
	def lwlr (self,xarr,yarr,k):
		#用于曲线的光滑化，只能不断尝试 
		m=np.shape(xarr)[0]
		yhat=np.zeros(m)
		for i in range(m):
			yhat[i]=self.__localRegression (k,xarr[i])
		#correlation coeffient compute
		coeff=self.__corrCoeff(yhat,yarr)

		#compute coeffient of determination R2
		yestimate=np.mat(yhat).T
		R2=self.__deterCoeff(yestimate,yarr)
		return (yhat,coeff,R2)
#######################################
	def mmToInches (self,mm):
		#mm transform to inches
		inches=mm*0.0393700787
		return inches
#######################################
	def plotLinearRegre (self,x1,y1,x2,y2):
		#color r,g,b...
		#linestyle,[ '-' | '--' | '-.' | ':' | 'steps' | …]
		#marker,[ '+' | ',' | '.' | '1' | '2' | '3' | '4' ]
		xscat=np.array(x1)
		yscat=np.array(y1)
		xliner=np.array(x2)
		yliner=np.array(y2)
		width=self.mmToInches(90)
		height=self.mmToInches(55)
		plt.figure(1, figsize=(width, height))
		plt.subplot(111)
		p1=plt.plot(xscat,yscat,'bo')
		p2=plt.plot(xliner,yliner)

		plt.setp(p1, linewidth=0.3)
		plt.setp(p2, color='r',linestyle='-',linewidth=2)


		plt.grid(c='k',linestyle='--',linewidth=0.3)
		plt.xlabel("time(s)")  
		plt.ylabel("value(m)") 
#		plt.xlim(0.0,  1.0)
#		plt.ylim(3,  4.75)
		plt.title("A simple plot")
		plt.legend()
		
		plt.savefig("easyplot1.png",dpi = 600,bbox_inches="tight")  

		plt.show()
#######################################
	def __ridgeRegression (self,xarr,yarr,lada):
		#岭回归的实质是在矩阵xTx上加一个ladaI从而使得矩阵非奇异
		#w=(XTX+ladaI).inv*XTy
		#input:xarr--[x0,x1,x2,..]Nx(m-1)  yarr--Nx1 需要转换为矩阵
		#output:ws--regression coeffient
		xarrmat=np.mat(xarr)
		yarrmat=np.mat(yarr)
		XTX=xarrmat.T*xarrmat
		xinter=XTX+np.eye(np.shape(xarrmat)[1])*lada
		if np.linalg.det(xinter)==0.0:
			print ("matrix singularity")
			return
		else:
			ws=xinter.I*(xarrmat.T*yarrmat)
		return ws
#######################################
	def ridgeTest (self,xarr,yarr,numTestPts):
		#ridgeRegression test
		#input:xarr--[x1,x2,..]numTestPtsx(m-1)  yarr--numTestPtsx1 需要转换为矩阵
		#numTestPts--total test number
		#output:
		xmat=np.mat(xarr)
		ymat=np.mat(yarr)
		ymean=np.mean(ymat,0)
		xmeans=np.mean(xmat,0)
		ymat=ymat-ymean
		xvar=np.var(xmat,0)
		xmat=(xmat-xmeans)/xvar
		wmat=np.zeros((numTestPts,np.shape(xmat)[1]))
		for i in range(numTestPts):
			ws=self.__ridgeRegression(xmat,ymat,np.exp(i*0.1-30))
			wmat[i,:]=ws.T
		return wmat
#######################################
	def crossValidation (self,xarr,yarr,numval=10):
		#cross validaton
		#input:xarr--[x1,x2,..]Nx(m-1)  yarr--Nx1 numval-20%times test
		m=len(yarr)
		indexlist=range(m)
		erromat=np.zeros((numval,500))
		trainErrorMat=np.zeros((numval,500))
		for i in range(numval):
			trainx=[]
			trainy=[]
			testx=[]
			testy=[]
			random.shuffle(indexlist)
			for j in range(m):
				a=indexlist[j]
				xarrr=xarr[a]
				yarrr=yarr[a]
				if j<m*0.8:
					trainx.append(list(xarrr))
					trainy.append(yarrr)
				else:
					testx.append(list(xarrr))
					testy.append(yarrr)
			testyy=np.mat(testy).T
			trainyy=np.mat(trainy).T
			wmat=self.ridgeTest(trainx,trainyy,500)
			np.savetxt("weight"+str(i+1)+".txt",wmat,fmt="%f")
			for k in range(500):
				mattestx=np.mat(testx)
				mattrainx=np.mat(trainx)
				meantrain=np.mean(mattrainx,0)
				vartrain=np.var(mattrainx,0)
				mattestx=(mattestx-meantrain)/vartrain
				mattrainx=(mattrainx-meantrain)/vartrain
				yest=mattestx*np.mat(wmat[k,:]).T+np.mean(trainyy)
				trainyest=mattrainx*np.mat(wmat[k,:]).T+np.mean(trainyy)
				erromat[i,k]=((yest.T.A-np.array(testyy))**2).sum()
				trainErrorMat[i,k]=((trainyest.T.A-np.array(trainyy))**2).sum()
		np.savetxt("erromat.txt",erromat.T/float(len(testyy)),fmt="%f")
		np.savetxt("trainerromat.txt",trainErrorMat.T,fmt="%f")
		meanerrors=np.mean(erromat,0)
		listerrors=list(meanerrors)
		minindex=listerrors.index(min(listerrors))
		print (minindex)
		minmean=float(min(meanerrors))
		bestweights=wmat[np.nonzero(meanerrors==minmean)]
		xmat=np.mat(xarr)
		ymat=np.mat(yarr).T
		meanx=np.mean(xmat,0)
		varx=np.var(xmat,0)
		unreg=bestweights/varx

		np.savetxt("bestweigthRidgreRegression.txt",unreg.T,fmt="%f")
		print ("the best model from ridge regression is:\n",unreg)
		return minindex
#######################################
	def forwardStep (self,xarr,yarr,eps,numIter):
		#前向逐步回归
		#input:input:xarr--[x1,x2,..]Nx(m-1)  yarr--Nx1 需要转换为矩阵
		#eps-误差限值  numIter-最大迭代次数
		xmat=np.mat(xarr)
		ymat=np.mat(yarr)
		ymean=np.mean(ymat,0)
		ymat=ymat-ymean
		xmeans=np.mean(xmat,0)
		xvar=np.var(xmat,0)
		xmat=(xmat-xmeans)/xvar
		m,n=np.shape(xmat)
		returnmat=np.zeros((numIter,n))
		ws=np.zeros((n,1)) 
		wstest=ws.copy()
		wsmax=ws.copy()
		for i in range(numIter):
#			print ws.T
			lowestError=np.inf
			for j in range(n):
				for sign in [-1,1]:
					wstest=ws.copy()
					wstest[j]+=eps*sign
					ytest=xmat*wstest
					rssE=((ymat.A-ytest.A)**2).sum()
					if rssE<lowestError:
						lowestError=rssE
						wsmax=wstest
			ws=wsmax.copy()
			returnmat[i,:]=ws.T
		return returnmat
#######################################
	def PCA (self,topNfeat):
		#principal component analysis,用于数据的预清洗
		matdata=np.mat(self.data)[:,:]
		meanvals=np.mean(matdata,axis=0)

		meanremoved=(matdata-meanvals)
		covmat=np.cov(meanremoved,rowvar=0)
		np.savetxt("covmat.txt",covmat,fmt="%f")
		eigvals,eigvects=np.linalg.eig(np.mat(covmat))
		eigvalind=np.argsort(eigvals)
		print ("eigvects",eigvects)
		eigvalind=eigvalind[:-(topNfeat+1):-1]
		print ("eigvalid",eigvalind)
		
		redeigvects=eigvects[:,eigvalind]
		print (redeigvects*redeigvects.T)

		lowddatamat=meanremoved*redeigvects
		reconmat=(lowddatamat*redeigvects.T)+meanvals
	

		np.savetxt('newdata.txt',lowddatamat,fmt="%.8f")
		np.savetxt('recoverdata.txt',reconmat,fmt="%.8f")

		return lowddatamat,reconmat
#######################################	
	def Smooth_5points (self,xarr,yarr,niterate):
		#
		m=len(xarr)
		yinput=yarr
		for k in range(niterate):
			yout=[]
			yout.append( (69*yinput[0] +4*(yinput[1] +yinput[3]) -6*yinput[2] -yinput[4]) /70.0)
			yout.append((2* (yinput[0] +yinput[4]) +27*yinput[1] +12*yinput[2] -8*yinput[3]) /35.0)
			for j in np.arange(2,m-2):
				yout.append((-3*(yinput[j-2] +yinput[j+2]) +12*(yinput[j-1] +yinput[j+1]) +17*yinput[j]) /35.0)
			yout.append((2*(yinput[m-1] +yinput[m-5]) +27*yinput[m-2] +12*yinput[m-3] -8*yinput[m-4]) /35.0)
			yout.append((69*yinput[m-1] +4* (yinput[m-2] +yinput[m-4]) -6*yinput[m-3] -yinput[m-5]) /70.0)
			yinput=yout
		return (xarr,yout)
#######################################
	def CrossValidationForward (self,xarr,yarr,numval,eps,numIter):
		#cross validaton
		#input:xarr--[x1,x2,..]Nx(m-1)  yarr--Nx1 numval-10%times test
		m=len(yarr)
		indexlist=range(m)
		erromat=np.zeros((10,numIter))
		for i in range(numval):
			trainx=[]
			trainy=[]
			testx=[]
			testy=[]
			random.shuffle(indexlist)
			for j in range(m):
				a=indexlist[j]
				xarrr=xarr[a]
				yarrr=yarr[a]
				if j<m*0.9:
					trainx.append(list(xarrr))
					trainy.append(yarrr)
				else:
					testx.append(list(xarrr))
					testy.append(yarrr)
			testyy=np.mat(testy).T
			trainyy=np.mat(trainy).T
			wmat=self.forwardStep (trainx,trainyy,eps,numIter)
			cwd=os.getcwd()
			pathsave=os.path.join(cwd,'bestweigh'+str(i)+".txt")
			np.savetxt(pathsave,wmat,fmt="%.8f")
			mm,nn=np.shape(wmat)
			for k in range(mm):
				mattestx=np.mat(testx)
				mattrainx=np.mat(trainx)
				meantrain=np.mean(mattrainx,0)
				vartrain=np.var(mattrainx,0)
				mattestx=(mattestx-meantrain)/vartrain
				yest=mattestx*np.mat(wmat[k,:]).T+np.mean(trainyy)
				erromat[i,k]=((yest.T.A-np.array(testyy))**2).sum()
		np.savetxt("errormat.txt",erromat/(len(testyy)),fmt="%.8f")
		d=np.amin(erromat)
		result=np.where(erromat==d)
		print (result)
		rownum=list(result[0])[0]
		columnum=list(result[1])[0]
		print (rownum,columnum)
#######################################
	def __lassoRegression (self,xarr,yarr,lambd=0.2,threshold=0.1):
		#Lasso的cost function的两部分求解
		#RSS(w)=sigma(i=1,m)(yi-sigma(j=1,n)[xijwj])^2
		#正则项的求导：diff(sigma|wj|)= -lambda, wk<0
		#                               [-lambda,lambda] wk=0
		#                               lambda          wk>0
		#
		# wkestimate=(pk+lambda/2)/zk, pk<-lambda/2
		#            0               , -lambda/2<=pk<=lambda/2
		# (pk-lambda/2)/zk			 , pk>lambda/2
		# pk=sigma(xik)(yi-sigma(j=1,n,j!=k)xijwj)
		# zk=sigma(mxik^2)
	
		rss=lambda xarr,yarr,w: (yarr-xarr*w).T*(yarr-xarr*w)
		m,n=xarr.shape
		w=np.matrix(np.zeros((n,1)))
		r=rss(xarr,yarr,w)
		niter=itertools.count(1)
		for it in niter:
			for k in range(n):
				z_k=(xarr[:,k].T*xarr[:,k])[0,0]
				p_k=0
				for i in range(m):
					p_k+=xarr[i,k]*(yarr[i,0]-sum([xarr[i,j]*w[j,0] for j in range(n) if j!=k]))
				if p_k<-lambd/2.0:
					w_k=(p_k+lambd/2.0)/float(z_k)
				if p_k>lambd/2.0:
					w_k=(p_k-lambd/2.0)/float(z_k)
				else:
					w_k=0
				w[k,0]=w_k
			r_prime=rss(xarr,yarr,w)
			delta=abs(r_prime-r)[0,0]
			r=r_prime
			print("Iteration:{},delta={}".format(it,delta))
			if delta<threshold:
				break
		return w
#######################################
	def lassoTest (self,xarr,yarr,numTestPts):
		xmeans=np.mean(xarr,0)
		ymean=np.mean(yarr,0)
		xmat=xarr-xmeans
		ymat=yarr-ymean
		xvar=np.var(xmat,0)
		yvar=np.var(ymat,0)
		xinput=(xarr-xmeans)/xvar
		yinput=(yarr-ymean)

		wmat=np.zeros((numTestPts,np.shape(xmat)[1]))
		for i in range(numTestPts):
			ws=self.__lassoRegression (xmat,ymat,lambd=np.exp(i-10))
			wmat[i,:]=ws.T
		np.savetxt("allWeights.txt",wmat,fmt="%f")
		return wmat
#######################################
	def lassoCrossValidation (self,xarr,yarr,numval=10):
		m=len(yarr)
		indexlist=range(m)
		erromat=np.zeros((numval,40))
		for i in range(numval):
			trainx=[]
			trainy=[]
			testx=[]
			testy=[]
			random.shuffle(indexlist)
			for j in range(m):
				a=indexlist[j]
				xarrr=xarr[a]
				yarrr=yarr[a]
				if j<m*0.8:
					trainx.append(list(xarrr))
					trainy.append(yarrr)
				else:
					testx.append(list(xarrr))
					testy.append(yarrr)
			testyy=np.mat(testy).T
			trainyy=np.mat(trainy).T
			wmat=self.lassoTest (np.mat(trainx),trainyy,40)
			for k in range(40):
				mattestx=np.mat(testx)
				mattrainx=np.mat(trainx)
				meantrain=np.mean(mattrainx,0)
				vartrain=np.var(mattrainx,0)
				mattestx=(mattestx-meantrain)/vartrain
				yest=mattestx*np.mat(wmat[k,:]).T+np.mean(trainyy)
				erromat[i,k]=((yest.T.A-np.array(testyy))**2).sum()
		meanerrors=np.mean(erromat,0)
		listerrors=list(meanerrors)
		minindex=listerrors.index(min(listerrors))
		print (minindex)
		minmean=float(min(meanerrors))
		bestweights=wmat[np.nonzero(meanerrors==minmean)]
		xmat=np.mat(xarr)
		ymat=np.mat(yarr).T
		meanx=np.mean(xmat,0)
		varx=np.var(xmat,0)
		unreg=bestweights/varx
		print ("the best model from ridge regression is:\n",unreg)
		return unreg
		




######################################
if __name__=='__main__':

	#Linear regression test
#	data=np.loadtxt("ex0.txt")
#	instance=Regression(data)
#	wtot=instance.linearRegression()
#	w=wtot[0]
#	corrcoef=wtot[1]
#	R2=wtot[2][0]
#	print w
#	print "corrcoef=:",corrcoef
#	print "R2=:",R2
#	a=np.mat(data)
#	x=a[:,1]
#	y=a[:,2]
#	datacopy=data.copy()
#	xCopy=np.mat(datacopy[:,0:2])
#	xCopy.sort(0)
#	yliner=np.mat(xCopy)*np.mat(w)
#	instance. plotLinearRegre (x,y,xCopy[:,1],yliner)

############################################################	
#	local regression test
#	data=np.loadtxt("smooth.txt")
#	instance=Regression(data)
#	a=np.mat(data)
#	x=a[:,0:2]
#	y=a[:,2]
#	result=instance.lwlr(x,y,1)
#	yhat=result[0] 
#	print "corrcoef=:",result[1]
#	print "R2=:",result[2][0]
#	x1=a[:,1]
#	y1=a[:,2]
#	x2=[]
#	y2=[]
#	x2index=np.array(x1).argsort(0)
#	for i in range(len(x2index)):
#		index=(x2index[i])[0]
#		x2.append(np.array(x1)[index][0])
#		y2.append(yhat[index])
#	instance. plotLinearRegre (x1,y1,x2,y2)

#############################################################
#	ridge regression test，开始部分的特征大小表明了主要特征
#	data=np.loadtxt("LogData/1.txt")
#	instance=Regression(data)
#	xarr=data[:,:-2]
#	yarr=np.mat(data[:,-2]).T
#	ws=instance.ridgeTest(xarr,yarr,30)
#	print np.shape(ws)
#	np.savetxt("wsridge.txt",ws,fmt="%f")
#	x=np.arange(-10,20,1)
#	width=instance.mmToInches(90)
#	height=instance.mmToInches(55)
#	plt.figure(1, figsize=(width,height))
#	plt.subplot(111)
#	p1=plt.plot(x,ws)
##	plt.legend(('1','2','3','4','5','6','7','8'))
#	plt.show()
##############################################################
# forward stape regression test,根据最终特征的大小确定主要特征
#	data=np.loadtxt("LogData/1.txt")
#	instance=Regression(data)
#	xarr=data[:,:-2]
#	yarr=np.mat(data[:,-2]).T
#
#	result=instance.forwardStep(xarr,yarr,0.0001,2000)
#	np.savetxt("ForwardStepWiseCoeffient/1.txt",result,fmt="%f")
###############################################################
####	cross validation
#	cwd=os.getcwd()
#	pathsave=os.path.join(cwd,'RidgeCoeffient/',"RidgeCoeffient.txt")
#	ridgeresult=[]
#		
#	for i1 in range(41):
#		path1=os.path.join(cwd,'LogData/',str(i1+1)+".txt")
#		data=np.loadtxt(path1)
#		instance=Regression(data)
#		xarr=data[:,:-2]
#		yarr=np.array(data[:,-2]).T
#		result=instance.crossValidation(xarr,yarr)
#		ridgeresult.append(np.array(result).T)
#	np.savetxt(pathsave,(ridgeresult),fmt="%.8f")

###############################################################
	#pca,lowdmat是降维后的矩阵   reconmat-根据主特征数重构数据
#	data=np.loadtxt("testSet.txt")
#	instance=Regression(data)
#	lowdmat,reconmat=instance.PCA(1)
#	plt.figure(1, figsize=(90*0.04, 55*0.04))
#	print lowdmat
#	print reconmat
#	plt.subplot(111)
#	plt.plot(data[:,0],data[:,1],"b^")
#	plt.plot(reconmat[:,0],reconmat[:,1],"ro")
##	plt.plot(lowdmat[:,0],lowdmat[:,1],"gs")
#	plt.savefig("pca.eps",dpi = 600,bbox_inches="tight")
#	plt.show()
##############################################################
	#smooth_5points
#	data=np.loadtxt("smooth.txt")
#	instance=Regression(data)
#	x=data[:,1]
#	y=data[:,2]
#	result=instance.Smooth_5points (x,y,20000)
#	instance.plotLinearRegre (x,y,result[0],result[1])
	
###################crossvalidationForward###########################

#	data=np.loadtxt("LogData/41.txt")
#	instance=Regression(data)
#	xarr=data[:,:-2]
#	yarr=np.array(data[:,-1]).T
#	instance.CrossValidationForward (xarr,yarr,10,0.0001,2000)

###################lassoRegressionTest##############################
	loadData=np.loadtxt("testSet.txt")
	instance=Regression(loadData)
	xarr=np.mat(loadData[:,:-1])
	yarr=np.mat(loadData[:,-1]).T
	lambd=0.2
	instance.lassoTest (xarr,yarr,lambd)