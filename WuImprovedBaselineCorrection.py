#-*-coding: UTF-8-*-
####################################################################################
"""
successfully executed in python 3.6
"""
####################################################################################
import numpy as np
import os
import math
import copy
from Statistic import *
from LHSamples import Sample
###############################################################################
def accToVelocity (dt,acc):
	#convert acceleration(cm/s2) to velocity(cm/s)
	vel=[0]
	num=len(acc)
	for i in range(num-1):
		velocity=((acc[i]+acc[i+1])*dt/float(2))+vel[-1]
		vel.append(velocity)
	return vel
##############################################################################
def velToDisp (dt,vel):
	#convert velocity(cm/s) to displacement(cm)
	disp=[0]
	for i in range(len(vel)-1):
		displacement=(vel[i]+vel[i+1])*dt/float(2)+disp[-1]
		disp.append(displacement)
	return disp
#######################################################################
def improvedMethod (accFilePath,velFilePath,dispFilePath,t,fileNamei,nIterate,saveAccPath,\
					saveVelPath,saveDispPath,T3,T1Self=None,T2=None):
	"""
	improved Wu et al method for basedline correction
	Wu Y-M, Wu C-F. Approximate recovery of coseismic deformation from Taiwan strong-motion records.
	 Journal of Seismology. 2007;11(2):159-70.
	:param accFilePath: the file path of acceleration
	:param velFilePath: the file path of velocity
	:param dispFilePath: the file path of displacement
	:param t: time interval of motion (s)
	:param fileNamei: fileName of the processed ground motion
	:param nIterate: sample numbers for t2 values
	:param saveAccPath: the save path of processed acceleration
	:param saveVelPath: the save path of processed velocity
	:param saveDispPath: the save path of processed displacement
	:param T3: T3 position in the motion
	:param T1Self: T1 position in the motion, if T1self is none,the program will automatically determine it
	:return: None
	"""
	cwd=os.getcwd()
	pathAccE=os.path.join(cwd,accFilePath,str(fileNamei)+".txt")
	txtopenAccE=np.loadtxt(pathAccE)
	copyAccE1 = copy.deepcopy(txtopenAccE)
	T3=int(T3)
	for i2 in range(len(txtopenAccE)):
		if copyAccE1[i2] * 981 > 50:
			T150 = i2
			break
	pga=max(abs(txtopenAccE))*9.81*100 # convert acceleration to cm/s2
	
	if pga>60: # only pga>60 cm/s2 it needs baseline correction
		cwd=os.getcwd()
		pathDispE=os.path.join(cwd,dispFilePath,str(fileNamei)+".txt") #upload displacement time history to process
		txtopenDispE=np.loadtxt(pathDispE).tolist()
		lengthTxt=len(txtopenDispE)
		copyDispE=copy.deepcopy(txtopenDispE)
		reversedDispE=copyDispE.reverse()
		# automatically determine T1 position
		for i2 in range(lengthTxt):
			if txtopenDispE[-1]>0:
				if copyDispE[i2]<=0:
					T1=(lengthTxt-i2-1)
					break
			else:
				if copyDispE[i2]>=0:
					T1=(lengthTxt-i2-1)
					break
		# if T1 position larger than that determined by Iwan etal (1985), then use Iwan's T1 position
		if T150>T1:
			T1=T150
		else:
			T1=T1
		# if users provide T1, and use it in the following process
		if T1Self!=None:
			T1=int(T1Self)
		cwd=os.getcwd()
		pathVelE=os.path.join(cwd,velFilePath,str(fileNamei)+".txt")
		txtopenVelE=np.loadtxt(pathVelE).tolist()
		v0=[]
		af=[]
		fValue=[]
		T22=[]
		maxfValue=None
		# randomly generage nIterate intergers between T3 and (lengthTxt-10)
		bounds = [(T3,lengthTxt-10)]
		instance = Sample(bounds, nIterate)
		samples = list(instance.LHSample()[:,0])
		if T2!=None:
			T2Index=T2
		else:
			for i3 in range(nIterate):
				# print(i3)
				T2=int(samples[i3])
				X0 = [1 for x in range(T2, lengthTxt)]
				X1 = [x * t for x in range(T2, lengthTxt)]
				Y = txtopenVelE[T2:lengthTxt]
				hZX11 = np.mat(X0).T
				hZX33 = np.mat(X1).T
				hZY33 = np.mat(Y).T
				velData = np.hstack((hZX11, hZX33, hZY33))
				instanceVel = Regression(velData)
				wtotvel = instanceVel.linearRegression()
				wvel = wtotvel[0].tolist()
				v0vel = wvel[0][0]
				Afvel = wvel[1][0]
				# y=v0vel+Afvel*t
				Amvel = (v0vel + Afvel * T2 * t) / float((T2 - T1) * t)
				accOrignal = [x * 981 for x in txtopenAccE]
				for i4 in range(T1, T2):
					accOrignal[i4] = accOrignal[i4] - Amvel
				for i5 in range(T2, len(accOrignal)):
					accOrignal[i5] = accOrignal[i5] - Afvel
				velBaseline = accToVelocity(t, accOrignal)
				dispBaseline = velToDisp(t, velBaseline)

				X10 = [1 for x in range(T3, lengthTxt)]
				X11 = [x * t for x in range(T3, lengthTxt)]
				Y11=dispBaseline[T3:lengthTxt]
				hZX1=np.mat(X10).T
				hZX2=np.mat(X11).T
				hZY=np.mat(Y11).T
				linerData=np.hstack((hZX1,hZX2,hZY))
				#call linear regression function to fit the displacement from T3 to end
				instance=Regression(linerData)
				wtot=instance.linearRegression()
				w=wtot[0].tolist() #regression coefficients
				corrcoef=wtot[1] #correlation coefficient
				var=wtot[4] #standard deviation
				v00=w[0][0] #constant of the regression line
				aff=w[1][0] #slope of the regression line
				r=corrcoef
				b=abs(aff)
				fvalue=abs(r)/float(b*var**2) #calculate f value
				fValue.append(fvalue)
				v0.append(v00)
				af.append(aff)
				T22.append(T2)
			maxIndex=fValue.index(max(fValue)) #find the maximum f index
			T2Index=T22[maxIndex] #obtain the optimized T2 position
			maxfValue=max(fValue)
		#conduct baseline correction based on the optimized T2
		X31=[1 for x in range(T2Index,lengthTxt)]
		X33=[x*t for x in range(T2Index,lengthTxt)]
		Y3=txtopenVelE[T2Index:lengthTxt]
		hZX11=np.mat(X31).T
		hZX33=np.mat(X33).T
		hZY33=np.mat(Y3).T
		velData=np.hstack((hZX11,hZX33,hZY33))
		instanceVel=Regression(velData)
		wtotvel=instanceVel.linearRegression()
		wvel=wtotvel[0].tolist()
		v0vel=wvel[0][0]
		Afvel=wvel[1][0]
		#y=v0vel+Afvel*t
		Amvel=(v0vel+Afvel*T2Index*t)/float((T2Index-T1)*t)
		accOrignal=[x*981 for x in txtopenAccE]
		for i4 in range(T1,T2Index):
			accOrignal[i4]=accOrignal[i4]-Amvel
		for i5 in range(T2Index,len(accOrignal)):
			accOrignal[i5]=accOrignal[i5]-Afvel
		velBaseline=accToVelocity (t,accOrignal)
		dispBaseline=velToDisp (t,velBaseline)
		accToG=[x/float(981) for x in accOrignal]
		cwd=os.getcwd()
		# save the baseline corrected motion
		pathaccBaseCorreE=os.path.join(cwd,saveAccPath,str(fileNamei)+".txt")
		np.savetxt(pathaccBaseCorreE,accToG,fmt="%f")

		pathvelBaseCorreE=os.path.join(cwd,saveVelPath,str(fileNamei)+".txt")
		np.savetxt(pathvelBaseCorreE,velBaseline,fmt="%f")

		pathdispBaseCorreE=os.path.join(cwd,saveDispPath,str(fileNamei)+".txt")
		np.savetxt(pathdispBaseCorreE,dispBaseline,fmt="%f")
		return T1,T2Index,T3,maxfValue
#######################################################################
#########################---main program---############################
#######################################################################
###provide the acceleration, velocity and displacement paths of the unprocessed motion
accFilePath='ChiChiEarthquakeAccg/E'
velFilePath='ChiChiEarthquakeVel/E'
dispFilePath='ChiChiEarthquakeDisp/E'
###provide the save paths for the processed acceleration, velocity and displacement
saveAccPath='accBaselineCorre/E'
saveVelPath='velBaselineCorre/E'
saveDispPath='dispBaselineCorre/E'
dt=0.005 #time interval (s)
nIterate=100 # sample size for T2 position from T3 to the end
fileNamei='TCU068' #file name of unprocessed motion
#########################################################################
#########################################################################
##automatically determine T1 and T3,T1=(4500,5500),T3=(5000,7000)
bounds = [(6000,7000),(7000,9000)]
NIter=10 #iterate number for T1 and T3
instance = Sample(bounds, NIter)
samples =instance.LHSample()
T1sample=samples[:,0]
T3sample=samples[:,1]
T1List=[]
T2List=[]
T3List=[]
fvalueList=[]
for j1 in range(10):
	print(j1)
	###call the improved Wu et al. method to conduct baseline correction
	T11,T22,T33,fvalue=improvedMethod (accFilePath,velFilePath,dispFilePath,dt,\
									   fileNamei,nIterate,saveAccPath,saveVelPath,saveDispPath,T3sample[j1],T1sample[j1])
	T1List.append(T11)
	T2List.append(T22)
	T3List.append(T33)
	fvalueList.append(fvalue)
maxIndex=fvalueList.index(max(fvalueList))
finalT1=T1List[maxIndex]
finalT2=T2List[maxIndex]
finalT3=T3List[maxIndex]
print("finalT1,T2,T3",finalT1,finalT2,finalT3)
#########################################################################
#########################################################################
T1=finalT1 #T1 position in the motion, if T1=None the program will automatically determine T1
T3=finalT3 # T3 position in the motion
T2=finalT2 # T2 position in the motion
T11,T22,T33,fvalue=improvedMethod (accFilePath,velFilePath,dispFilePath,dt,\
									   fileNamei,nIterate,saveAccPath,saveVelPath,saveDispPath,T3,T1,T2)

