#-*-coding: UTF-8-*-
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from LHSamples import Sample
import os
############################################
def dispPulse (t,Dp,t1,Tp):
	"""
	Baker et al. displacement pulse point value (cm)
	:param t: time point(s)
	:param Dp: permanent displacement(cm)
	:param t1:pulse start time(s)
	:param Tp:pulse period (s)
	:return:discrete displacement value(cm)
	"""
	discretePulse=Dp/float(2)*np.sin(math.pi*(t-t1-Tp/float(2))/float(Tp))+Dp/float(2)
	return discretePulse
############################################
def velPulse (t,Dp,t1,Tp):
	"""
	Baker et al. velocity pulse point value (cm/s)
	:param t: time point(s)
	:param Dp: permanent displacement(cm)
	:param t1: pulse start time(s)
	:param Tp: pulse period (s)
	:return: discrete velocity value(cm/s)
	"""
	discretePulse=(Dp*math.pi)/(Tp*float(2))*np.cos(math.pi*(t-t1-Tp/float(2))/float(Tp))
	return discretePulse
############################################acc(g)
def accPulse (t,Dp,t1,Tp):
	"""
	Baker et al. acceleration pulse point value (g)
	:param t: time point(s)
	:param Dp: permanent displacement(cm)
	:param t1: pulse start time(s)
	:param Tp: pulse period (s)
	:return: discrete acceleration value(g)
	"""
	discretePulse=-(Dp*math.pi**2)/(float(2)*Tp**2)*np.sin(math.pi*(t-t1-Tp/float(2))/float(Tp))
	return discretePulse/float(981)
############################################
def iterateFunction (Tp,t1,dt,disp):
	"""
	Calculate error between the original displaceemnt and the displacement pulse model
	:param Tp: pulse period (s)
	:param t1: pulse start time(s)
	:param dt: time interval (s)
	:param disp: original displacement time history(cm)
	:return: deltaError-mean error
			 Tp-pulse period (s)
			 t1=pulse start time(s)
			 Dp-permanent displacement(cm)
	"""
	startNumber=int(t1/float(dt))
	endNumber=int((t1+Tp)/float(dt))
	intervalT=[x*dt for x in range(startNumber,endNumber)]
	pulseMotion=[]
	counterMotion=disp[startNumber:endNumber]
	Dp=sum(disp[-50:])/float(50)
	for each  in intervalT:
		discretePoint=dispPulse(each,Dp,t1,Tp)
		pulseMotion.append(discretePoint)
	deltaError=sum([(x-y)**2 for x,y in zip(counterMotion,pulseMotion)])/float(endNumber-startNumber)
	return deltaError,Tp,t1,Dp
####################################################################
#####################---main program---#############################
####################################################################
fileName="TCU075.txt"
disp=np.loadtxt('dispBaselineCorre/E/'+fileName)
vel=np.loadtxt('velBaselineCorre/E/'+fileName)
acc=np.loadtxt('accBaselineCorre/E/'+fileName)
dt=0.005
numPoint=len(disp)
errorList=[10**10]
TpList=[0]
t1List=[0]
DpList=[0]
nIterate=2000
bounds = [(0.01,10),(20,33)]
instance = Sample(bounds, nIterate)
Tp = instance.LHSample()[:,0]
t1 = instance.LHSample()[:,1]
for i1 in range(nIterate):
	errorTerm,TpTerm,t1Term,DpTerm=iterateFunction(Tp[i1],t1[i1],dt,disp)
	if errorTerm<errorList[0]:
		errorList[0]=errorTerm
		TpList[0]=TpTerm
		t1List[0]=t1Term
		DpList[0]=DpTerm
print ("Tp:",TpList[0])
print ("t1:",t1List[0])
print ("Dp:",DpList[0])
####################################################################
startNumber=int(t1List[0]/float(dt))
endNumber=int((t1List[0]+TpList[0])/float(dt))
intervalT=[x*dt for x in range(startNumber,endNumber)]
counterMotion=disp[startNumber:endNumber]
dispMotion=[]
velMotion=[]
accMotion=[]
for each  in intervalT:
		discretePoint=dispPulse(each,DpList[0],t1List[0],TpList[0])
		dispMotion.append(discretePoint)
		velDiscretePoint=velPulse (each,DpList[0],t1List[0],TpList[0])
		velMotion.append(velDiscretePoint)
		accDiscretePoint=accPulse(each,DpList[0],t1List[0],TpList[0])
		accMotion.append(accDiscretePoint)
####################################################################
revisedStart=[0 for x in range(startNumber)]
revisedEnd=[dispMotion[-1] for x in range(endNumber,numPoint)]
revisedEndVel=[0 for x in range(endNumber,numPoint)]
revisedEndAcc=[0 for x in range(endNumber,numPoint)]
revisedDisp=revisedStart+dispMotion+revisedEnd
revisedVel=revisedStart+velMotion+revisedEndVel
revisedAcc=revisedStart+accMotion+revisedEndAcc
residualDisp=[num2-num3 for num2,num3 in zip(disp,revisedDisp)]
residualVel=[num2-num3 for num2,num3 in zip(vel,revisedVel)]
residualAcc=[num2-num3 for num2,num3 in zip(acc,revisedAcc)]
####################################################################
np.savetxt("accPulse/E/"+fileName,revisedAcc,fmt="%f")
np.savetxt("velPulse/E/"+fileName,revisedVel,fmt="%f")
np.savetxt("dispPulse/E/"+fileName,revisedDisp,fmt="%f")
np.savetxt("accResidual/E/"+fileName,residualAcc,fmt="%f")
np.savetxt("velResidual/E/"+fileName,residualVel,fmt="%f")
np.savetxt("dispResidual/E/"+fileName,residualDisp,fmt="%f")
####################################################################
times=[dt*num1 for num1 in range(numPoint)]
plt.subplot(331)
# plt.plot(times,acc,"-k")
plt.plot(times,revisedAcc,"--r")

plt.subplot(334)
plt.plot(times,vel,"-k")
plt.plot(times,revisedVel,"--r")

plt.subplot(337)
plt.plot(times,disp,"-k")
plt.plot(times,revisedDisp,"--r")

plt.show()







	

		










	


	
	
	
		
	
	



	
	
		
	



	
	