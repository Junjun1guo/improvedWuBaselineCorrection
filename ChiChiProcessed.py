#-*-coding: UTF-8-*-
import numpy as np
import os
import math

# 1gal=0.001019716213 g
factor=0.001019716213
if not os.path.exists('ChiChiEarthquakeAccg/'):
	os.makedirs('ChiChiEarthquakeAccg')
pathDir =  os.listdir('ChiChiGroundMotionRawAcc/')
fileName=[]
for allDir in pathDir:
	child =  os.path.splitext(allDir)
	fileName.append(child[0])
deltat=[]


for i1 in range(len(fileName)):
	cwd=os.getcwd()
	pathall=os.path.join(cwd,'ChiChiGroundMotionRawAcc/',str(fileName[i1])+".txt")
	txtopen=np.loadtxt(pathall)
	deltat.append(txtopen[:,0][1]-txtopen[:,0][0])
	Uper=[x*factor for x in txtopen[:,1]]
	North=[x*factor for x in txtopen[:,2]]
	East=[x*factor for x in txtopen[:,3]]

	pathsave=os.path.join(cwd,'ChiChiEarthquakeAccg/U',str(fileName[i1])+".txt")
	np.savetxt(pathsave,Uper,fmt="%f")
#
	pathsave=os.path.join(cwd,'ChiChiEarthquakeAccg/N',str(fileName[i1])+".txt")
	np.savetxt(pathsave,North,fmt="%f")
#
	pathsave=os.path.join(cwd,'ChiChiEarthquakeAccg/E',str(fileName[i1])+".txt")
	np.savetxt(pathsave,East,fmt="%f")

pathsaveT=os.path.join(cwd,'ChiChiEarthquakeAccg/',"deltaT.txt")
hZtransName=np.mat(fileName).T
hZtransDelta=np.mat(deltat).T
aa=np.hstack((hZtransName,hZtransDelta))

np.savetxt(pathsaveT,aa,fmt="%s")

	






		



	