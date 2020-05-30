# improvedWuBaselineCorrection
conduct ground motion baseline correction based on Wu et al. [2007] method

Wu Y-M, Wu C-F. Approximate recovery of coseismic deformation from Taiwan strong-motion records. Journal of Seismology. 2007;11(2):159-70.
## basic example:
```python
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


```
