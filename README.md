# improvedWuBaselineCorrection
conduct ground motion baseline correction based on Wu et al. [2007] method

Wu Y-M, Wu C-F. Approximate recovery of coseismic deformation from Taiwan strong-motion records. Journal of Seismology. 2007;11(2):159-70.

basic example:
```python
###provide the acceleration, velocity and displacement paths of the unprocessed motion
accFilePath='ChiChiEarthquakeAccg/N'
velFilePath='ChiChiEarthquakeVel/N'
dispFilePath='ChiChiEarthquakeDisp/N'
###provide the save paths for the processed acceleration, velocity and displacement
saveAccPath='accBaselineCorre/N'
saveVelPath='velBaselineCorre/N'
saveDispPath='dispBaselineCorre/N'
dt=0.005 #time interval (s)
nIterate=500 # sample size for T2 position from T3 to the end
fileNamei='TCU068' #file name of unprocessed motion
T1=6050 #T1 position in the motion, if T1=None the program will automatically determine T1
T3=10000 # T3 position in the motion
###call the improved Wu et al. method to conduct baseline correction
improvedMethod (accFilePath,velFilePath,dispFilePath,dt,fileNamei,nIterate,saveAccPath,saveVelPath,saveDispPath,T3,T1)
```
