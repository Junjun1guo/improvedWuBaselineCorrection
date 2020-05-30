#coding=utf-8
import numpy as np
import matplotlib.pyplot as pl

class Sample():
	def __init__(self,bounds,N):
		"""
		initial data
		bounds: variable bounds, such as bounds=[(10,21),(21,34)]
		N:sample numbers
		"""
		self.bounds=bounds
		self.N=N
		self.D=len(self.bounds)
	
	def LHSample( self):
		result = np.empty([self.N, self.D])
		temp = np.empty([self.N])
		d = 1.0 / self.N

		for i in range(self.D):

			for j in range(self.N):
				temp[j] = np.random.uniform(
					low=j * d, high=(j + 1) * d, size = 1)[0]

			np.random.shuffle(temp)

			for j in range(self.N):
				result[j, i] = temp[j]

		b = np.array(self.bounds)
		lower_bounds = b[:,0]
		upper_bounds = b[:,1]
		if np.any(lower_bounds > upper_bounds):
			print ("The parameters low bound value should smaller than the upper bound!")
			return 

		#   sample * (upper_bound - lower_bound) + lower_bound
		np.add(np.multiply(result,
                       (upper_bounds - lower_bounds),
                       out=result),
           lower_bounds,
           out=result)
		
		return result

#if __name__ =='__main__':
##
###	N = 100
###	bounds=[(20,90),(5,30),(2,6),(1,3),(1,2),\
###			(317400,372600),(23370,30230),(25,50),
###			(0.025,0.1),(0.15,0.3),(20,150)]
###	instance=Sample(bounds,N)
###	results=instance.LHSample()
###	np.savetxt("SampleResults.txt",results,fmt="%f")
###	XY = np.array(results)
###	X = XY[:,0]
###	Y = XY[:,1]
###	pl.scatter(X,Y)
###	pl.show()
#	N =2000
#	bounds=[(0.01,100)]
#	instance=Sample(bounds,N)
#	results=instance.LHSample()
#	np.savetxt("beta.txt",results,fmt="%f")
#	XY = np.array(results)
#	X = XY[:,0]
#	Y = XY[:,1]
#	pl.scatter(X,Y)
#	pl.show()