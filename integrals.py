import numpy as np
import scipy as sp




	
def Int(l,m,lmb,k,f):
	H = np.zeros(shape = (l+1,abs(m)+2),dtype = 'complex64')
	H[:,0] = 1
	for i in range(l+1):
		H[i,-1] = lmb - f*1j+1.0/(i+1)
	
	H[0,1] = 0.5*(1-np.exp(1j*np.log((lmb**2+(k-f)**2)/(lmb**2+(k+f)**2))/(2*k)-(np.arctan((k-f)/lmb)+np.arctan((k+f)/lmb))/k))
	
	for j in range(1,l+1):
		H[j,1] = j*(j+0.5)/(1.0+k*k*j*j)*( (k*k+(lmb - f*1j)**2)*H[j-1,1] - 2*(lmb - f*1j)*H[j-1,0] + H[j-1,-1])
	
	for i in range(2,m+1):
		H[l,i] = (2*H[l,i-1]*((-1+(lmb-f*1j)*(i-2-l))/(i-1)) + H[l,i-2]*(2*l+3-i)/(i-1))/(k**2+(lmb-f*1j)**2)
	
	return H[l,m]
