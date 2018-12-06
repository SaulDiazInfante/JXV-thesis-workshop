'''
	This class solve implicitly and explicitly the SIS reduce equation proposal in:
	A Stochastic Differential Equation SIS Epidemic Model
	
	
	\begin{equation}
		dI(t) = I(t)
			\left[
				\left(
					\beta N - \mu -\sigma - \beta I(t) 
				\right)
				dt
				+ \sigma (N -I(t)) dB(t)
			\right] 
	\end{equation}
	Gray, A Greenhalgh, D,Hu, L,Mao, X, Pan, J
'''

import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint
class StochastisSIS:
	'''
		Set the parameters and terms of the Stochastic Duffin Van der Pol equation.
	'''
	def InitializeMesh(self, k, p, r, T0, T):
		'''
			Set setencil parameters
		'''
#Stensil of the mesh
		self.k = k
		self.p = p
		self.r = r
		self.T0 = T0
		self.N=2*10.0**k
		self.P=10.0**p
		self.R=10.0**r
		#self.N = 2.0 ** k
		#self.P = 2.0 ** p
		#self.R = 2.0 ** r
		self.T = T
#
		self.dt = self.T / np.float(self.N)
		self.IndexN = np.arange(self.N + 1)
#set of index to Ito integral
		self.tau = self.IndexN[0:self. N + 1:self.P]
		self.t = np.linspace(0, self.T, self.N + 1)
#
		self.Dt = np.float(self.R) * self.dt
		self.L = self.N / self.R
#diffusion part
		self.DistNormal1 = np.random.randn(np.int(self.N))
		self.DistNormal2 = np.random.randn(np.int(self.N))
		self.dW = np.sqrt(self.dt) * self.DistNormal1
		self.dZ = 0.5 * (self.dt ** 1.5) * (self.DistNormal1 + (3.0 ** (-.5)) * self.DistNormal2)
		self.Z = np.cumsum(self.dZ)
		self.Z = np.concatenate(([0], self.Z))
		self.W = np.cumsum(self.dW)
		self.W = np.concatenate(([0], self.W))

	def SetParametersStochasticSIS(self, beta, mu ,gamma, sigma, M, Xzero):
		'''
			Set parameters of SDE Duffin Van Der Pol.
		'''
		self.beta = beta
		self.mu = mu 
		self.gamma = gamma
		self.sigma = sigma
		self.M = M
		self.Xzero = Xzero

	def a(self, X):
		'''
			The drifft term of the SDE.
		'''
		beta = self.beta
		gamma = self.gamma
		mu = self.mu
		M = self.M
		return X * (beta * M -mu - gamma - beta * X) 

	def b(self, X):
		'''
			The diffusion term.
		'''
		sigma = self.sigma
		M = self.M
		return sigma * (M - X) * X

class NumericsStochasticSIS(StochastisSIS):
	'''
		Numerics for the Stochastic SIS Equation.
	'''
	def EM(self, flag=0):
#%  Preallocate Xem for efficiency.
		Dt =self.Dt
		L = self.L
		R = self.R
		
		if flag == 1:
			Dt = self.dt
			L = self. N
			R =1
		self.Xem = np.zeros(L + 1)
		self.Xem[0] = self.Xzero
		for j in np.arange(L):
			self.Winc = np.sum(self.dW[R * (j):R * (j + 1)])
			increment = self.Xem[j] + Dt * self.a(self.Xem[j]) + self.Winc * self.b(self.Xem[j])
			self.Xem[j + 1] = increment
		xem = self.Xem
		return xem

	def SteklovMatus(self):
		h = self.Dt
		self.Xstkm = np.zeros([self.L + 1])
		self.Xstkm[0] = self.Xzero
		beta = self.beta
		gamma = self.gamma
		mu = self.mu
		M = self.M
		a = beta * M - (mu + gamma)
		b = -beta		
		for j in np.arange(self.L):
			self.Winc =np.sum(self.dW[self.R*(j):self.R*(j+1)])
			xj = self.Xstkm[j]
			#Ahx = (np.exp(a  *h) * xj)/(a + b * xj)
			num = a * np.exp(a * h) * xj
			den = a + b * xj * (1.0 - np.exp(a * h))
			self.Xstkm[j + 1] = num / den + self.b(xj) * self.Winc
		xstkm = self.Xstkm
		return xstkm
	def fdet(self, y,t):
		beta = self.beta
		gamma = self.gamma
		mu = self.mu
		M = self.M
		I = y * (beta * M -mu - gamma - beta * y)
		return I
	def DetSol(self):
		y0 = self.Xzero
		t = self.t
		self.SolDet = odeint(self.fdet, y0, t)
		solDet = self.SolDet
		return solDet
	def SaveData(self):
		'''
			Method to save the numerical solutions and parameters of Van der Pol ODE.
		'''
		#t=self.t[0:-1:self.R].reshape([self.t[0:-1:self.R].shape[0],1])
		t = self.dt * self.tau
		n =self.Xem.shape[0]
		Uem = self.Xem[n-10001:-1]
		Ustkm = self.Xstkm[n-10001:-1]
		tagPar = np.array([
		'k=',
		'r=',
		'T0=',
		'N=',
		'R=',
		'T=',
		'dt=',
		'Dt=',
		'L=',
		'M',
		'beta=',
		'mu=',
		'gamma',
		'x0='
		])
		ParValues = np.array([
		self.k,
		self.r,
		self.T0,
		self.N,
		self.R,
		self.T,
		self.dt,
		self.Dt,
		self.L,
		self.M,
		self.beta,
		self.mu,
		self.gamma,
		self.Xzero
	])
		strPrefix = str(self.Dt)
		name1 = 'Parameters' + strPrefix + '.txt'
		name2 = 'Solution'+str(self.sigma) + strPrefix + '.txt'
		PARAMETERS = np.column_stack((tagPar, ParValues))
		np.savetxt(name1, PARAMETERS, delimiter=" ", fmt="%s")
		np.savetxt(name2,
			np.transpose(
				(
					t[n-10001:-1], Uem, Ustkm
				)
			), fmt='%1.12f', delimiter='\t')
