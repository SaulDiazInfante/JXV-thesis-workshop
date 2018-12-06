import numpy as np
import matplotlib.pyplot as plt
from StochasticSIS import NumericsStochasticSIS
# Stochastic SIS Parameters Fig1

beta = 0.5
mu = 20.0
gamma = 25.0
sigma = np.array([0.035,0.08,.03,.01])
M = 100.0
Xzero = np.array([90.0, 1.0])
k = 6
p = 1
r = p
T0 = 0.0
T = 10.0
#kk = np.arange(11, 21)
#
#
'''
	w=1.0
	fig_width_pt =256.075*w					# Get this from LaTeX using \showthe\columnwidth
	inches_per_pt = 1.0/72.27				# Convert pt to inch
	golden_mean = (np.sqrt(5) - 1.0) / 2.0		# Aesthetic ratio
	fig_width = fig_width_pt * inches_per_pt	# width in inches
	fig_height = fig_width * golden_mean		# height in inches
	fig_size =  [fig_width, fig_height]
	params = {'backend': 'ps',
			'axes.labelsize': 10,
			'text.fontsize': 10,
			'legend.fontsize': 8,
			'xtick.labelsize': 8,
			'ytick.labelsize': 8,
			'text.usetex': True,
			'figure.figsize': fig_size}
	plt.rcParams.update(params)
'''
#
for j in np.arange(sigma.shape[0]):
	StoSIS = NumericsStochasticSIS()
	StoSIS.InitializeMesh(k, p, r, T0, T)
	#Xzero = 90
	StoSIS.SetParametersStochasticSIS(beta, mu, gamma, sigma[j], M, Xzero[0])
	h = StoSIS.Dt
	strEl = 'Euler h=' + str(h)
	strSt = 'SSLS h =' + str(h)
	xem = StoSIS.EM(flag=0)
	xst = StoSIS.SteklovMatus()
	xdet = StoSIS.DetSol()
	StoSIS.SaveData()
	strFileName ='IStochastiSIS' + str(j)  + '.eps'
#
	fig,(ax0, ax1) = plt.subplots(ncols=2, figsize=(10, 6))
	#plt.axes([0.15,0.2,.95-0.15,.95-0.2])
	ax0.set_xlim(-1, 11)
	ax0.set_ylim(-2, 92)
	R0D = beta*M/(mu+gamma)
	R0S = R0D - (sigma[j]**2 * M**2)/(2*(mu+sigma[j]))
#
	ax0.plot(StoSIS.t, xdet,
		color='red',
		marker='.',
		alpha=1,
		lw=1,
		ms=1,
		label="Deterministic Solution"+'\n'
		+r'$R_0^D=$' + str(R0D) + ',\n'
		
	)
	ax0.plot(
		StoSIS.dt*StoSIS.tau,
		xst,
		color='green',
		marker='.',
		alpha=1,
		lw=1,
		ms=1,
		label=strSt + ',\n'
		+r'$R_0^S=$' + str(R0S) + ',\n'
		+r'$\sigma^2=$'+str(sigma[j]**2) + ',\n'
		+r'$\beta / N=$'+str(beta/M) +',\n'
	)
	
	if j >1:
		xi = (1.0/sigma[j]**2) * (np.sqrt(beta**2 - 2*(sigma[j]**2)*(mu+gamma)) -(beta - sigma[j]**2*M) )
		ax0.plot(StoSIS.dt*StoSIS.tau, xi * np.ones(xst.shape[0]),
			color='orange', 
			marker='.',
			alpha=1, 
			lw=0, 
			ms=1, 
			label=r'$\xi=$' + str(xi)
		)
	ax0.set_xlabel(r'$t$')
	ax0.set_ylabel(r'$I(t)$')
	ax0.grid()
	handles, labels = ax0.get_legend_handles_labels()
	ax0.legend(loc='best', fancybox=True, framealpha=0.5)
	#
	#Xzer0 = 1
	#
	StoSIS.SetParametersStochasticSIS(beta, mu, gamma, sigma[j], M, Xzero[1])
	h = StoSIS.Dt
	strEl = 'Euler h=' + str(h)
	strSt = 'SM h =' + str(h)
	xem = StoSIS.EM(flag=0)
	xst = StoSIS.SteklovMatus()
	xdet = StoSIS.DetSol()
	StoSIS.SaveData()
	#
	#
	ax1.set_xlim(-1, 11)
	ax1.set_ylim(-1, np.max([xem.max()+2,12]))
	ax1.plot(StoSIS.dt*StoSIS.tau,
		xst,
		color='green',
		marker='.',
		alpha=1,
		lw=0,
		ms=1,
		label=strSt
	)
	ax1.plot(StoSIS.t, xdet,
		color='red',
		marker='.',
		alpha=1.0,
		lw=1,
		ms=1,
		label="Deterministic Solution"  
	)
	#
	if j >1:
		xi = (1.0/sigma[j]**2) * (np.sqrt(beta**2 - 2*(sigma[j]**2)*(mu+gamma)) -(beta - sigma[j]**2*M) )
		ax1.plot(StoSIS.dt*StoSIS.tau, xi * np.ones(xst.shape[0]),
			color='orange', 
			marker='.',
			alpha=1, 
			lw=1, 
			ms=1
		)
	ax1.set_xlabel(r'$t$')
	ax1.set_ylabel(r'$I(t)$')
	ax1.grid()
	strFileName ='IStochastiSIS' + str(j)  + '.png'
	#
#
#	ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#	plt.tight_layout()
	plt.savefig(strFileName)
	plt.show()
	
'''
	Nmax = 10000
	DisXem=np.zeros(Nmax)
	DisXst=np.zeros(Nmax)
	
	DisXem[0] = xem[-1]
	DisXst[0] = xst[-1]
	
	for i in np.arange(Nmax-1):
		StoSIS.InitializeMesh(k, p, r, T0, T)
		if np.mod(k,10)== 0:
			print str(i+1)+'/1000'
			print '---------------------------------------------'
		DisXem[i] = StoSIS.EM(flag=0)[-1]
		DisXst[i] = StoSIS.SteklovMatus()[-1]
		print DisXem[i], DisXst[i]
	
	strFileName = 'DataDistribution' + str(j) + '.txt'
	np.savetxt(strFileName,
		np.transpose(
			(
				DisXem, DisXst
			)
		), fmt='%1.15f', delimiter='\t')
'''