import numpy as np
import matplotlib.pyplot as plt
from cmath import phase

################################################################################ 

class FresnelCoefficients:

	def __init__(self, n=1, ni=3.47-1.4j):

		self.n = n
		self.ni = ni
		self.phi = np.linspace(0,np.pi/2,200)
		self.phii = [np.arcsin(np.sin(phi)*(self.n/self.ni)) for phi in self.phi ]

		if type(ni) == complex:
			self.realCoeff(self.phii)
			self.ImagCoeff()
			self.plotImgCoeff()
		elif type(n) == complex:
			self.realCoeff(self.phii)
			self.ImagCoeff()
			self.plotImgCoeff()
		else:
			self.realCoeff(self.phii)
			self.plotrealCoeff()


	def ImagCoeff(self):

		self.rPhsParl = []
		self.rPhsOrt = []
		self.rModParl = []
		self.rModOrt = []

		for i in range(len(self.rParl)):

			self.rPhsParl.append(phase(self.rParl[i]))
			self.rPhsOrt.append(phase(self.rOrt[i]))

			self.rModParl.append(abs(self.rParl[i]))
			self.rModOrt.append(abs(self.rOrt[i]))


	def realCoeff(self, phii):

		self.rParl = [(self.n-self.ni) / (self.n + self.ni)]
		self.rOrt = [self.rParl[0]]
		self.RParl = [self.rParl[0]**2]
		self.ROrt = [self.rParl[0]**2]
		
		self.TParl = [1-self.RParl[0]]
		self.TOrt = [1-self.ROrt[0]]

		self.tParl = [np.sqrt(self.TParl[0]*(self.n/self.ni))]
		self.tOrt = [self.tParl[0]]


		for i in range(1,len(phii)):

			rOrt = -np.sin(self.phi[i]-phii[i]) / np.sin(self.phi[i]+phii[i]) 
			
			if not np.isnan(rOrt): 
				self.rOrt.append(rOrt)
			else:
				self.rOrt.append(1)

			rParl = -np.tan(self.phi[i]-phii[i]) / np.tan(self.phi[i]+phii[i])

			if not np.isnan(rParl):
				self.rParl.append(rParl)
			else:
				self.rParl.append(-1)

			tOrt = ((2*self.n*np.cos(self.phi[i])) /
						                            (self.n*np.cos(self.phi[i])+self.ni*np.cos(phii[i])))

			if not np.isnan(tOrt): 
				self.tOrt.append(tOrt)
			else:
				self.tOrt.append(0)

			tParl = ((2*self.n*np.cos(self.phi[i])) /
						                            (self.n*np.cos(phii[i])+self.ni*np.cos(self.phi[i])))

			if not np.isnan(tParl):
				self.tParl.append(tParl)
			else:
				self.tParl.append(0)

			self.ROrt.append(self.rOrt[i]**2)
			self.RParl.append(self.rParl[i]**2)
			self.TParl.append(1-self.RParl[i])
			self.TOrt.append(1-self.ROrt[i])

		TParl = [ (self.n/self.ni)*(np.cos(self.phi[i])/np.cos(phii[i]) * self.tParl[i]**2) for i in range(len(self.tParl))]
		TOrt = [ (self.n/self.ni)*(np.cos(self.phi[i])/np.cos(phii[i]) * self.tOrt[i]**2) for i in range(len(self.tOrt))]

	def plotrealCoeff(self):

		f, plots = plt.subplots(2, 2)

		plots[0,0].axhline(0, color='black').set_dashes([1,1,1,1])
		plots[0,0].plot(self.phi, self.rParl, label="$r_{||}$")
		plots[0,0].plot(self.phi, self.rOrt, label="$r_{\perp}$")
		
		plots[0,0].set_xticks([0,np.pi/2])
		plots[0,0].set_xlim([0,np.pi/2])
		
		plots[0,0].set_yticks([-1,0,1])
		plots[0,0].set_ylim([-1,1])

		plots[0,0].legend()

		plots[0,1].plot(self.phi, self.tParl, label="$t_{||}$")
		plots[0,1].plot(self.phi, self.tOrt, label="$t_{\perp}$")
		
		plots[0,1].set_xticks([0,np.pi/2])
		plots[0,1].set_xlim([0,np.pi/2])
		
		plots[0,1].set_yticks([min([min(self.tParl), min(self.tOrt)]),
			                   max([max(self.tParl), max(self.tOrt)])])
		plots[0,1].set_ylim([min([min(self.tParl), min(self.tOrt)]),
			                   max([max(self.tParl), max(self.tOrt)])])

		plots[0,1].legend()

		plots[1,0].plot(self.phi, self.RParl, label="$R_{||}$")
		plots[1,0].plot(self.phi, self.ROrt, label="$R_{\perp}$")
		
		plots[1,0].set_xticks([0,np.pi/2])
		plots[1,0].set_xlim([0,np.pi/2])
		
		plots[1,0].set_yticks([0,1])
		plots[1,0].set_ylim([0,1])

		plots[1,0].legend()

		plots[1,1].plot(self.phi, self.TParl, label="$T_{||}$")
		plots[1,1].plot(self.phi, self.TOrt, label="$T_{\perp}$")
		
		plots[1,1].set_xticks([0,np.pi/2])
		plots[1,1].set_xlim([0,np.pi/2])
		
		plots[1,1].set_yticks([0,1])
		plots[1,1].set_ylim([0,1])

		plots[1,1].legend()

		plt.show()

	def plotImgCoeff(self):

		f, (ax1,ax2) = plt.subplots(1, 2)

		ax1.plot(self.phi, self.rModParl, label="$|r_{||}|$")
		ax1.plot(self.phi, self.rModOrt, label="$|r_{\perp}|$")
		
		ax1.set_xticks([0,np.pi/2])
		ax1.set_xlim([0,np.pi/2])
		
		ax1.set_yticks([min([min(self.rModParl), min(self.rModOrt)]),
			            max([max(self.rModParl), max(self.rModOrt)])])
		ax1.set_ylim([min([min(self.rModParl), min(self.rModOrt)]),
			          max([max(self.rModParl), max(self.rModOrt)])])

		ax1.legend()

		ax2.plot(self.phi, self.rPhsParl, label="$\phi_{||}$")
		ax2.plot(self.phi, self.rPhsOrt, label="$\phi_{\perp}$")
		
		ax2.set_xticks([0,np.pi/2])
		ax2.set_xlim([0,np.pi/2])
		
		ax2.set_yticks([min([min(self.rPhsParl), min(self.rPhsOrt)]),
			            max([max(self.rPhsParl), max(self.rPhsOrt)])])
		ax2.set_ylim([min([min(self.rPhsParl), min(self.rPhsOrt)]),
			          max([max(self.rPhsParl), max(self.rPhsOrt)])])

		ax2.legend()

		plt.show()