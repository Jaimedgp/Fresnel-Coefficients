import numpy as np
import matplotlib.pyplot as plt

################################################################################

class FresnelCoefficients:

	def __init__(self, n=1, ni=3.47-1.4j):

		self.n = n
		self.ni = ni
		self.phi = np.linspace(0,np.pi/2,200)

		if 'complex' in [type(n), type(ni)]:
			print 'gilipollas'
		else:
			if ni > n:
				phii = [np.arcsin(np.sin(phi)*(self.n/self.ni)) for phi in self.phi ]
				self.realCoeff_nipn(phii)
			elif n > ni:
				self.realCoeff_npni()

		self.plotCoeff()

	def realCoeff_nipn(self, phii):

		self.rParl = []
		self.rOrt = []
		self.tParl = []
		self.tOrt = []

		self.RParl = []
		self.ROrt = []
		self.TParl = []
		self.TOrt = []

		for i in range(len(phii)):
			self.rOrt.append(-np.sin(self.phi[i]-phii[i]) / 
				                                 np.sin(self.phi[i]+phii[i]))

			self.rParl.append(-np.tan(self.phi[i]-phii[i]) /
			                                     np.tan(self.phi[i]+phii[i]))

			self.tOrt.append((2*np.sin(phii[i])*np.cos(self.phi[i])) / 
				                                 np.sin(self.phi[i]+phii[i]))

			self.tParl.append((2*np.sin(phii[i])*np.cos(self.phi[i])) / 
				     np.sin(self.phi[i]+phii[i])*np.cos(self.phi[i]-phii[i]))

			self.ROrt.append(self.rOrt[i]**2)
			self.RParl.append(self.rParl[i]**2)
			self.TParl.append(1-self.RParl[i])
			self.TOrt.append(1-self.ROrt[i])

	def realCoeff_npni(self):

		phiL = np.arcsin(self.ni/self.n)
		print phiL

		phiR = np.linspace(0, phiL)
		phiiR = [np.arcsin(np.sin(phi)*(self.n/self.ni)) for phi in phiR ]

		self.realCoeff_nipn(phiiR)

		while len(self.rParl) != len(self.phi):
			self.rParl.append(-1)
		while len(self.tParl) != len(self.phi):
			self.tParl.append(0)
		while len(self.rOrt) != len(self.phi):
			self.rOrt.append(1)
		while len(self.tOrt) != len(self.phi):
			self.tOrt.append(0)

	def plotCoeff(self):

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
		
		plots[0,1].set_yticks([0,1])
		plots[0,1].set_ylim([0,1])

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