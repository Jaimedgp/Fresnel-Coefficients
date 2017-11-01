import numpy as np
import matplotlib.pyplot as plt

################################################################################

class FresnelCoefficients:

	def __init__(self, n=1, ni=3.47-1.4j):

		self.n = n
		self.ni = ni
		self.phi = np.linspace(0,np.pi/2,200)

	def realCoeff(self):

		phii = [np.arcsin(np.sin(phi)*(self.n/self.ni)) for phi in self.phi ]

		rParl = []
		rOrt = []

		for i in range(len(phii)):
			rOrt.append(-np.sin(self.phi[i]-phii[i]) / 
				                                 np.sin(self.phi[i]+phii[i]))

			rParl.append(-np.tan(self.phi[i]-phii[i]) /
			                                     np.tan(self.phi[i]+phii[i]))

		plt.axhline(0, color='black').set_dashes([1,1,1,1])
		plt.plot(self.phi, rParl, label="$r_{||}$")
		plt.plot(self.phi, rOrt, label="$r_{\perp}$")
		
		plt.xticks([0,np.pi/2])
		plt.xlim([0,np.pi/2])
		
		plt.yticks([-1,0,1])
		plt.ylim([-1,1])

		plt.legend()
		plt.show()
