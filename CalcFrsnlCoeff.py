import numpy as np
import matplotlib.pyplot as plt

################################################################################ 

class FresnelCoefficients:

	def __init__(self, n=1, ni=3.47-1.4j):

		self.n = n
		self.ni = ni
		self.phi = np.linspace(0,np.pi/2,200)

		if complex == type(n):
			self.imaginCoeff(n.real, n.imag)
			self.plotImgCoeff()
		elif complex == type(ni):
			self.imaginCoeff(ni.real, ni.imag)
			self.plotImgCoeff()
		else:
			phii = [np.arcsin(np.sin(phi)*(self.n/self.ni)) for phi in self.phi ]
			self.realCoeff(phii)
			self.plotrealCoeff()

	def imaginCoeff(self, a, b):

		self.rModParl = []
		self.rModOrt = []
		self.rphsParl = []
		self.rphsOrt = []

		for i in range(len(self.phi)):

			ncos = self.ni*np.sqrt(1-(((self.n*np.sin(self.phi[i]))**2)/(self.ni)**2))
			cos = np.sqrt(1-(((self.n*np.sin(self.phi[i]))**2)/(self.ni)**2))

			c = ncos.real
			d = ncos.imag

			A = self.n*cos - a*np.cos(self.phi[i])
			G = self.n*cos + a*np.cos(self.phi[i])
			B = b*np.cos(self.phi[i])

			ReParl = (A*G - B**2) / (G**2 + B**2)
			ImgParl = (-B*G - A*B) / (G**2 + B**2)

			self.rModParl.append(np.sqrt(ReParl**2 + ImgParl**2))
			self.rphsParl.append(2*np.arctan(ImgParl/ReParl))

			X = self.n*np.cos(self.phi[i]) - c
			Y = self.n*np.cos(self.phi[i]) + c

			ReOrt = (X*Y - d**2) / (Y**2 + d**2)
			ImgOrt = -(X+Y)*d / (Y**2 + d**2)

			#A = self.n*np.cos(self.phi[i]) - a*cos
			#G = self.n*np.cos(self.phi[i]) + a*cos
			#B = b*cos

			#ReOrt = (A*G - B**2) / (G**2 + B**2)
			#ImgOrt = -(B*G + A*B) / (G**2 + B**2)

			self.rModOrt.append(np.sqrt(ReOrt**2 + ImgOrt**2))
			self.rphsOrt.append(2*np.arctan(ImgOrt/ReOrt))

	def plotImgCoeff(self):

		f, (ax1,ax2) = plt.subplots(1, 2)

		ax1.axhline(0, color='black').set_dashes([1,1,1,1])
		ax1.plot(self.phi, self.rModParl, label="$|r_{||}|$")
		ax1.plot(self.phi, self.rModOrt, label="$|r_{\perp}|$")
		
		ax1.set_xticks([0,np.pi/2])
		ax1.set_xlim([0,np.pi/2])
		
		ax1.set_yticks([min([min(self.rModParl), min(self.rModOrt)]),
			            max([max(self.rModParl), max(self.rModOrt)])])
		ax1.set_ylim([min([min(self.rModParl), min(self.rModOrt)]),
			          max([max(self.rModParl), max(self.rModOrt)])])

		ax1.legend()

		ax2.axhline(0, color='black').set_dashes([1,1,1,1])
		ax2.plot(self.phi, self.rphsParl, label="$|\phi_{||}|$")
		ax2.plot(self.phi, self.rphsOrt, label="$|\phi_{\perp}|$")
		
		ax2.set_xticks([0,np.pi/2])
		ax2.set_xlim([0,np.pi/2])
		
		ax2.set_yticks([min([min(self.rphsParl), min(self.rphsOrt)]),
			            max([max(self.rphsParl), max(self.rphsOrt)])])
		ax2.set_ylim([min([min(self.rphsParl), min(self.rphsOrt)]),
			          max([max(self.rphsParl), max(self.rphsOrt)])])

		ax2.legend()

		plt.show()

	def realCoeff(self, phii):

		self.rParl = [(self.n-self.ni) / (self.n+self.ni)]
		self.rOrt = [self.rParl[0]]

		self.RParl = [self.rParl[0]**2]
		self.ROrt = [self.rParl[0]**2]
		self.TParl = [1-self.RParl[0]]
		self.TOrt = [self.TParl[0]]
		self.tParl = [ np.sqrt((self.n*np.cos(self.phi[0])*self.TParl[0] / (self.ni*np.cos(phii[0]))))]
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

			tOrt = ((2*np.sin(phii[i])*np.cos(self.phi[i])) / 
				                                 np.sin(self.phi[i]+phii[i]))

			if not np.isnan(tOrt): 
				self.tOrt.append(tOrt)
			else:
				self.tOrt.append(0)

			tParl = ((2*np.sin(phii[i])*np.cos(self.phi[i])) / 
				     np.sin(self.phi[i]+phii[i])*np.cos(self.phi[i]-phii[i]))

			if not np.isnan(tParl):
				self.tParl.append(tParl)
			else:
				self.tParl.append(0)

			self.ROrt.append(self.rOrt[i]**2)
			self.RParl.append(self.rParl[i]**2)
			self.TParl.append(1-self.RParl[i])
			self.TOrt.append(1-self.ROrt[i])

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