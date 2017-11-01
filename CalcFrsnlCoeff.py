import numpy as np
from math import pi
import matplotlib.pyplot


class FresnelCoefficients:

	def __init__(self, n=1, ni=3.47-1.4j):

		self.n = n
		self.ni = ni
		self.phi = np.linspace(0,pi/2,200)

	def realCoeff(self):

		