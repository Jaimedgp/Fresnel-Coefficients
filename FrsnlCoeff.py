import numpy as np
from cmath import phase

################################################################################ 

class FresnelCoefficients:

    def __init__(self, n=1, ni=3.47-1.4j):

        self.n = n
        self.ni = ni
        self.phi = np.linspace(0,np.pi/2,200)
        self.phii = [np.arcsin(np.sin(phi)*(self.n/self.ni)) for phi in self.phi ]

    def reflexCoeff(self):

        rParl = [(self.n-self.ni) / (self.n + self.ni)]
        rOrt = [rParl[0]]

        for i in range(1,len(self.phii)):

            rOrtogonal = - (np.sin(self.phi[i]-self.phii[i]) /
                                               np.sin(self.phi[i]+self.phii[i])) 

            if not np.isnan(rOrtogonal): 
                rOrt.append(rOrtogonal)
            else:
                rOrt.append(1)

            rParalel = - (np.tan(self.phi[i]-self.phii[i]) / 
                                                np.tan(self.phi[i]+self.phii[i]))

            if not np.isnan(rParalel):
                rParl.append(rParalel)
            else:
                rParl.append(-1)

        return [rParl, rOrt]

    def transCoeff(self):

        tParl = [(2*self.n) / (self.n + self.ni)]
        tOrt = [(2*self.n) / (self.n + self.ni)]

        for i in range(1,len(self.phii)):

            tOrtogonal = ((2*self.n*np.cos(self.phi[i])) /
                     (self.n*np.cos(self.phi[i])+self.ni*np.cos(self.phii[i])))

            if not np.isnan(tOrtogonal): 
                tOrt.append(tOrtogonal)
            else:
                tOrt.append(0)

            tParalel = ((2*self.n*np.cos(self.phi[i])) /
                     (self.n*np.cos(self.phii[i])+self.ni*np.cos(self.phi[i])))

            if not np.isnan(tParalel):
                tParl.append(tParalel)
            else:
                tParl.append(0)

        return [tParl, tOrt]

    def reflexFact(self):

        rParl, rOrt = self.reflexCoeff()

        RParl = []
        ROrt = []

        for i in rParl:

            RParl.append(i**2)
        for n in rOrt:
            ROrt.append(n**2)

        return [RParl, ROrt]

    def transFact(self):

        RParl, ROrt = self.reflexFact()

        TParl = []
        TOrt = []

        for i in RParl:
            TParl.append(1-i)
        for n in ROrt:
            TOrt.append(1-n)

        return [TParl, TOrt]

    def imagCoeff(self):
    
        rParl, rOrt = self.reflexCoeff()

        rPhsParl, rPhsOrt, rModParl, rModOrt = [], [], [], []

        for i in range(len(rParl)):

            rPhsParl.append(phase(rParl[i]))
            rPhsOrt.append(phase(rOrt[i]))

            rModParl.append(abs(rParl[i]))
            rModOrt.append(abs(rOrt[i]))

        return [rPhsParl, rPhsOrt, rModParl, rModOrt]

    def plotImg(self):

        f, (ax1,ax2) = plt.subplots(1, 2)

        ax1.axhline(0, color='black').set_dashes([1,1,1,1])
        ax1.plot(self.phi, rModParl, label="$|r_{||}|$")
        ax1.plot(self.phi, rModOrt, label="$|r_{\perp}|$")

        ax1.set_xticks([0,np.pi/2])
        ax1.set_xlim([0,np.pi/2])

        ax1.set_yticks([min([min(rModParl), min(rModOrt)]),
                        max([max(rModParl), max(rModOrt)])])
        ax1.set_ylim([min([min(rModParl), min(rModOrt)]),
                      max([max(rModParl), max(rModOrt)])])

        ax1.legend()

        ax2.axhline(0, color='black').set_dashes([1,1,1,1])
        ax2.plot(self.phi, rPhsOrt, label="$\phi_{\perp}$")
        ax2.plot(self.phi, rPhsParl, label="$\phi_{||}$")

        ax2.set_xticks([0,np.pi/2])
        ax2.set_xlim([0,np.pi/2])

        ax2.set_yticks([min([min(rPhsParl), min(rPhsOrt)]),
                        max([max(rPhsParl), max(rPhsOrt)])])
        ax2.set_ylim([min([min(rPhsParl), min(rPhsOrt)]),
                      max([max(rPhsParl), max(rPhsOrt)])])

        ax2.legend()

        plt.show()