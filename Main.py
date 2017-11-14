from FrsnlCoeff import FresnelCoefficients as FC
import matplotlib.pyplot as plt
import numpy as np

class Main:

    def __init__(self, n=1, ni=3.47-1.4j):

        self.n = n
        self.ni = ni

        if complex == type(n):
            self.plotImg()
        elif complex == type(ni):
            self.plotImg()
        else:
            self.plotreal()

    def plotreal(self):

        case1 = FC(self.n, self.ni)

        rParl, rOrt = case1.reflexCoeff()
        tParl, tOrt = case1.transCoeff()
        RParl, ROrt = case1.reflexFact()
        TParl, TOrt = case1.transFact()


        f, plots = plt.subplots(2, 2)

        plots[0,0].axhline(0, color='black').set_dashes([1,1,1,1])
        plots[0,0].plot(case1.phi, rParl, label="$r_{||}$")
        plots[0,0].plot(case1.phi, rOrt, label="$r_{\perp}$")

        plots[0,0].set_xticks([0,np.pi/2])
        plots[0,0].set_xlim([0,np.pi/2])

        plots[0,0].set_xlabel("$\\varphi/rad$")
        plots[0,0].set_ylabel("$r$")

        plots[0,0].set_yticks([-1,0,1])
        plots[0,0].set_ylim([-1,1])

        plots[0,0].legend()

        plots[0,1].plot(case1.phi, tParl, label="$t_{||}$")
        plots[0,1].plot(case1.phi, tOrt, label="$t_{\perp}$")

        plots[0,1].set_xticks([0,np.pi/2])
        plots[0,1].set_xlim([0,np.pi/2])

        plots[0,1].set_xlabel("$\\varphi/rad$")
        plots[0,1].set_ylabel("$t$")

        plots[0,1].set_yticks([min([min(tParl), min(tOrt)]),
                               max([max(tParl), max(tOrt)])])
        plots[0,1].set_ylim([min([min(tParl), min(tOrt)]),
                               max([max(tParl), max(tOrt)])])

        plots[0,1].legend()

        plots[1,0].plot(case1.phi, RParl, label="$R_{||}$")
        plots[1,0].plot(case1.phi, ROrt, label="$R_{\perp}$")

        plots[1,0].set_xticks([0,np.pi/2])
        plots[1,0].set_xlim([0,np.pi/2])

        plots[1,0].set_xlabel("$\\varphi/rad$")
        plots[1,0].set_ylabel("$R$")

        plots[1,0].set_yticks([0,1])
        plots[1,0].set_ylim([0,1])

        plots[1,0].legend()

        plots[1,1].plot(case1.phi, TParl, label="$T_{||}$")
        plots[1,1].plot(case1.phi, TOrt, label="$T_{\perp}$")

        plots[1,1].set_xticks([0,np.pi/2])
        plots[1,1].set_xlim([0,np.pi/2])

        plots[1,1].set_xlabel("$\\varphi/rad$")
        plots[1,1].set_ylabel("$T$")

        plots[1,1].set_yticks([0,1])
        plots[1,1].set_ylim([0,1])

        plots[1,1].legend()

        plt.show()

    def plotImg(self):

        case1 = FC(self.n, self.ni)

        rPhsParl, rPhsOrt, rModParl, rModOrt = case1.imagCoeff()

        f, (ax1,ax2) = plt.subplots(1, 2)

        ax1.plot(case1.phi, rModParl, label="$|r_{||}|$")
        ax1.plot(case1.phi, rModOrt, label="$|r_{\perp}|$")

        ax1.set_xlabel("$\\varphi/rad$")
        ax1.set_ylabel("$|r|$")

        ax1.set_xticks([0,np.pi/2])
        ax1.set_xlim([0,np.pi/2])

        ax1.set_yticks([0,
                        max([max(rModParl), max(rModOrt)])])
        ax1.set_ylim([0,
                      max([max(rModParl), max(rModOrt)])])

        ax1.legend()

        ax2.plot(case1.phi, rPhsParl, label="$\phi_{||}$")
        ax2.plot(case1.phi, rPhsOrt, label="$\phi_{\perp}$")

        ax2.set_xticks([0,np.pi/2])
        ax2.set_xlim([0,np.pi/2])

        ax2.set_xlabel("$\\varphi/rad$")
        ax2.set_ylabel("$\phi$")

        ax2.set_yticks([min([min(rPhsParl), min(rPhsOrt)]),
                        max([max(rPhsParl), max(rPhsOrt)])])
        ax2.set_ylim([min([min(rPhsParl), min(rPhsOrt)]),
                      max([max(rPhsParl), max(rPhsOrt)])])

        ax2.legend()

        plt.show()

    def plotMoreno(self):

        a, b = self.ni.real, self.ni.imag
        k = int(abs(b*10))
        incr = b/abs(b)

        for n in range(0,k):
            s = complex(a, (b-incr*(n*0.1)) )

            casen = FC(self.n, s)
            rPhsParl, rPhsOrt, rModParl, rModOrt = casen.imagCoeff()

            plt.plot(casen.phi, rPhsParl, label="$\phi_{||} n'= $"+ str(s))

        s = complex(a, 0 )

        casen = FC(self.n, s)
        rPhsParl, rPhsOrt, rModParl, rModOrt = casen.imagCoeff()

        rPhsP = [np.pi]

        for r in range(1,len(rPhsParl)):
        	rPhsP.append(-1*rPhsParl[r])

        plt.plot(casen.phi, rPhsP, label="$\phi_{||} n'= $"+ str(s))

        plt.xlabel("$\\varphi/rad$")
        plt.ylabel("$\phi$")
        plt.legend()
        plt.show()