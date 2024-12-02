import numpy as np
import matplotlib.pyplot as plt
from scipy import special as sps

'''''''''''''''''''''''''''''''''Lorentz-Drude Model Implementation for Au'''''''''''''''''''''''''''''''''''''''''''''
# import numpy as np
# import DielectricFunction as DF
# wl = np.linspace(400E-9, 1000E-9, 601)
# a = DF.LorentzDrude(9.030, [(0.760,0.053,0), (0.024,0.241,0.415), (0.010,0.345,0.830), (0.071,0.870,2.969), (0.601,2.494,4.304), (4.384,2.214,13.320)], wl)
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

''''''''''''''''''''''''''''''''Drude-Sommerfeld Model Implementation for Au'''''''''''''''''''''''''''''''''''''''
# import numpy as np
# import DielectricFunction as DF
# wl = np.linspace(400E-9, 1000E-9, 601)
# a = DF.DrudeSommerfeld(9.84, 9.01, 0.072, wl)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''''''''''''''''''''''''''''Modified Drude-Sommerfeld Model Implementation for Au'''''''''''''''''''''''''''''''''''
# import numpy as np
# import DielectricFunction as DF
# wl = np.linspace(400E-9, 1000E-9, 601)
# a = DF.ModifiedDrudeSommerfeld(9.84, 9.01, 0.072, 5.6, 0.17, 2.4, wl)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''''''''''''''''''''''''''''Constant Refractive Index Material'''''''''''''''''''''''''''''''''''
# import numpy as np
# wl = np.linspace(400E-9, 1000E-9, 601)
# e = 2.2*np.ones(wl.size)
# em = 1*np.ones(wl.size)
# import MieScattering as MS
# M = MS.Mie(em, e, wl, 3.5E-6)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''''''''''''''''''''''''''''Refractive Index from Clipboard'''''''''''''''''''''''''''''''''''
# import pandas as pd
# import numpy as np
# d = np.array(pd.read_clipboard())
# wl = d[:, 0]
# n = d[:, 1]
# k = d[:, 2]
# e = n**2 - k**2 -1j*2*n*k
# em = 1.33*np.ones(wl.size)
# import MieScattering as MS
# M = MS.Mie(em, e, wl, 3.5E-9)
# M.MieParam(10)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



# a.DielectricFunctionCalc()
# em = 1.77*np.ones(601)
# import MieScattering as MS
# M = MS.Mie(em, a.e, a.wl, 50E-9)
# M.MieParam(10)

# import pandas as pd
# df = pd.DataFrame(M.Cs)
# df.to_clipboard(excel=True)

class Mie:
    def __init__(self, e_m, e_p, wl, a):

        self.mu_m = 1 # Medium
        self.mu_p = 1 # Particle

        self.e_m = e_m # Medium
        self.e_p = e_p # Particle

        self.n_m = np.sqrt((np.abs(self.e_m)+np.real(self.e_m))/2) + 1j*np.sqrt((np.abs(self.e_m)-np.real(self.e_m))/2)
        self.n_p = np.sqrt((np.abs(self.e_p)+np.real(self.e_p))/2) + 1j*np.sqrt((np.abs(self.e_p)-np.real(self.e_p))/2)

        self.wl = wl
        self.k0 = 2*np.pi/self.wl
        self.k_m = self.n_m * self.k0
        self.k_p = self.n_p * self.k0

        self.a = a
        self.rho_m = self.a * self.k_m
        self.rho_p = self.a * self.k_p

        self.an = np.zeros(self.wl.size)
        self.bn = np.zeros(self.wl.size)
        self.cn = np.zeros(self.wl.size)
        self.dn = np.zeros(self.wl.size)

        self.Cs = np.zeros(self.wl.size)
        self.Ce = np.zeros(self.wl.size)
        self.Ca = np.zeros(self.wl.size)

    def CalcCs(self, a, b, k):

        return 2*np.pi*(2*k + 1)*(np.power(np.abs(a), 2) + np.power(np.abs(b), 2)) / (np.power(self.k_m, 2))

    def CalcCe(self, a, b, k):

        return 2*np.pi*(2*k + 1)*np.real(a + b) / (np.power(self.k_m, 2))

    def CalcCa(self, a, b, k):

        return (self.CalcCe(a, b, k) - self.CalcCs(a, b, k))

    def MieParam(self, n):

        fig = plt.figure()
        ax = fig.subplots()
        cx, = ax.plot([], [], 'r')
        color = iter(plt.cm.rainbow(np.linspace(0, 1, n)))
        plt.show(block=False)

        for ll in range(n):
            k = ll + 1

            ak = (self.mu_m * np.power(self.n_p, 2) * self.SphBes(k, self.rho_p, False, 'j') * (self.rho_m * self.SphBes(k, self.rho_m, True, 'j') + self.SphBes(k, self.rho_m, False, 'j')) - \
                  self.mu_p * np.power(self.n_m, 2) * self.SphBes(k, self.rho_m, False, 'j') * (self.rho_p * self.SphBes(k, self.rho_p, True, 'j') + self.SphBes(k, self.rho_p, False, 'j'))) \
               / (self.mu_m * np.power(self.n_p, 2) * self.SphBes(k, self.rho_p, False, 'j') * (self.rho_m * self.SphBes(k, self.rho_m, True, 'h') + self.SphBes(k, self.rho_m, False, 'h')) - \
                  self.mu_p * np.power(self.n_m, 2) * self.SphBes(k, self.rho_m, False, 'h') * (self.rho_p * self.SphBes(k, self.rho_p, True, 'j') + self.SphBes(k, self.rho_p, False, 'j')))

            bk = (self.mu_p * self.SphBes(k, self.rho_p, False, 'j') * (self.rho_m * self.SphBes(k, self.rho_m, True, 'j') + self.SphBes(k, self.rho_m, False, 'j')) - \
                  self.mu_m * self.SphBes(k, self.rho_m, False, 'j') * (self.rho_p * self.SphBes(k, self.rho_p, True, 'j') + self.SphBes(k, self.rho_p, False, 'j'))) \
               / (self.mu_p * self.SphBes(k, self.rho_p, False, 'j') * (self.rho_m * self.SphBes(k, self.rho_m, True, 'h') + self.SphBes(k, self.rho_m, False, 'h')) - \
                  self.mu_m * self.SphBes(k, self.rho_m, False, 'h') * (self.rho_p * self.SphBes(k, self.rho_p, True, 'j') + self.SphBes(k, self.rho_p, False, 'j')))

            self.an = self.an + ak
            self.bn = self.bn + bk

            self.Cs = self.Cs + self.CalcCs(ak, bk, k)
            self.Ce = self.Ce + self.CalcCe(ak, bk, k)
            self.Ca = self.Ca + self.CalcCa(ak, bk, k)

            c = next(color)
            cx.set_xdata(self.wl)
            cx.set_ydata(np.abs(self.Ce))
            ax.plot(self.wl, np.abs(self.CalcCe(ak, bk, k)), c=c)
            plt.draw()
            plt.pause(0.1)

            asdf = 1


    def SphBes(self, nu, x, d, ftype):

        if ftype == 'j':
            return sps.spherical_jn(nu, x, d)

        if ftype == 'y':
            return sps.spherical_yn(nu, x, d)

        if ftype == 'h':
            return (sps.spherical_jn(nu, x, d) + 1j*sps.spherical_yn(nu, x, d))

    # def MieParam(self, n):

    #     for ll in range(n):
    #         k = ll + 1
    #
    #         ak = (self.mu_m * np.power(self.n_p, 2) * self.SphBes(k, self.rho_p, False, 'j') * (self.rho_m * self.SphBes(k, self.rho_m, True, 'j') + self.SphBes(k, self.rho_m, False, 'j')) - \
    #              self.mu_p * np.power(self.n_m, 2) * self.SphBes(k, self.rho_m, False, 'j') * (self.rho_p * self.SphBes(k, self.rho_p, True, 'j') + self.SphBes(k, self.rho_p, False, 'j'))) \
    #            / (self.mu_m * np.power(self.n_p, 2) * self.SphBes(k, self.rho_p, False, 'j') * (self.rho_m * self.SphBes(k, self.rho_m, True, 'h') + self.SphBes(k, self.rho_m, False, 'h')) - \
    #              self.mu_p * np.power(self.n_m, 2) * self.SphBes(k, self.rho_m, False, 'h') * (self.rho_p * self.SphBes(k, self.rho_p, True, 'j') + self.SphBes(k, self.rho_p, False, 'j')))
    #
    #         bk = (self.mu_p * self.SphBes(k, self.rho_p, False, 'j') * (self.rho_m * self.SphBes(k, self.rho_m, True, 'j') + self.SphBes(k, self.rho_m, False, 'j')) - \
    #              self.mu_m * self.SphBes(k, self.rho_m, False, 'j') * (self.rho_p * self.SphBes(k, self.rho_p, True, 'j') + self.SphBes(k, self.rho_p, False, 'j'))) \
    #            / (self.mu_p * self.SphBes(k, self.rho_p, False, 'j') * (self.rho_m * self.SphBes(k, self.rho_m, True, 'h') + self.SphBes(k, self.rho_m, False, 'h')) - \
    #              self.mu_m * self.SphBes(k, self.rho_m, False, 'h') * (self.rho_p * self.SphBes(k, self.rho_p, True, 'j') + self.SphBes(k, self.rho_p, False, 'j')))
    #
    #         self.an = self.an + ak
    #         self.bn = self.bn + bk
    #
    #         self.Cs = self.Cs + self.CalcCs(ak, bk, k)
    #         self.Ce = self.Ce + self.CalcCe(ak, bk, k)
    #         self.Ca = self.Ca + self.CalcCa()
    #
    # def SphBes(self, nu, x, d, ftype):
    #
    #     if d:
    #         return self.dSphBes(nu, x, ftype)
    #
    #     if ftype == 'j':
    #         return np.sqrt(np.pi/(2*x))*sps.jv(nu+0.5, x)
    #
    #     if ftype == 'y':
    #         return np.sqrt(np.pi/(2*x))*sps.yv(nu+0.5, x)
    #
    #     if ftype == 'h':
    #         return np.sqrt(np.pi/(2*x))*sps.hankel1(nu+0.5, x)
    #
    # def dSphBes(self, nu, x, ftype):
    #
    #     return (nu * self.SphBes(nu - 1, x, False, ftype) - (nu + 1) * self.SphBes(nu + 1, x, False, ftype)) / (2 * nu + 1)
