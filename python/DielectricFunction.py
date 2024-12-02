import numpy as np
import pandas as pd

c = 3E8
h = 6.626E-34
q = 1.602E-19

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


class LorentzDrude:
    def __init__(self, wp, res, wl):
        self.wp = wp # [eV]
        self.resonance = res # [, eV, eV]
        self.wl = wl  # [m]
        self.e = 1
        self.w_eV = (c*h/q) / self.wl # [eV]

    def DielectricFunctionCalc(self):

        for (fn, gn, wn) in self.resonance:
            self.e += np.power(self.wp, 2) * fn \
                /(np.power(wn, 2) - np.power(self.w_eV, 2) + 1j*self.w_eV*gn)

        self.n = np.sqrt((np.abs(self.e) + np.real(self.e)) / 2)
        self.k = np.sqrt((np.abs(self.e) - np.real(self.e)) / 2)

        df = pd.DataFrame(self.e)
        df.to_clipboard(excel=True)


class DrudeSommerfeld:
    def __init__(self, ec, wp, gamma, wl):
        self.e = ec # [eV]
        self.wp = wp # [eV]
        self.gamma = gamma # [eV]
        self.wl = wl # [m]

        self.w_eV = (c*h/q) / self.wl # [eV]

    def DielectricFunctionCalc(self):

        self.e = self.e - \
                 np.power(self.wp, 2) / (np.power(self.w_eV, 2) + np.power(self.gamma, 2)) - \
                 1j*self.gamma*np.power(self.wp, 2) / (self.w_eV * (np.power(self.w_eV, 2) + np.power(self.gamma, 2)))

        self.n = np.sqrt((np.abs(self.e) + np.real(self.e)) / 2)
        self.k = np.sqrt((np.abs(self.e) - np.real(self.e)) / 2)

        df = pd.DataFrame(self.e)
        df.to_clipboard(excel=True)


class ModifiedDrudeSommerfeld:
    def __init__(self, ec, wp, gamma, A, delta, wc, wl):

        self.e = ec # [eV]
        self.wp = wp # [eV]
        self.gamma = gamma # [eV]

        self.A = A
        self.delta = delta # [eV]
        self.wc = wc # [eV]

        self.wl = wl # [m]

        self.w_eV = (c*h/q) / self.wl # [eV]

    def DielectricFunctionCalc(self):

        self.e = self.e - \
                 np.power(self.wp, 2) / (np.power(self.w_eV, 2) + np.power(self.gamma, 2)) - \
                 1j*self.gamma*np.power(self.wp, 2) / (self.w_eV * (np.power(self.w_eV, 2) + np.power(self.gamma, 2))) - \
                 1j*self.A / (1 + np.exp(-(self.w_eV - self.wc)/self.delta))

        self.n = np.sqrt((np.abs(self.e) + np.real(self.e)) / 2)
        self.k = np.sqrt((np.abs(self.e) - np.real(self.e)) / 2)

        df = pd.DataFrame(self.e)
        df.to_clipboard(excel=True)
