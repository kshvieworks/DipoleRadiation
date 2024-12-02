import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import tkinter
import pathlib

e0 = 8.854E-12 # [F/m]
c = 3E8 # [m/s]
fd = pathlib.Path(__file__).parent.resolve()

# import numpy as np
# import DielectricFunction as DF
# wl = np.linspace(400E-9, 700E-9, 31)
# a = DF.ModifiedDrudeSommerfeld(9.84, 9.01, 0.072, 5.6, 0.17, 2.4, wl)
# a.DielectricFunctionCalc()

# import CoupledDipoleRadiation as CDR
# DDA = CDR.DiscreteDipoleApproximation(np.ones(wl.size), a.e, wl, 15E-9, 1E-9, 32E-9, 0)
# DDA.DrawFieldEnhancement('Both')

# import pandas as pd
# df = pd.DataFrame(np.abs(DDA.EffRadPattern[0, :]))
# df.to_clipboard(excel=True)

class DiscreteDipoleApproximation:

    def __init__(self, em, ep, wl, a, r, d, theta):

        self.e_m = em
        self.e_p = ep
        self.wl = wl
        self.a = a
        self.r = r
        self.d = d
        self.theta = theta
        self.phi = self.theta + np.deg2rad(0)

        self.n_m = np.sqrt((np.abs(self.e_m)+np.real(self.e_m))/2) + 1j*np.sqrt((np.abs(self.e_m)-np.real(self.e_m))/2)
        self.n_p = np.sqrt((np.abs(self.e_p)+np.real(self.e_p))/2) + 1j*np.sqrt((np.abs(self.e_p)-np.real(self.e_p))/2)
        self.wn0 = 2*np.pi/self.wl
        self.wn_m = self.n_m * self.wn0
        self.wn_p = self.n_p * self.wn0

        self.PolarTheta = np.tile(np.deg2rad(np.linspace(0,360,361)), (self.wl.size, 1)).T
        self.EffRadPattern = np.zeros((self.wl.size, int(self.PolarTheta.size/self.wl.size)), dtype=np.complex128)


    def DrawFieldEnhancement(self, rtype):

        alpha_eff = self.CalcEffectivePolarizability(self.d, self.phi, rtype)
        p_eff = alpha_eff * e0 * self.e_m

        fig = plt.figure(figsize=(12, 10))
        ax = fig.gca(projection='polar')
        plt.ion()
        plt.show()

        color = iter(plt.cm.rainbow(np.linspace(0, 1, self.wl.size)))
        ax.set_title("15nm Radius Coupled Au NanoParticle Radiation Pattern @ r = 1nm", fontsize=18)
        ax.grid(True)
        ax.set_theta_zero_location("N")
        ax.tick_params(axis='both', which='major', labelsize=20)

        self.EffRadPattern = self.CalcDipoleRad(self.r, self.PolarTheta, p_eff, rtype)

        ymax = 1.1 * np.round(int(np.power(np.max(np.abs(self.EffRadPattern)), 2)), decimals=-6)
        ax.set_ylim(0, 15E-13)

        for k in range(self.wl.size):

            c = next(color)
            asd = ax.plot(self.PolarTheta[:, k], np.power(np.abs(self.EffRadPattern[:, k]), 2), c=c, alpha = 0.7, label=f"{int(np.round((1E9*self.wl[k]), decimals=-1))}nm")
            ax.legend(loc=1, bbox_to_anchor=(1.2, 0.9), ncol=1, fancybox=True, shadow=True)
            plt.pause(0.001)
            # input("Press Enter to Continue.")

    def update(self, frame, ax, t, lines, x, y, wl):

        pre_line = lines[frame-1]
        pre_line.set_data(0, 0)

        line = lines[frame]
        line.set_data(x, np.power(np.abs(y[:, frame]), 2))
        t.set_text(f"{int(np.round((1E9*wl[frame]), decimals=-1))}nm")
        return pre_line, line, t

    def savegif(self):

        fig = plt.figure(figsize=(12, 10))
        ax = fig.gca(projection='polar')

        self.lines = []
        color = iter(plt.cm.rainbow(np.linspace(0, 1, self.wl.size)))

        for k in range(self.wl.size):
            c = next(color)
            lobj = ax.plot(0, 0, c=c, alpha = 0.7, lw=3, label=f"{int(np.round((1E9*self.wl[k]), decimals=-1))}nm")[0]
            self.lines.append(lobj)

        for line in self.lines:
            line.set_data(0, 0)

        ymax = 1.1 * np.round(int(np.power(np.max(np.abs(self.EffRadPattern)), 2)), decimals=-6)
        ax.set_ylim(0, 15E-13)
        # ax.set_ylim(0, ymax)
        ax.legend(loc=1, bbox_to_anchor=(1.2, 0.9), ncol=1, fancybox=True, shadow=True)
        ax.set_title("15nm Radius Coupled Au NanoParticle Radiation Pattern @ r = 1nm", fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.grid(True)
        ax.set_theta_zero_location("N")
        plt.show()

        t = ax.text(0, 0, "Ready", fontsize=20)

        ani = FuncAnimation(fig, self.update, frames=self.wl.size, fargs=(ax, t, self.lines, self.PolarTheta[:, 0], self.EffRadPattern, self.wl), interval=500, blit=False, repeat=True)

        filepath = tkinter.filedialog.asksaveasfilename(initialdir=f"{fd}/",
                                                        title="Save as",
                                                        filetypes=(("GIF Files", ".gif"),
                                                                   ("all files", "*")))
        filepath = f"{filepath}.gif"

        input("Press Enter to Continue.")

        ani.save(filepath, writer='imagemagick', fps=5)



    def CalcEffectivePolarizability(self, r, angle, rtype):

        return self.CalcPolarizability() / (1 - self.CalcDipoleRad(r, angle, self.CalcDipoleMoment(), rtype))

    def CalcPolarizability(self):

        return self.CalcDipoleMoment() / (e0 * self.e_m)

    def CalcDipoleMoment(self):

        return 4 * np.pi * e0 * self.e_m * np.power(self.a, 3) * (self.e_p-self.e_m) \
               / (self.e_p + 2*self.e_m)

    def CalcDipoleRadFar(self, r, angle, p):

        return -np.power(self.wn_m, 2) * p * np.power(np.e, 1j*self.wn_m*r) * np.power(np.sin(angle), 2) \
               / (4 * np.pi * e0 * self.e_m * r)

    def CalcDipoleRadNear(self, r, angle, p):

        return p * np.power(np.e, 1j*self.wn_m*r) * (1-1j*self.wn_m*r) * (3*np.power(np.cos(angle), 2) - 1) \
               / (4 * np.pi * e0 * self.e_m * np.power(r, 3))

    def CalcDipoleRad(self, r, angle, p, rtype):

        if rtype == 'Both':

            return self.CalcDipoleRadFar(r, angle, p) + self.CalcDipoleRadNear(r, angle, p)

        if rtype == 'Far':

            return self.CalcDipoleRadFar(r, angle, p)

        if rtype == 'Near':

            return self.CalcDipoleRadNear(r, angle, p)