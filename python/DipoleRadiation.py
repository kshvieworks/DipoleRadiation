import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import tkinter
import pathlib
e0 = 8.854E-12 # [F/m]
c = 3E8 # [m/s]

# import numpy as np
# import DielectricFunction as DF
# wl = np.linspace(400E-9, 700E-9, 31)
# a = DF.ModifiedDrudeSommerfeld(9.84, 9.01, 0.072, 5.6, 0.17, 2.4, wl)
# a.DielectricFunctionCalc()

# import DipoleRadiation as DR
# e = DR.DipoleRad(np.ones(wl.size), a.e, wl, 15E-9)
# e.DrawPolarPlot(1E-9, 'Both')

# import pandas as pd
# df = pd.DataFrame(np.abs(e.DipoleRadPTN[:, 0]))
# df.to_clipboard(excel=True)

fd = pathlib.Path(__file__).parent.resolve()

class DipoleRad:

    def __init__(self, em, ep, wl, a):

        self.e_m = em
        self.e_p = ep
        self.wl = wl
        self.a = a

        self.n_m = np.sqrt((np.abs(self.e_m)+np.real(self.e_m))/2) + 1j*np.sqrt((np.abs(self.e_m)-np.real(self.e_m))/2)
        self.n_p = np.sqrt((np.abs(self.e_p)+np.real(self.e_p))/2) + 1j*np.sqrt((np.abs(self.e_p)-np.real(self.e_p))/2)

        self.wn0 = 2*np.pi/self.wl
        self.wn_m = self.n_m * self.wn0
        self.wn_p = self.n_p * self.wn0

        self.theta = np.deg2rad(np.linspace(0,360,361))

        self.DipoleRadPTN = np.zeros((self.wl.size, self.theta.size), dtype=np.complex128)

    def CalcDipoleMoment(self, k):

        return 4 * np.pi * e0 * self.e_m[k] * np.power(self.a, 3) * (self.e_p[k]-self.e_m[k]) \
               / (self.e_p[k] + 2*self.e_m[k])

    def CalcDipoleRadFar(self, r, k):

        p = self.CalcDipoleMoment(k)

        return -np.power(self.wn_m[k], 2) * p * np.power(np.e, 1j*self.wn_m[k]*r) * np.sin(self.theta) \
               / (4 * np.pi * e0 * self.e_m[k] * r)

    def CalcDipoleRadNearRcomp(self, r, k):

        p = self.CalcDipoleMoment(k)

        return p * np.power(np.e, 1j*self.wn_m[k]*r) * (1-1j*self.wn_m[k]*r) * (2*np.cos(self.theta)) \
               / (4 * np.pi * e0 * self.e_m[k] * np.power(r, 3))

    def CalcDipoleRadNearTcomp(self, r, k):

        p = self.CalcDipoleMoment(k)

        return p * np.power(np.e, 1j*self.wn_m[k]*r) * (1-1j*self.wn_m[k]*r) * (-1*np.sin(self.theta)) \
               / (4 * np.pi * e0 * self.e_m[k] * np.power(r, 3))

    def CalcDipoleRad(self, r, k, type):

        if type == 'Both':

            return np.power((self.CalcDipoleRadFar(r, k) + self.CalcDipoleRadNearTcomp(r, k)), 2) + np.power(self.CalcDipoleRadNearRcomp(r, k), 2)

        if type == 'Far':

            return self.CalcDipoleRadFar(r, k)

        if type == 'Near':

            return self.CalcDipoleRadNear(r, k)

    def DrawPolarPlot(self, r, type):

        fig = plt.figure(figsize=(12, 10))
        ax = fig.gca(projection='polar')
        plt.ion()
        plt.show()

        color = iter(plt.cm.rainbow(np.linspace(0, 1, self.wl.size)))
        ax.set_title("15nm Radius Au NanoParticle Radiation Pattern", fontsize=20)
        ax.grid(True)
        ymax = 1.1 * np.round(int(np.max(np.abs(self.DipoleRadPTN[:, :]))), decimals=-6)
        # ax.set_ylim(0, 15E-13)
        ax.set_ylim(0, ymax)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.set_theta_zero_location("E")


        for k in range(self.wl.size):
            self.DipoleRadPTN[k, :] = self.CalcDipoleRad(r, k, type)
            c = next(color)
            asd = ax.plot(self.theta, np.abs(self.DipoleRadPTN[k, :]), c=c, alpha = 0.7, label=f"{int(np.round((1E9*self.wl[k]), decimals=-1))}nm")
            ax.legend(loc=1, bbox_to_anchor=(1.2, 0.9), ncol=1, fancybox=True, shadow=True)
            # ax.set_rmin(0)
            ax.set_rmax(np.max(np.abs(self.DipoleRadPTN[:, :])))
            plt.pause(0.001)
            asdf = 1
            # input("Press Enter to Continue.")

    def update(self, frame, ax, t, lines, x, y, wl):

        # ax.set_ylim(0, 1.5E8)
        pre_line = lines[frame-1]
        pre_line.set_data(0, 0)

        line = lines[frame]
        line.set_data(x, np.abs(y[frame, :]))
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

        # ax.set_ylim(0, 15E-13)
        ymax = 1.1 * np.round(int(np.max(np.abs(self.DipoleRadPTN[:, :]))), decimals=-6)
        ax.set_ylim(0, ymax)
        ax.set_rmax(np.max(np.abs(self.DipoleRadPTN[:, :])))
        ax.legend(loc=1, bbox_to_anchor=(1.2, 0.9), ncol=1, fancybox=True, shadow=True)
        ax.set_title("15nm Radius Au NanoParticle Radiation Pattern", fontsize=20)
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.grid(True)
        ax.set_theta_zero_location("E")
        plt.show()

        t = ax.text(0, 0, "Ready", fontsize=20)

        ani = FuncAnimation(fig, self.update, frames=self.wl.size, fargs=(ax, t, self.lines, self.theta, self.DipoleRadPTN, self.wl), interval=500, blit=False, repeat=True)

        filepath = tkinter.filedialog.asksaveasfilename(initialdir=f"{fd}/",
                                                        title="Save as",
                                                        filetypes=(("GIF Files", ".gif"),
                                                                   ("all files", "*")))
        filepath = f"{filepath}.gif"

        input("Press Enter to Continue.")

        ani.save(filepath, writer='imagemagick', fps=5)

