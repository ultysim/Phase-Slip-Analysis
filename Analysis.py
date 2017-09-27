from scipy.signal import savgol_filter as SG
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import FitRoutine


# Some constants:
#----------------------------------------------------------------
bval2014 = [5.5845,5.5840,5.5835,5.5830,5.5825,5.5820,5.5815,5.5810,5.5805,5.5800,5.5795,5.5790,5.5785,5.5780,5.5775]
bval73 = [5.7485,5.7480,5.7475,5.7470,5.7465,5.7460,5.7455,5.7450,5.7445,5.7440,5.7435,5.7430,5.7425,5.7420,5.7415,5.7410,5.7405]
bval2011 = [5.536,5.5355,5.535,5.5345,5.534,5.5335,5.533,5.5325,5.532,5.5315,5.531,5.5305,5.53,5.5295,5.529,5.5285,5.528,5.5275,5.527,5.5265,5.526,5.5255]

cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 1.0, 1.0)),

         'blue': ((0.0, 1.0, 1.0),
                  (0.5, 0.0, 0.0),
                  (1.0, 0.0, 0.0))
         }

yellowredblue = LinearSegmentedColormap('YellowRedBlue', cdict)
#-----------------------------------------------------------------

class Fit:
    def __init__(self,fit_array):
        '''
        Data structure to hold data fit values
        :param fit_array: Array of values assuming standard definition protocol
        '''
        self.freq = fit_array[0]
        self.amp = fit_array[1]
        self.phase_l = fit_array[2]
        self.phase_r = fit_array[3]
        self.height = fit_array[4]
        self.slope = fit_array[5]
        self.fit_quality = fit_array[6]
        self.pos_left = fit_array[7]
        self.pos_right = fit_array[8]
        self.delta_phi = fit_array[9]
        self.v_gate_index = fit_array[10]
        self.v0 = fit_array[11]


    def output_array(self):
        fit_array = [0 for i in range(12)]
        fit_array[0] = self.freq
        fit_array[1] = self.amp
        fit_array[2] = self.phase_l
        fit_array[3] = self.phase_r
        fit_array[4] = self.height
        fit_array[5] = self.slope
        fit_array[6] = self.fit_quality
        fit_array[7] = self.pos_left
        fit_array[8] = self.pos_right
        fit_array[9] = self.delta_phi
        fit_array[10] = self.v_gate_index
        fit_array[11] = self.v0

        return fit_array

    def plot(self,full = True):
        '''
        Plots the fit
        :param full: If True plots fit across entire data, else across computed fit width
        :return: None, generates plots
        '''
        if full:
            xpos = np.linspace(-50,0,10000)
        else:
            xpos = np.linspace(self.pos_left, self.pos_right, 10000)
        ytargl = self.height + self.amp * np.sin(self.freq * xpos + self.phase_l) + self.slope * (xpos - self.v0)
        ytargr = self.height + self.amp * np.sin(self.freq * xpos + self.phase_r) + self.slope * (xpos - self.v0)
        plt.plot(xpos, ytargl, 'k', linewidth=2)
        plt.plot(xpos, ytargr, 'r', linewidth=2)

    def yValue(self,x):
        yl = self.height + self.amp * np.sin(self.freq * x + self.phase_l) + self.slope * (x - self.v0)
        yr = self.height + self.amp * np.sin(self.freq * x + self.phase_r) + self.slope * (x - self.v0)
        return yl,yr

    def plot_select(self,parity,color,linewidth,xmin,xmax):
        xpos = np.linspace(xmin, xmax, 10000)

        if parity == 'b':
            ytargl = self.height + self.amp * np.sin(self.freq * xpos + self.phase_l) + self.slope * (xpos - self.v0)
            ytargr = self.height + self.amp * np.sin(self.freq * xpos + self.phase_r) + self.slope * (xpos - self.v0)
            plt.plot(xpos, ytargl, color[0], linewidth=linewidth[0])
            plt.plot(xpos, ytargr, color[1], linewidth=linewidth[1])
        elif parity == 'l':
            ytargl = self.height + self.amp * np.sin(self.freq * xpos + self.phase_l) + self.slope * (xpos - self.v0)
            plt.plot(xpos, ytargl, color, linewidth=linewidth)
        else:
            ytargr = self.height + self.amp * np.sin(self.freq * xpos + self.phase_r) + self.slope * (xpos - self.v0)
            plt.plot(xpos, ytargr, color, linewidth=linewidth)



class Run:
    def __init__(self, raw_data, prune_data, bvals, best_fits_array):
        self.bvals = bvals
        self.raw_data = list()
        self.prune_data = list()

        for i in range(len(raw_data)):
            self.raw_data.append(Sweep(raw_data[i], bvals[i], best_fits_array[i]))
            self.prune_data.append(Sweep(prune_data[i], bvals[i], best_fits_array[i]))

    def slip_events(self, prune=True, dV=0.025, minv=1, maxv=1, mindiff=0.03, minsep=0.2):

        if prune:
            datahold = self.prune_data
        else:
            datahold = self.raw_data

        hold = list()

        for i in datahold:
            hold.append(i.slip_events(dV, minv, maxv, mindiff, minsep))

        out = [item for sublist in hold for item in sublist]

        return out


class Sweep:

    def __init__(self,data,bval,best_fits):
        '''
        Hold a single sweep file with data and accompanying fits
        :param data: Assumes array where [x_values,y_values]
        '''
        self.x = data[0]
        self.y = data[1]
        self.b_val = bval
        self.fits = []
        for i in best_fits:
            self.fits.append(Fit(i))

    def plot(self):
        plt.plot(self.x,self.y)

    def plot_fit(self, fitnum, full=True):
        plt.plot(self.x, self.y)
        self.fits[fitnum].plot(full)

    def slip_events(self, dV=0.025, minv=1, maxv=1, mindiff=0.03, minsep=0.2):
        return FitRoutine.countEdge(self.x, self.y, self.fits, self.b_val, dV, minv, maxv, mindiff, minsep)

def generate_2014():
    run = Run(np.load('2014_Raw.npy'), np.load('2014_Prune.npy'), bval2014, np.load('52bf.npy'))
    return run


def generate_2011():
    run = Run(np.load('2011_Raw.npy'), np.load('2011_Prune.npy'), bval2011, np.load('2011bf.npy'))
    return run

def generate_73():
    run = Run(np.load('73_Raw.npy'), np.load('73_Prune.npy'), bval73, np.load('73bf.npy'))
    return run

