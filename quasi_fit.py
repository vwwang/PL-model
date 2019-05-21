import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets

def quasi(Cbx,tXr,T,tTnr,tTr,Nb):   #quasi is for interact function at the end
    G = []
    NxY = []
    Ghigh = 110        #high bound of Generation rate 10^(Ghighx0.2)
    Glow = 75          #low bound of Generation rate 10^(Glowx0.2)
    Nxh = [8*10**14]   #high bound of exciton concentration used for brute force algorithm
    Nxl = [0]          #low bound of exciton concentration used for brute force algorithm
    tT = 1/(1/tTnr+1/tTr)
    Generation = [10**(0.2*p) for p in range(Glow,Ghigh)]
    
    for g in range(len(Generation)):    #start of brute force for getting exciton concentration
        for i in range(20):
            Nx = [t for t in np.arange(Nxl[0], Nxh[0], (Nxh[0] - Nxl[0])/10)]
            Nx.append(Nxh[0])
            EValue = [t/tXr+t*t*Cbx+(T*t/(1+T*t))*(Nb/tT) - Generation[g] for t in Nx]
            abEValue = [abs(t) for t in EValue]
            minindex = [np.argmin(abEValue)]
            if minindex[0] == 10:
                Nxh = [Nxh[0]]
                Nxl = [Nxl[0] + (minindex[0]-1) * (Nxh[0] - Nxl[0]) / 10]
            elif minindex[0] == 0:
                Nxh = [Nxl[0] + (minindex[0]+1) * (Nxh[0] - Nxl[0]) / 10]
                Nxl = [Nxl[0]]
            else:
                Nxh = [Nxl[0] + (minindex[0] + 1) * (Nxh[0] - Nxl[0]) / 10]
                Nxl = [Nxl[0] + minindex[0] * (Nxh[0] - Nxl[0]) / 10]

        NxY.append(Nx[minindex[0]])
        Nxh = [8 * 10 ** 14]
        Nxl = [0]
    QY = [(NxY[i]/tXr+(T*NxY[i]*Nb/(1+T*NxY[i])/tTr))/Generation[i]*100 for i in range(len(NxY))]
    plt.figure(figsize=(11, 5))
    plt.subplot(121)
    plt.plot(Generation, QY, 'go', markersize=5)
    plt.yscale('log');
    plt.ylim(0.01, 100);
    plt.xscale('log');
    plt.xlabel('Generation rate (cm^-2)', size = 15);
    plt.ylabel('PL QY (%)', size = 15);
    
    RXr = [NxY[i]/tXr for i in range(len(NxY))]
    RXX = [NxY[i]*NxY[i]*Cbx for i in range(len(NxY))]
    RTr = [(T*NxY[i]*Nb/(1+T*NxY[i])/tTr) for i in range(len(NxY))]
    RTnr = [(T*NxY[i]*Nb/(1+T*NxY[i])/tTnr) for i in range(len(NxY))]
    plt.subplot(122)
    plt.plot(Generation, RXr, 'rs', label = 'Xr')
    plt.plot(Generation, RXX, 'r--', label = 'XX')
    plt.plot(Generation, RTr, 'bs',label = 'Tr')
    plt.plot(Generation, RTnr, 'b--', label = 'Tnr')
    plt.yscale('log');
    plt.xscale('log');
    plt.xlabel('Generation rate (cm^-2)', size = 15);
    plt.ylabel('Recombination rate (cm^-2)', size = 15);
    plt.legend(framealpha=1, frameon=True, prop={'size': 13});
    plt.show();

def run_widget():
    widgets.interact(quasi, 
         Cbx = widgets.FloatText(value=0.06,description='Cbx (cm^2/s):',disabled=False),
         tXr = widgets.FloatText(value=1.9*10**-9,description='tXr (s):',disabled=False),
         tTr = widgets.FloatText(value=1*10**-6,description='tTr (s):',disabled=False),
         tTnr = widgets.FloatText(value=1*10**-9,description='tTnr (s):',disabled=False),
         T = widgets.FloatText(value=3.5*10**-10,description='T(cm^2):',disabled=False),
         Nb = widgets.FloatText(value=2*10**10,description='N (cm^-2):',disabled=False));

if __name__ == '__main__':
    run_widget()