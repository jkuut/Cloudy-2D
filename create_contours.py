import os
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize as opt
from math import sqrt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from scipy import integrate
import matplotlib.cm as cm
import scipy.ndimage as ndimage
from matplotlib.colors import LogNorm


def hist2d(ax, var1, var2, var3, binsize, fluxline):

    smooth = 0.0
    colorbar_flag = 1
    levels = [1e32]
   # binsize = [14,14]
    hdata, xedges, yedges = np.histogram2d(var1, var2, bins=binsize, weights=var3)
    count, xedges2, yedges2 = np.histogram2d(var1, var2, bins=binsize) 

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    hdata_masked = np.ma.masked_where(hdata <= 0.0, hdata)
    count_masked = np.ma.masked_where(count <= 0.0, count)
    hdata_masked = hdata_masked/count_masked

    if smooth == 0:
        hdata_smooth = hdata_masked
    else:
        hdata_smooth = ndimage.gaussian_filter(hdata_masked, sigma=smooth, order=0)

    hdata_smooth = np.ma.masked_where(hdata <= 0.0, hdata_smooth)
    hdata_smooth = hdata_smooth.T

    im = ax.imshow(hdata_smooth,
                    interpolation='none',
                    origin='lower',
                    extent=extent,
                    cmap='seismic',
                    vmin=1e30,
                    vmax=1e34,
                    aspect='auto',
                    norm=LogNorm())
    icount = 0
    for level in levels:
        cs = ax.contour(hdata_smooth,
                        levels,
                        colors = 'k',
                        linestyles = ['--'],
                        linewidths = 2.0,
                        origin='lower',
                        extent=extent)
        zc = cs.collections[0]
        plt.setp(zc, linewidth=1)
        p = cs.collections[0].get_paths()[0]
        v = p.vertices
        xx = v[:,0]
        yy = v[:,1]
        ff = open('figs/contour_' + fluxline + '.txt', 'w')

        for kj in range(len(xx)):
            ff.write(str(yy[kj]) + " " + str(xx[kj]) + "\n")
        ff.close()

        icount += 1
    return im





if __name__ == '__main__':


    masslossLabel = np.around(np.linspace(4.9, 6.3, 15), 1)
    separations = np.linspace(2, 14, 13)

    handles = []
    zerohandles = []

    for se in separations:
        for de in masslossLabel:
            filePath = "dens_" + str(de) + "_sep_" + str(se)
            if os.path.exists("lum_limit/importantOutput_" + filePath + ".txt"):
                handles.append(filePath)
            #else:
            #    zerohandles.append(filePath)
                
                
                
    wavelengths = []
    lineIDs = []
    occs = []
    
    for i in range(len(handles)):
        f = open("lum_limit/importantOutput_" + handles[i] + ".txt", 'r')
        for line in f:
            ac = line.rstrip().split()
            if len(ac) == 3:
                emLine = ac[1]
            elif len(ac) == 4:
                if ac[1] == "H" or ac[1] == "He":
                    continue
                else:
                    emLine = ac[1] + " " + ac[2]
                    
            wavel = ac[0]
            
            if not wavel in wavelengths:
                wavelengths.append(wavel)
                lineIDs.append(emLine)
                occs.append(1)
            else:
                indx = wavelengths.index(wavel)
                occs[indx] += 1
                
    wavelengths = np.array(wavelengths)
    lineIDs = np.array(lineIDs)
    occs = np.array(occs)
    
    limFactor = 80
    wavelengths = wavelengths[np.where(occs > limFactor)]
    lineIDs = lineIDs[np.where(occs > limFactor)]
    occs = occs[np.where(occs > limFactor)]
    
    lineLums = np.zeros((len(handles)+len(zerohandles), len(occs))) + 1e29

    #print lineLums[:, 0] means one column, i.e. one line


    seps = []
    mass = []
    ii = 0
    for i in range(len(handles)):
        aaa = handles[i].split('_')
        seps.append(float(aaa[3]))
        mass.append(-float(aaa[1]))
        f = open("lum_limit/importantOutput_" + handles[i] + ".txt", 'r')
        for line in f:
            ac = line.rstrip().split()
            wavel = ac[0]
            lum = ac[-1]
            if wavel in wavelengths:
                indx = np.where(wavelengths == wavel)[0][0]
                lineLums[i, indx] = lum 
        f.close()
        ii = i
        
    for j in range(len(zerohandles)):
        aaa = zerohandles[j].split('_')
        seps.append(float(aaa[3]))
        mass.append(-float(aaa[1]))
                
                
    #wavelengths = wavelengths[::-1]
    #lineIDs = lineIDs[::-1]
    #lineLums = lineLums[::-1]

    doneIDs = []
    for a in range(len(wavelengths)):   
        if not lineIDs[a] in doneIDs:
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            im1 = hist2d(ax1, seps, mass, lineLums[:, a], [len(separations), len(masslossLabel)], wavelengths[a] + lineIDs[a])
            s1 = str(int(np.around(float(wavelengths[a]), 0))) + " " + lineIDs[a]
            cb1 = plt.colorbar(im1, ax=ax1, label=s1 + r' Line luminosity (erg s$^{-1}$ cm$^{-2}$)')
            fig.savefig("figs/" + s1.replace(' ', '_') + ".jpg", dpi=200)
            plt.close()
            doneIDs.append(lineIDs[a])
                




                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
