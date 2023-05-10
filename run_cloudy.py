import os
import numpy as np
import math
import matplotlib.lines as mlines
from multiprocessing import Pool
import subprocess
from scipy import interpolate
import shutil



cloudyTemplate = """title theta_{theta}
sphere expanding
blackbody 5.880018092416749, luminosity = 35.87719510152741
blackbody 5.117827279972283, disk = 1e4 K
luminosity total 35.87719510152741
radius {innerR}, power -2
stop radius {outerR}
cosmic ray background
background
abundances GASS 
metals 0.2 linear
iterate to convergence
stop temperature 1e3 K
stop neutral column density 25
SET PRES IONIZ 10000
print line faint -2
print line sort wavelength range 3500 to 7000 
save averages, file=".dep" last
column density hydrogen 1
column density hydrogen 2
column density helium 1
column density helium 2
column density helium 3
end of averages
#save lines ".lines", array, units Angstrom, last
save overview ".ovr" last
save continuum ".cont" last
save element oxygen ".oxy"
save element oxygen density ".oxy.den"
save lines emissivity ".ems", last, emergent
O  6  1031.91
O  6  1037.62
H  1  4340.46
He 2  4685.64
H  1  4861.33
O  3  4958.91
O  3  5006.84
Fe10  6374.54
H  1  6562.81
end of lines
dlaw table
"""




def denslaw(theta, dens, r_c=3.67 ):
    n_c = 10**dens
    radius = np.logspace(-3, 4.5, 3000, endpoint=False)  # radius in logarithmic from 10^-2 to 10^4.5 AU
    densarray = (np.power(radius, 2) + r_c**2 + 2 * radius * r_c * np.cos(theta))
    densarray = n_c * np.power(densarray/r_c**2, -1)
    temp = 0.8278 / (np.power(radius, 2) + r_c**2 + 2 * radius * r_c * np.cos(theta))
    temp = np.power(( 1.0 - temp), -1)
    densarr = np.log10(densarray*temp)
    for i in range(len(densarr)):
        if np.isnan(densarr[i]) or densarr[i] > 15.0:
            densarr[i] = 15.0
    for i in range(len(densarr)):
        if densarr[i] == 15.0:
            densarr[i] = (densarr[i-2] + densarr[i-1] + densarr[i] + densarr[i+1])/4.0
    radius2 = np.log10(radius * 1.4959787e13)  # transfrom to cm
    return radius2, densarr #return radius array in logarithmic cm and the corresponding density



def create_cloudy_input(theta, dens, innerR=11.0, outerR=16.0):
    
    thetaR = np.radians(theta)

    #calculate the density law with given theta 
    radius, densarray = denslaw(thetaR, dens)
    context = {
        "innerR" : innerR, 
        "outerR" : outerR,
        "theta": theta
    } 
    

    with  open("dens_" + str(dens) + "_theta_" + str(theta) + ".in", "w") as myfile:
        myfile.write(cloudyTemplate.format(**context))
        np.savetxt(myfile, np.c_[radius, densarray], fmt='%1.5f')
        myfile.write("end of dlaw")
    return


def run_cloudy(inputFile2):
    inputFile = inputFile2.rstrip('.in')
    subprocess.call('cloudy {}'.format(inputFile), shell=True)
    #subprocess.call('cloudy {}'.format(inputFile), shell=True)
    handle = inputFile.split('_')[1]
    endDir = "dens_" + handle + "/"
    endFiles = [inputFile + ij for ij in [".in", ".ems", ".out", ".ovr", ".dep", ".cont"]]
    for fffile in endFiles:
        shutil.move(fffile, endDir+fffile)
    return inputFile



def integrate(thetas, dens, cwd):
    workdir = cwd + "/dens_"+str(dens)+"/"
    os.chdir(workdir)
    outputFiles = ["dens_" + str(dens) + "_theta_" + str(j) + ".ems" for j in thetas]
    totals = np.zeros(9)
    th = 0
    for ffile in outputFiles:
        try:
            alldata = np.loadtxt(ffile, dtype='float')
            depth = alldata[:, 0] + 10**11
        except:
            print("Empty input file: ")# + ffile
            alldata = np.loadtxt(outputFiles[th-1], dtype='float')
            depth = alldata[:, 0] + 10**11 

        dth = np.radians(np.concatenate(([1], thetas[1:] - thetas[0:-1])))
        N = 1e5
        x = np.logspace(np.log10(min(depth)*1.01), np.log10(max(depth)*0.99), N)
        dx = x[1:] - x[0:-1]
        for i in range(1, 10):
            data = alldata[:, i]      
            ff = interpolate.interp1d(depth, data, kind='linear')
            newdata = ff(x)
            SUM = np.sum( newdata[0:-1] * np.power(x[0:-1], 2) * dx ) * 2.0 * math.pi * anglesin * dth[th]
            totals[i-1] = totals[i-1] + SUM
        th = th + 1
    print(totals)

    np.savetxt("RESULTS_" + str(dens) + ".txt", totals)
    return



if __name__ == '__main__':

    densities = np.around(np.linspace(8.3, 9.1, 17),2)
    
    for de in densities:
        if not os.path.isdir("dens_" + str(de)):
            os.mkdir("dens_" + str(de))
    
    #thetas1 = np.linspace(0, 120, 31)
    #thetas2 = np.linspace(121, 180, 60)
    #thetas = np.concatenate((thetas1, thetas2))
    thetas = np.linspace(0, 180, 181)
    
    for den in densities:
        for th in thetas:
            create_cloudy_input(th, den)    # Create cloudy inputs, with theta given in degrees
        
    inputFiles = ["dens_" + str(j) + "_theta_" + str(i) for j in densities for i in thetas]
    

    flag = 0
    p = Pool(4)
    for inFile in p.imap(run_cloudy, inputFiles):
        flag = flag + 1
        if flag % 10 == 0:
            print(inFile)

    cwd = os.path.dirname(os.path.realpath(__file__))            
    for dd in densities:
        integrate(thetas, dd, cwd)







