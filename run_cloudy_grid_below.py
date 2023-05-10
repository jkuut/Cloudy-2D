import os
import numpy as np
import math
import matplotlib.lines as mlines
from multiprocessing import Pool
import subprocess
from scipy import interpolate
import shutil

# Fixed parameters:
M_WD = 1.0
M_giant = 1.66666
massRatio = M_WD / M_giant
# Mass ratio is WD mass divided by giant mass
# Mass ratio in Eggletons formula is defined as giants mass divided by WD mass
# i.e. in this case as q^-1
T_star = 3000 # Donor's temp in K
R_star = 200  # radius in solar radii
# Dust condensation radius for C-rich star:
Rd = R_star/2.0 * (T_star/1500.0)**(5.0/2.0)
R_d = Rd/215.0  # convert to AU


cloudyTemplate_steady = """title theta_{theta}
sphere expanding
blackbody {temp}
luminosity total {lum}
radius {innerR}, power -2
stop radius {outerR}
cosmic ray background
background
abundances GASS
iterate to convergence
stop temperature off
stop neutral column density 24
SET PRES IONIZ 10000
set nend 3500
print line faint -2
print line sort wavelength range 3500 to 7000 
save continuum ".cont" last
save averages, file=".dep" last
column density hydrogen 1
column density hydrogen 2
column density helium 1
column density helium 2
column density helium 3
end of averages
save lines ".lines", array, units Angstrom, last
save overview ".ovr" last
save lines emissivity ".ems", last, emergent
H  1  4340.46
He 2  4685.64
H  1  4861.33
O  3  4958.91
O  3  5006.84
He 2  5411.37
He 1  5875.64
O  1  6300.30
O  1  6363.78
Fe10  6374.54
H  1  6562.81
N  2  6583.45
He 1  6678.15
S  2  6716.44
S  2  6730.82
end of lines
dlaw table
"""

cloudyTemplate_below = """title theta_{theta}
sphere expanding
blackbody {temp_BL}, luminosity = {lum}
blackbody {temp}, disk = 1e4 K
luminosity total {lum}
radius {innerR}, power -2
stop radius {outerR}
cosmic ray background
background
abundances GASS
iterate to convergence
stop temperature off
stop neutral column density 24
SET PRES IONIZ 10000
set nend 3500
print line faint -2
print line sort wavelength range 3500 to 7000 
save continuum ".cont" last
save averages, file=".dep" last
column density hydrogen 1
column density hydrogen 2
column density helium 1
column density helium 2
column density helium 3
end of averages
save lines ".lines", array, units Angstrom, last
save overview ".ovr" last
save lines emissivity ".ems", last, emergent
H  1  4340.46
He 2  4685.64
H  1  4861.33
O  3  4958.91
O  3  5006.84
He 2  5411.37
He 1  5875.64
O  1  6300.30
O  1  6363.78
Fe10  6374.54
H  1  6562.81
N  2  6583.45
He 1  6678.15
S  2  6716.44
S  2  6730.82
end of lines
dlaw table
"""



def get_TempLum(loss, a):
    massloss = np.power(10, -float(loss))
    # Note minus sign in exponent because of definition of mass ratio
    RL = 0.49 * massRatio**(-2.0/3.0) * a / (0.6 * massRatio**(-2.0/3.0) + math.log(1.0 + massRatio**(-1.0/3.0)))
    RdRL = R_d / RL
    beta_wrlof = 25.0/9.0 * massRatio**2.0 * ( -0.284 * RdRL**2 + 0.918 * RdRL - 0.234)

    # Grav constant * 1 solar mass / 1 AU / (15 km/s)2 = 3.943
    beta_bhl = 1.5/2.0*(massRatio * M_giant * 3.943/a)**2 * ( 1.0 + (1.0 + massRatio)* M_giant * 3.943/a)**(-3.0/2.0)
    maxbetabhl = 1.5/2.0*(massRatio * M_giant * 3.943/1.0)**2 * ( 1.0 + (1.0 + massRatio)* M_giant * 3.943/1.0)**(-3.0/2.0)
    beta = max(beta_wrlof, beta_bhl)
    if beta > 0.5:
        beta = 0.5
    elif beta < 0.0001:
        beta = 0.0001
    elif a < 4.0 and beta < maxbetabhl:
        beta = maxbetabhl

    massacc = beta * massloss


    Htemp, Hmdot = np.loadtxt('Hachisu/Hachisu_' + str(M_WD) + '.txt', unpack=True) 
    Hrad, Hmdot2 = np.loadtxt('Hachisu/radius_' + str(M_WD) + '.txt', unpack=True) 
    Htemp = np.power(10, Htemp)
    Hmdot = np.power(10, Hmdot)
    Hrad = np.power(10, Hrad)
    Hmdot2 = np.power(10, Hmdot2)

    
    if massacc > min(max(Hmdot), max(Hmdot2)):
        return 0, 0, 0
    
    elif massacc < max(min(Hmdot), min(Hmdot2)):
        
        # Accretion disk luminosity:
        # L = 0.5 * G * 1 Msun * 1 Msun/year / 1 Rsun
        minradius = min(Hrad)
        lum = 0.5 * 1.2e41 * massacc * M_WD / minradius
        # Accretion disk max temp = 0.488 * T_star
        temp = 4.782e5 * np.power(massacc*M_WD/minradius**3, 0.25) 
        lum2 = np.log10(lum)
        temp2 = np.log10(0.488 * temp)
        # Boundary layer temp 
        t_s = 1.214e7 * M_WD / minradius
        temp_BL = temp * np.power(t_s / temp, 0.125)
        accFlag = np.log10(temp_BL)

    else:
        f1 = interpolate.interp1d(Hmdot, Htemp, kind='linear')
        f2 = interpolate.interp1d(Hmdot2, Hrad, kind='linear')
        temp = f1(massacc)
        radius = f2(massacc)
        lum = np.power(radius, 2) * np.power(temp, 4) * 3.449e18
        lum2 = np.log10(lum)
        temp2 = np.log10(temp)
        accFlag = 1

        smoothedSeparations, smoothedMassloss, smoothedLums = np.loadtxt('Hachisu/smoothedLums_' + str(M_WD) + '.txt', unpack=True)
        newLum = 0
        for kl in range(len(smoothedLums)):
            if smoothedSeparations[kl] == a and smoothedMassloss[kl] == loss:
                newLum = smoothedLums[kl]
                #if newLum/lum < 1.001 and newLum/lum > 0.999:
                #    return 0, 0, 0
                #else:
                lum2 = np.log10(newLum)
    return temp2, lum2, accFlag


def denslaw(theta, dens, r_c):
    Rstar = R_star*0.00465 # Convert from solar radii to AU
    n_c = 10**dens
    radius = np.logspace(-2.5, 3.9, 3000, endpoint=False)  # radius in logarithmic from 10^-2 to 10^4.5 AU
    densarray = (np.power(radius, 2) + r_c**2 + 2 * radius * r_c * np.cos(theta))
    densarray = n_c * np.power(densarray/r_c**2, -1)
    temp = Rstar / (np.power(radius, 2) + r_c**2 + 2 * radius * r_c * np.cos(theta))
    temp = np.power(( 1.0 - temp), -1)
    densarr = np.log10(densarray*temp)
    for i in range(len(densarr)):
        if np.isnan(densarr[i]) or densarr[i] > 15.0:
            densarr[i] = 15.0
    for i in range(len(densarr)):
        if densarr[i] == 15.0:
            densarr[i] = (densarr[i-2] + densarr[i-1] + densarr[i] + densarr[i+1])/4.0
    radius2 = np.log10(radius * 1.4959787e13)  # transfrom to cm
    return radius2, densarr


def create_cloudy_input(theta, handle, innerR=11.0, outerR=17.0):

    loss = float(handle.split('_')[1])
    massloss = np.power(10, -loss)
    separation = float(handle.split('_')[3])

    dens = np.log10(massloss/(1.568e-16*np.power(separation, 2)))

    temp, lum, accFlag = get_TempLum(loss, separation)
    
    thetaR = np.radians(theta)
    #calculate the density law with given theta 
    radius, densarray = denslaw(thetaR, dens, separation)
        
    context = {
        "temp": temp,
        "lum": lum,
        "innerR" : innerR, 
        "outerR" : outerR,
        "theta": theta
    } 
    context_disk = {
        "temp": temp,
        "temp_BL": accFlag,
        "lum": lum,
        "innerR" : innerR, 
        "outerR" : outerR,
        "theta": theta
    } 
    

    with  open(handle + "-theta_" + str(theta) + ".in", "w") as myfile:

        if accFlag == 1: # steady burning
            myfile.write(cloudyTemplate_steady.format(**context))
        elif accFlag > 1:  # below steady
            myfile.write(cloudyTemplate_below.format(**context_disk))
        else:
            #print "shouldnt end here"
            exit()
        np.savetxt(myfile, np.c_[radius, densarray], fmt='%1.5f')
        myfile.write("end of dlaw")
    return


def run_cloudy(inputFile2):
    inputFile = inputFile2.rstrip('.in')
    subprocess.call('/afs/mpa/temp/kuuttila/Cloudy/c17.02/source/cloudy.exe -p  {}'.format(inputFile), shell=True)
    #subprocess.call('cloudy {}'.format(inputFile), shell=True)
    handle = inputFile.split('-')[0]
    endDir = handle + "/"
    endFiles = [inputFile + ij for ij in [".in", ".ems", ".out", ".ovr", ".dep", ".lines", ".cont"]]
    for fffile in endFiles:
        shutil.move(fffile, endDir+fffile)
    return inputFile



def integrate(thetas, handle, cwd):
    workdir = cwd + "/" + handle +"/"
    os.chdir(workdir)
    outputFiles = [handle + "-theta_" + str(j) + ".ems" for j in thetas]
    totals = np.zeros(15)
    th = 0
    for ffile in outputFiles:
        try:
            alldata = np.loadtxt(ffile, dtype='float')
            depth = alldata[:, 0] + 10**11
        except:
            alldata = np.loadtxt(outputFiles[th-1], dtype='float')
            depth = alldata[:, 0] + 10**11 
        anglesin = np.sin(np.radians(thetas[th]))
        dth = np.radians(np.concatenate(([1], thetas[1:] - thetas[0:-1])))
        N = 1e5
        x = np.logspace(np.log10(min(depth)*1.01), np.log10(max(depth)*0.99), N)
        dx = x[1:] - x[0:-1]
        for i in range(1, 16):
            data = alldata[:, i]      
            ff = interpolate.interp1d(depth, data, kind='linear')
            newdata = ff(x)
            SUM = np.sum( newdata[0:-1] * np.power(x[0:-1], 2) * dx ) * 2.0 * math.pi * anglesin * dth[th]
            totals[i-1] = totals[i-1] + SUM
        th = th + 1
    #print totals

    np.savetxt("RESULTS_" + handle + ".txt", totals)
    return




if __name__ == '__main__':


    massloss1 = np.around(np.linspace(4.9, 8.0, 32), 1)
    massloss2 = np.around(np.linspace(4.95, 7.95, 31), 2)
    masslossLabel = np.sort(np.concatenate((massloss1, massloss2)))
    separations = np.around(np.linspace(2, 10, 9), 0)

    handles = []
    for de in masslossLabel:
        for se in separations:
            flag1, flag2, flag3 = get_TempLum(de, se)
            if flag1 and flag2:
                if not os.path.isdir("dens_" + str(de) + "_sep_" + str(se)):
                    os.mkdir("dens_" + str(de) + "_sep_" + str(se))
                    handles.append("dens_" + str(de) + "_sep_" + str(se))

    
    #thetas1 = np.linspace(0, 120, 31)
    #thetas2 = np.linspace(122, 180, 30)
    #thetas = np.concatenate((thetas1, thetas2))
    thetas = np.linspace(0, 180, 181)

    for den in handles:
        for th in thetas:
            create_cloudy_input(th, den)


    inputFiles = [j + "-theta_" + str(i) for j in handles for i in thetas]

    flag = 0
    p = Pool(4)
    for inFile in p.imap(run_cloudy, inputFiles):
        flag = flag + 1


    cwd = os.path.dirname(os.path.realpath(__file__))            
    for dd in handles:
        integrate(thetas, dd, cwd)





