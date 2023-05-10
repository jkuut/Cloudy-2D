import os
import numpy as np
import math





def create_out_files(inputFiles):

    outFiles = [ifile.strip('.lines') + "_OUT.lines" for ifile in inputFiles]
    nums = [str(c) for c in range(0, 10)]
    k = 0
    for infile in inputFiles:
        with open(outFiles[k], 'w') as out:
            lines = open(infile, 'r').readlines()[1:]
            try:
                for line2 in sorted(lines, key=lambda line: float(line.split()[0])):
                    line = line2.rstrip()
                    if line.endswith("t"):
                        ac = line.split()
                        if float(ac[-2]) > 30:
                            if ac[-4].endswith("A"):
                                if len(ac) == 6:
                                    lineID = ac[1]
                                    if lineID[-1] in nums:
                                        out.write(line + "\n")
                                else:
                                    out.write(line + "\n")
            except:
                print infile
        k = k + 1
    return outFiles



def create_fullOutput(density, sep, thetas, outFiles):

    wavelengths = []
    emissionLines = []
    lums = []
    ki = 0
    with open(outFiles[0], 'r') as tt:
        for line in tt:
            ac = line.rstrip().split()
            if len(ac) == 6:
                emLine = ac[1]
            elif len(ac) == 7:
                if len(ac[1]) == 2 or len(ac[2]) == 2:
                    emLine = ac[1] + " " + ac[2]
                elif len(ac[1]) == 1 and len(ac[2]) == 1:
                    emLine = ac[1] + "  " + ac[2]
            wavel = ac[-4]
        
            if ki > 0:
                if wavel != wavelengths[-1] or emLine != emissionLines[-1]:
                    wavelengths.append(wavel)
                    emissionLines.append(emLine)
                    lums.append(0.0)
                #else:
                #    print "Duplicate:" + line
            else:
                wavelengths.append(wavel)
                emissionLines.append(emLine)
                lums.append(0.0)
            
            ki = ki + 1


    i = 1
    for fil in outFiles[1:]:
        dth = np.radians(thetas[i] - thetas[i-1])
        anglesin = np.sin(np.radians(thetas[i]))
        with open(fil, 'r') as qq:
            for line in qq:
                ac = line.rstrip().split()
                wavel = ac[-4]
                lum = np.power(10, float(ac[-2]))/(4.0*3.141592)
                lum = lum * dth * anglesin * 2.0 * 3.141592
                
                if len(ac) == 6:
                    emLine = ac[1]
                elif len(ac) == 7:
                    if len(ac[1]) == 2 or len(ac[2]) == 2:
                        emLine = ac[1] + " " + ac[2]
                    elif len(ac[1]) == 1 and len(ac[2]) == 1:
                        emLine = ac[1] + "  " + ac[2]
                
                if wavel in wavelengths:
                    indx = wavelengths.index(wavel)
                    if emLine == emissionLines[indx]:
                        lums[indx] = lums[indx] + lum
                    else:
                        allIndexes = np.where(np.array(wavelengths) == wavel)[0]
                        flag = 1
                        for ii in allIndexes:
                            if emLine == emissionLines[ii]:
                                lums[ii] = lums[ii] + lum
                                flag = 0
                                break
                        if flag:
                            wavelengths.append(wavel)
                            lums.append(lum)
                            emissionLines.append(emLine)
                else:
                    wavelengths.append(wavel)
                    lums.append(lum)
                    emissionLines.append(emLine)
                
        i = i + 1

    wavelengths2 = [float(ij.strip('A')) for ij in wavelengths]
    wavelengths2, emissionLines, lums = (list(t) for t in zip(*sorted(zip(wavelengths2, emissionLines, lums))))
    
    Hbeta = 0.0
    filePath = "dens_" + str(density) + "_sep_" + str(sep)
    with open(filePath + "/fullOutput_" + filePath + ".txt", 'w') as out:
         for j in range(len(lums)):
            out.write("%.2f    " % (wavelengths2[j]) + emissionLines[j] + "    %.5e \n" % (lums[j])) 
            if wavelengths2[j] == 4861.33:
                Hbeta = lums[j]
            
    with open("lum_limit/importantOutput_" + filePath +".txt", 'w') as imp:
        for j in range(len(lums)):
            if lums[j] > 1.0e33:
                if wavelengths2[j] > 3500 and wavelengths2[j] < 7000:
                    imp.write("%.2f    " % (wavelengths2[j]) + emissionLines[j] + "    %.5e \n" % (lums[j]))
        
        
    with open("Hbeta_limit/importantOutput_" + filePath +".txt", 'w') as imp:
        for j in range(len(lums)):
            ratio = lums[j]/Hbeta
            if ratio > 0.01:
                if wavelengths2[j] > 3500 and wavelengths2[j] < 7000:
                    imp.write("%.2f    " % (wavelengths2[j]) + emissionLines[j] + "    %.5e \n" % (ratio))
                    
    for rfile in outFiles:
        os.remove(rfile)
                    
    return 



if __name__ == '__main__':


    masslossLabel = np.around(np.linspace(5.0, 6.3, 14), 1)
    separations = np.linspace(1, 14, 14)
    
    thetas1 = np.linspace(0, 120, 31)
    thetas2 = np.linspace(122, 180, 30)
    thetas = np.concatenate((thetas1, thetas2))

    for se in separations:
        for de in masslossLabel:
            filePath = "dens_" + str(de) + "_sep_" + str(se)
            if os.path.isdir(filePath):
                inputFiles = [filePath + "/" + filePath + "-theta_" + str(i) + ".lines" for i in thetas]
                outFiles = create_out_files(inputFiles)
                create_fullOutput(de, se, thetas, outFiles)
        










