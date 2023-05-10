# This script creates some useful files for an observing run with WiFeS
#
# - a map with the instrument fields
# - a list of the fields position to be read by TAROS
#

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20,})
rc('text', usetex=False)

import aplpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pyfits
from scipy import ndimage


#from matplotlib import rcParams
#rcParams['font.family']='sans-serif'
#rcParams['font.sans-serif']=['Helvetica']
#plt.rcParams['font.size']=13

# define some useful functions
def dectodeg (dec) :
    return np.sign(dec[0])*(abs(dec[0]) + dec[1]/60. + dec[2]/3600.)

def ratodeg (ra) :
    return (ra[0] + ra[1]/60. + ra[2]/3600.)/24.*360

# Converte degrees back to hh/min/ss
def degtora (deg) :
    aa = deg/360.0 * 24.0
    ra1 = int(aa)
    mins = (aa - ra1)*60.0
    ra2 = int(mins)
    ra3 = ((aa - ra1)*60.0 - ra2)*60.0
    ra3 = np.round(ra3, 4)
    return [ra1, ra2, ra3]

def myrotmatrix(angle):
    rotmatrix = np.zeros((2,2))
    rotmatrix[0,0] = np.cos(np.radians(angle))
    rotmatrix[0,1] = -np.sin(np.radians(angle))
    rotmatrix[1,0] = np.sin(np.radians(angle))
    rotmatrix[1,1] = np.cos(np.radians(angle))
    return rotmatrix


# Calculate the RA of the second point using the cosine law
# cos(angsep) = cos(90 - dec1) * cos(90 - dec2) + sin(90 - dec1) * sin(90 - dec2) * cos(ra1 - ra2)
# Angsep = angular separation of the two points, given in arcseconds. Positive angsep means RA2 > RA1
def cosinelaw(angsep, dec1, dec2, ra1):
    twopi = 2.0 * 3.14159265358979
    pitwo = 3.14159265358979 / 2.0
    de1 = dectodeg(dec1)/360.0*twopi
    de2 = dectodeg(dec2)/360.0*twopi
    ra1 = ratodeg(ra1)/360.0*twopi
    angsep = angsep/3600.0/360.0*twopi
    cosDa = (np.cos(angsep) - np.cos(pitwo - de1)*np.cos(pitwo - de2)) / (np.sin(pitwo - de1)*np.sin(pitwo - de2))
    if angsep > 0:
        ra2 = ra1 + np.arccos(cosDa)
    else:
        ra2 = ra1 - np.arccos(cosDa)
    ra2 = degtora(ra2/twopi*360.0)
    return ra2



# define the instrument field of view
wifes_field = [25,38]
wifes_corners = np.zeros((4,2))
wifes_corners[0] = [wifes_field[0]/2.,wifes_field[1]/2.]
wifes_corners[1] = [-wifes_field[0]/2.,wifes_field[1]/2.]
wifes_corners[2] = [-wifes_field[0]/2.,-wifes_field[1]/2.]
wifes_corners[3] = [wifes_field[0]/2.,-wifes_field[1]/2.]

muse_corners = np.zeros((4,2)) # From the MUSE manual, march 2014
muse_corners[0] = [0.0086683878,   0.0084439565]
muse_corners[1] = [359.99128-360., 0.008499508]
muse_corners[2] = [359.99172-360.,-0.0082782215]
muse_corners[3] = [0.0080579666,  -0.008389310]
muse_corners *= 3600 # make it in arcsec

# ----------- USER INPUT --------------------
# Make it pretty with overlaid color image
#do_pretty = True
do_pretty = False


### Which source to use?
# Sources: 
#   CAL87
#   RXJ0513.9-6951
#   RXJ0439.8-6809
#   RXJ0527.8-6954
#   RXJ0537.7-7034
#   1E0035.4-7230
#   RXJ0048.4-7332
#   CAL83 (4 fields)
#   Lin358
source = 'CAL87'


# First, load the DSS image of the target
# Update dss2-r file here
dssname = './' + source + '/' + source + '.fits'

# Orientation
#alpha = RA; delta = DEC 
#alpha0 and delta0 for centering finder chart
#alpha1, delta1 to alphaN, deltaN are centers for the N fields of the mosaic

fields_location = {} # centers of the field

if source == 'CAL87':
    sourceid = 'CAL87'
    alpha0 = [05.0, 46.0, 46.54]
    delta0 = [-71.0, 08.0, 53.9]

    delta1 = [-71.0, 09.0, 03.9]
    alpha1 = cosinelaw(-17, delta0, delta0, alpha0) #first parameter is angular separation in arcseconds (positive or negative)

    delta2 = [-71.0, 08.0, 43.9]
    alpha2 = cosinelaw(-17, delta0, delta0, alpha0)

    delta3 = [-71.0, 09.0, 03.9]
    alpha3 = cosinelaw(17, delta0, delta0, alpha0)

    delta4 = [-71.0, 08.0, 43.9]
    alpha4 = cosinelaw(17, delta0, delta0, alpha0)


    ref_star = ['ref_'+sourceid,[05.0, 46.0, 50.1],[-71.0, 10.0, 16.67]]
    fields_location.__setitem__(1,['CAL87_1',alpha1,delta1,90,'wifes'])
    fields_location.__setitem__(2,['CAL87_2',alpha2,delta2,90,'wifes'])
    fields_location.__setitem__(3,['CAL87_3',alpha3,delta3,90,'wifes'])
    fields_location.__setitem__(4,['CAL87_4',alpha4,delta4,90,'wifes'])

elif source == 'CAL83': 
    sourceid = 'CAL83'
    alpha0 = [05.0, 43.0, 34.16]
    delta0 = [-68.0, 22.0, 22.19]


    delta1 = [-68.0, 22.0, 32.69]
    alpha1 = cosinelaw(-17, delta1, delta1, alpha0)

    delta2 = [-68.0, 22.0, 11.69]
    alpha2 = cosinelaw(-17, delta2, delta2, alpha0)

    delta3 = [-68.0, 22.0, 32.69]
    alpha3 = cosinelaw(17, delta3, delta3, alpha0)

    delta4 = [-68.0, 22.0, 11.69]
    alpha4 = cosinelaw(17, delta4, delta4, alpha0)

    ref_star = ['ref_'+sourceid,[05.0, 43.0, 21.49],[-68.0, 22.0, 15.79]]
    fields_location.__setitem__(1,['CAL83_1',alpha1,delta1,90,'wifes'])
    fields_location.__setitem__(2,['CAL83_2',alpha2,delta2,90,'wifes'])
    fields_location.__setitem__(3,['CAL83_3',alpha3,delta3,90,'wifes'])
    fields_location.__setitem__(4,['CAL83_4',alpha4,delta4,90,'wifes'])


elif source == 'RXJ0513.9-6951': 
    sourceid = source.split('-')[1]
    alpha0 = [05.0, 13.0, 50.81]
    delta0 = [-69.0, 51.0, 47.30]
    alpha1 = cosinelaw(-15, delta0, delta0, alpha0)
    delta1 = [-69.0, 51.0, 47.30]
    ref_star = ['ref_'+sourceid,[05.0, 13.0, 39.726],[-69.0, 51.0, 23.13]]
    fields_location.__setitem__(1,[sourceid,alpha1,delta1,90,'wifes'])

elif source == 'RXJ0439.8-6809': 
    sourceid = source.split('-')[1]
    alpha0 = [04.0, 39.0, 49.64]
    delta0 = [-68.0, 09.0, 01.4]
    alpha1 = cosinelaw(-15, delta0, delta0, alpha0)
    delta1 = [-68.0, 09.0, 01.4]
    ref_star = ['ref_'+sourceid,[04.0, 39.0, 52.546],[-68.0, 08.0, 40.243]]
    fields_location.__setitem__(1,[sourceid,alpha1,delta1,90,'wifes'])

elif source == 'RXJ0527.8-6954':  
    sourceid = source.split('-')[1]
    alpha0 = [05.0, 27.0, 48.6]
    delta0 = [-69.0, 54.0, 02.0]
    alpha1 = cosinelaw(-15, delta0, delta0, alpha0)
    delta1 = [-69.0, 54.0, 02.0]
    ref_star = ['ref_'+sourceid,[05.0, 28.0, 01.95],[-69.,54.,57.498]]
    fields_location.__setitem__(1,[sourceid,alpha1,delta1,90,'wifes'])

elif source == 'RXJ0537.7-7034':  
    sourceid = source.split('-')[1]
    alpha0 = [05.0, 37.0, 43.0]
    delta0 = [-70.0, 34.0, 15.0]
    alpha1 = cosinelaw(15, delta0, delta0, alpha0)
    delta1 = [-70.0, 34.0, 15.0]
    ref_star = ['ref_'+sourceid,[05.0, 37.0, 52.475],[-70.0, 34.0, 41.03]]
    fields_location.__setitem__(1,[sourceid,alpha1,delta1,90,'wifes'])

elif source == '1E0035.4-7230': 
    sourceid = source.split('-')[1]
    alpha0 = [00.0, 37.0, 19.0]
    delta0 = [-72.0, 14.0, 14.0]
    alpha1 = cosinelaw(-15, delta0, delta0, alpha0)
    delta1 = [-72.0, 14.0, 14.0]
    ref_star = ['ref_'+sourceid,[00.,37.,23.395],[-72.,14.,31.579]]
    fields_location.__setitem__(1,[sourceid,alpha1,delta1,90,'wifes'])

elif source == 'RXJ0048.4-7332': 
    sourceid = source.split('-')[1]
    alpha0 = [00.0, 48.0, 20.02]
    delta0 = [-73.0, 31.0, 52.07]
    alpha1 = cosinelaw(-15, delta0, delta0, alpha0)
    delta1 = [-73.0, 31.0, 52.07]
    ref_star = ['ref_'+sourceid,[00.0, 48.0, 22.945],[-73.0, 30.0, 56.90]]
    fields_location.__setitem__(1,[sourceid,alpha1,delta1,90,'wifes'])



elif source == 'LIN358':
    sourceid = 'LIN358'
    alpha0 = [00.0, 59.0, 12.25]
    delta0 = [-75.0, 05.0, 17.6164]
    alpha1 = cosinelaw(15, delta0, delta0, alpha0) #first parameter is angular separation in arcseconds (positive or negative)
    delta1 = [-75.0, 05.0, 17.6164]
    ref_star = ['ref_'+sourceid,[00.0, 59.0, 20.86],[-75.0, 05.0, 34.854]]
    fields_location.__setitem__(1,[sourceid,alpha1,delta1,90,'wifes'])

else:
    print "INVALID SOURCE\n"
    exit()



# Observing fields :



# 


# ------------ start the program here ----------------------

# Try to fix the header of the DSS fits file
f = pyfits.open(dssname)
scidata = f[0].data
sciheader = f[0].header
f.close()

dss_fixed_hdu = pyfits.PrimaryHDU(scidata, header = sciheader)
outfits = pyfits.HDUList([dss_fixed_hdu])
# SAVE IT
outfits.writeto('./' + source + '/' + source + '_headfixed.fits', output_verify='fix',
                        clobber=True)

dssname = './' + source + '/' + source + '_headfixed.fits'


# Defines the field points for the Polygon function

fields = [] #it's a list

for i in fields_location :

    if fields_location[i][4] == 'wifes':
        corners = np.zeros_like(wifes_corners)+wifes_corners
    elif fields_location[i][4] == 'muse':
        corners = np.zeros_like(muse_corners)+muse_corners
    else:
        print ' Instrument unknown:',fields_location[i][4]   
    
    # Rotate the different corners by the correct amount
    for (j,corner) in enumerate(corners):    
        corners[j] = np.dot(corners[j],myrotmatrix(fields_location[i][3]))
        
    # Converts from arc sec to degrees/hours
    corners = corners /3600.
    
    # Finds the field center
    thisalpha = ratodeg(fields_location[i][1])
    thisdelta = dectodeg(fields_location[i][2])
    
    # Accounts for the declination in the RA size of the field
    corners[:,0] = corners[:,0]/np.cos(np.radians(thisdelta))

    #Store the field corners
    fields.append([])
    fields[i-1].append(corners[0] + [thisalpha,thisdelta])
    fields[i-1].append(corners[1] + [thisalpha,thisdelta])
    fields[i-1].append(corners[2] + [thisalpha,thisdelta])
    fields[i-1].append(corners[3] + [thisalpha,thisdelta])
    
    fields[i-1] = np.array(fields[i-1])

try:
    fields_xray = np.zeros((len(x_ray_zone),2))
    for (k,point) in enumerate(x_ray_zone):
        fields_xray[k] = np.array([ratodeg(point[1]),dectodeg(point[2])])
except:
    pass


# Start the plotting
plt.close(2)
fig1 = plt.figure(2, figsize = (9,7.5))
if not(do_pretty):
    fig1.suptitle(source)

ax1 = aplpy.FITSFigure(dssname, figure = fig1, north=True)
ax1.show_grayscale(invert = True, stretch = 'log')

ax1.add_grid()
ax1.grid.set_color('black')
ax1.set_tick_color('k')
#ax1.set_tick_xspacing(60./3600)

if do_pretty:
    ax1.add_scalebar(30./3600) #Length in degrees
    ax1.scalebar.show(30./3600., label='30"', corner = 'top left', frame= True)
else:
    ax1.add_scalebar(30./3600) #Length in degrees
    ax1.scalebar.show(30./3600., label='30"', corner = 'top left', frame= True)
    
ax1.show_markers(ratodeg(ref_star[1]),dectodeg(ref_star[2]),
                     marker='o',edgecolor='red',s=100)
if not(do_pretty):
    ax1.add_label(ratodeg(ref_star[1])-0.01,
    dectodeg(ref_star[2])-0.003,ref_star[0],color='r')


ax1.show_polygons(fields, edgecolor='r', linewidth = 3, linestyle= '-')
#try:
#    ax1.show_polygons([fields_xray], edgecolor=(0,1,0), linewidth = 2,ls='--')
#except:
#    pass


# Plot marker for the source on top of the fields
ax1.show_markers(ratodeg(alpha0),dectodeg(delta0),
                     marker='o',edgecolor='magenta',s=100)
if not(do_pretty):
    ax1.add_label(ratodeg(alpha0)-0.01,
    dectodeg(delta0)-0.003,source,color='m')

    
ax1.set_system_latex(False)

if not(do_pretty):
    for i in fields_location :
        ax1.add_label(ratodeg(fields_location[i][1]),
                    dectodeg(fields_location[i][2]),fields_location[i][0], 
                    color='b')
    
        this_label = fields_location[i][0] + ': '+ \
            str(int(fields_location[i][1][0])) + \
            'h' + str(int(fields_location[i][1][1])) +\
            'm' + str(fields_location[i][1][2]) + \
            's ; ' + str(int(fields_location[i][2][0])) + \
            r'$^{\circ}$' + str(int(fields_location[i][2][1])) + \
            "'" + str(fields_location[i][2][2]) + \
            '" ; P.A.='+str(fields_location[i][3]) + r'$^{\circ}$'
        if i < 10 :
            ax1.add_label(0.6,1.075-i/32.,this_label,relative=True,color = 'k',size=12)
        else :
            ax1.add_label(0.6,0.48-i/32.,this_label,relative=True,color = 'k',size=12)  

# I want to add the X-ray/HST image to this
#img = mpimg.imread('SNR0104-72.3.jpg')
#img = mpimg.imread('SNR0104-72.3_Ha.png')
img = mpimg.imread('E0102.jpg')
img_size = np.array(np.shape(img)[:2])
dss_size = np.array(ax1.image.get_size())
dss_center = dss_size/2.+1
#scale_factor = 1.02
#cheat = [19,-2]
#scale_factor = 0.36
#cheat = [22,-24]
scale_factor = 0.102
cheat = [-2.8,3.]
#img_rot=ndimage.rotate(img,2,mode='constant',cval=100)
#img_rot=ndimage.rotate(img,0,mode='constant',cval=100)
img_rot=ndimage.rotate(img,4,mode='constant',cval=100)
if do_pretty:
    for i in fields_location :
        ax1.add_label(ratodeg(fields_location[i][1]),
                    dectodeg(fields_location[i][2]),fields_location[i][0], 
                    color='b')
    plt.imshow(img_rot, alpha=0.5, extent=[dss_center[1]-img_size[1]*scale_factor/2.+cheat[0],
                                        dss_center[1]+img_size[1]*scale_factor/2.+cheat[0],
                                        dss_center[0]-img_size[0]*scale_factor/2.+cheat[1],
                                        dss_center[0]+img_size[0]*scale_factor/2.+cheat[1],])
else:
    plt.imshow(img_rot, alpha=0.0, extent=[dss_center[1]-img_size[1]*scale_factor/2.+cheat[0],
                                        dss_center[1]+img_size[1]*scale_factor/2.+cheat[0],
                                        dss_center[0]-img_size[0]*scale_factor/2.+cheat[1],
                                        dss_center[0]+img_size[0]*scale_factor/2.+cheat[1],])

# Plot some orientation arrows for this one
if do_pretty:
    aloc = np.array([ratodeg(alpha0)+17./3600./np.cos(np.radians(dectodeg(delta0))),
                     dectodeg(delta0)+66./3600])
    ax1.show_arrows([aloc[0], aloc[0]], [aloc[1],aloc[1]], [12./3600/np.cos(np.radians(aloc[1])),0], [0,12./3600], 
                    head_length=4,head_width=2, color='k')
    ax1.add_label(aloc[0],aloc[1]+16./3600.,'N', verticalalignment='center',horizontalalignment='center',
                  color = 'k', weight='heavy', size='large',
                  bbox=dict(facecolor='w',ec='none', alpha=0.75, boxstyle="round,pad=.3"))
    ax1.add_label(aloc[0]+16./3600./np.cos(np.radians(aloc[1])),aloc[1],'E', verticalalignment='center',horizontalalignment='center',
                  color = 'k', weight = 'heavy', size='large',
                  bbox=dict(facecolor='w',ec='none', alpha=0.75, boxstyle="round,pad=.3"))
else:
    aloc = np.array([ratodeg(alpha0)+68./3600./np.cos(np.radians(dectodeg(delta0))),
                     dectodeg(delta0)+54./3600])
    ax1.show_arrows([aloc[0], aloc[0]], [aloc[1],aloc[1]], [12./3600/np.cos(np.radians(aloc[1])),0], [0,12./3600], 
                    head_length=4,head_width=2, color='k')
    ax1.add_label(aloc[0],aloc[1]+16./3600.,'N', verticalalignment='center',horizontalalignment='center',
                  color = 'k', weight='heavy', size='large',
                  bbox=dict(facecolor='w',ec='none', alpha=0.75, boxstyle="round,pad=.3"))
    ax1.add_label(aloc[0]+16./3600./np.cos(np.radians(aloc[1])),aloc[1],'E', verticalalignment='center',horizontalalignment='center',
                  color = 'k', weight = 'heavy', size='large',
                  bbox=dict(facecolor='w',ec='none', alpha=0.75, boxstyle="round,pad=.3"))

# Zoom in a bit ...
if do_pretty:
    ax1.recenter(ratodeg(alpha0),dectodeg(delta0),radius=86./3600)
    plt.savefig('./' + source + '/' + source + '_pretty.eps')
    plt.savefig('./' + source + '/' + source + '_pretty.pdf')
    #    plt.savefig(dssname[:-5]+'_map.eps')
#    plt.savefig('./69-0_target_map.png')
#    plt.savefig('./69-0_target_map.eps')
else :
#    plt.savefig('./69-0_target_map.png')
#ax1.recenter(ratodeg(alpha0),dectodeg(delta0),radius=45./3600)
    plt.savefig('./' + source + '/' + source + '_map.png')#,bbox_inches='tight')
    plt.savefig('./' + source + '/' + source + '_map.pdf')#,bbox_inches='tight')
fig1.show()

if not(do_pretty):
    #Now, let's create a target file tobe read in TAROS
    file = open('./' + source + '/' + source + '_taros_list','w')
#    file = open('./taros_list','w')
    file.write('# This files contains the fields location for Taros\n')
    file.write('# Field_ID RA(J2000) Dec(J2000) PA\n#\n')

    this_label = '"' + ref_star[0] + '" '+ str(int(ref_star[1][0])) + \
        ' ' + str(int(ref_star[1][1])) +\
        ' ' + str(ref_star[1][2]) + \
        ' ' + str(int(ref_star[2][0])) + \
        ' ' + str(int(ref_star[2][1])) + \
        " " + str(ref_star[2][2]) + '\n'

    file.write(this_label)

    for i in fields_location : 
        this_label = '"' + fields_location[i][0] + '" '+ str(int(fields_location[i][1][0])) + \
            ' ' + str(int(fields_location[i][1][1])) +\
            ' ' + str(fields_location[i][1][2]) + \
            ' ' + str(int(fields_location[i][2][0])) + \
            ' ' + str(int(fields_location[i][2][1])) + \
            " " + str(fields_location[i][2][2]) + ' ! '+str(fields_location[i][3]) + '\n'

        file.write(this_label)

    file.close()

    #Let's also create a jskycalc fields file

    file = open('./' + source + '/' + source + '_jskycalc_list','w')

    for i in fields_location : 
        this_label = fields_location[i][0] + ' ' +\
            str(int(fields_location[i][1][0])) + \
            ' ' + str(int(fields_location[i][1][1])) +\
            ' ' + str(fields_location[i][1][2]) + \
            ' ' + str(int(fields_location[i][2][0])) + \
            ' ' + str(int(fields_location[i][2][1])) + \
            " " + str(fields_location[i][2][2]) + ' 2000\n'

        file.write(this_label)

    file.close()
