#
# Import data from NGC1300 fits image and mask foreground stars,
# crop image and save as appropriate file format for Gaia
#

import numpy as np
from astropy.io import fits
from AstroPython import Display

from matplotlib import pyplot as plot
from matplotlib import cm
import matplotlib as mpl
mpl.rcParams['figure.facecolor'] = 'w'
plot.ion()

from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interpolate


data = {
        # retrieve data from FITS files
        'ir'  : fits.getdata("Raw/ngc1300-dss2-ir.fits"),
        'red' : fits.getdata("Raw/ngc1300-dss2-red.fits"),
        'blue': fits.getdata("Raw/ngc1300-dss2-blue.fits")
    }

# model background level for all filters for each aperture
background = {

    'ir': [
        # mean,  std dev
        [4499.30, 175.587],
        [4549.20, 195.074],
        [4816.91, 203.928],
        [4906.38, 207.398],
        [5039.52, 201.842],
        [4778.03, 156.374],
        [4670.17, 182.195],
        [4801.97, 200.862],
        [6584.68, 267.167],
        [4682.94, 200.903],
        [4459.27, 199.213],
        [5019.66, 194.805],
        [4499.09, 166.825],
        [8001.50, 585.436]
    ],

    'red': [

        [5178.63, 222.264],
        [5376.69, 191.723],
        [6178.94, 282.316],
        [6466.31, 144.493],
        [6746.27, 302.845],
        [6052.45, 195.014],
        [5893.26, 185.470],
        [6637.61, 204.494],
        [9470.95, 269.433],
        [5909.86, 162.993],
        [5483.68, 175.466],
        [6379.69, 311.760],
        [5232.94, 192.232],
        [13992.8, 899.225]

    ],

    'blue': [
        [4073.02, 203.639],
        [4150.83, 226.415],
        [4531.65, 356.344],
        [5150.93, 222.480],
        [5168.56, 220.488],
        [4557.21, 249.786],
        [4422.69, 226.622],
        [4482.67, 243.598],
        [7854.97, 493.657],
        [4647.30, 219.719],
        [4268.31, 239.204],
        [4432.22, 250.525],
        [4098.84, 185.002],
        [9332.95, 2984.58]
    ]
}

regions = {

    'ir': [
        # x center,  y center,  radius
        # ----------------------------
        [640.24773, 584.59215, 7.92246],#
        [653.06269, 512.42372, 7.92246],#
        [603.82628, 534.68127, 10.6165],#
        [551.21752, 525.23867, 6.81093],#
        [543.12387, 512.42372, 9.20780],#
        [342.80589, 349.20166, 9.20780],#
        [557.96224, 337.73565, 8.92244],#
        [388.66986, 345.83264, 6.84294],#
        [415.64871, 471.95877, 6.84294],#
        [290.19705, 486.79714, 6.84294],#
        [253.10113, 493.88738, 6.84294],#
        [430.48708, 590.33678, 8.46956],#
        [603.15173, 638.89871, 11.2354],#
        [517.49388, 503.32998, 4.78753]
    ],

    'red': [
        # x center,  y center,  radius
        # ----------------------------
        [640.24773, 584.59215, 7.92246],#
        [653.06269, 512.42372, 7.92246],#
        [603.82628, 534.68127, 10.6165],#
        [551.21752, 525.23867, 6.81093],#
        [543.12387, 512.42372, 9.20780],#
        [342.80589, 349.20166, 9.20780],#
        [557.96224, 337.73565, 8.92244],#
        [388.66986, 345.83264, 6.84294],#
        [415.64871, 471.95877, 6.84294],#
        [290.19705, 486.79714, 6.84294],#
        [253.10113, 493.88738, 6.84294],#
        [430.48708, 590.33678, 8.46956],#
        [603.15173, 638.89871, 11.2354],#
        [517.49388, 503.32998, 4.78753]
    ],

    'blue': [
        # x center,  y center,  radius
        # ----------------------------
        [380.54816, 347.93060, 4.7041517],#
        [388.20763, 305.17151, 4.7041517],#
        [358.96165, 318.32767, 6.3037932],#
        [327.72501, 312.68918, 4.0441544],#
        [322.93634, 305.08948, 5.4673530],#
        [204.11658, 208.22245, 5.4673530],#
        [331.87979, 201.58613, 5.2979116],#
        [231.34754, 206.25982, 4.0631584],#
        [247.26867, 281.01793, 4.0631584],#
        [172.77320, 289.71700, 4.0631584],#
        [150.73721, 293.89279, 4.0631584],#
        [255.99006, 351.18030, 5.0290022],#
        [358.48154, 380.08577, 6.6712821],
        [307.72436, 299.68562, 2.8426906]
    ]

}

cropped = {

    'ir': [
        # parameters for `window`
        456.48938, # center x pixel
        456.79009, # center y pixel
        434.95953, # width in x
        355.24026  # width in y
    ],

    'red': [
        # parameters for `window`
        456.48938, # center x pixel
        456.79009, # center y pixel
        434.95953, # width in x
        355.24026  # width in y
    ],

    'blue': [
        # parameters for `window`
        271.53182, # center x pixel
        272.06184, # center y pixel
        258.26772, # width in x
        210.51011  # width in y
    ]
}

# create and index map of the pixels
ilength, jlength = np.shape(data['red'])
imap = np.arange(ilength)
jmap = np.arange(jlength)

# display monitor to view progress
display = Display.Monitor()

print('\n Applying masks to `ngc1300-dss2-red` ...')
# replace ds9 regions with appropriate background levels
for i in imap:
    for j in jmap:

        # display progress
        display.progress(jlength*i + j, ilength*jlength)

        for a, aperture in enumerate(regions['red']):

            x0, y0, r = aperture

            if ( j - x0 )**2 + ( i - y0 )**2 < r**2 :
                mu, sigma = background['red'][a]
                data['red'][i,j] = np.random.normal(mu, sigma)

display.complete()

# create and index map of the pixels
ilength, jlength = np.shape(data['ir'])
imap = np.arange(ilength)
jmap = np.arange(jlength)

print('\n Applying masks to `ngc1300-dss2-ir` ...')
# replace ds9 regions with appropriate background levels
for i in imap:
    for j in jmap:

        # display progress
        display.progress(jlength*i + j, ilength*jlength)

        for a, aperture in enumerate(regions['ir']):

            x0, y0, r = aperture

            if ( j - x0 )**2 + ( i - y0 )**2 < r**2 :
                mu, sigma = background['ir'][a]
                data['ir'][i,j] = np.random.normal(mu, sigma)

display.complete()

# create and index map of the pixels
ilength, jlength = np.shape(data['blue'])
imap = np.arange(ilength)
jmap = np.arange(jlength)

print('\n Applying masks to `ngc1300-dss2-blue` ...')
# replace ds9 regions with appropriate background levels
for i in imap:
    for j in jmap:

        # display progress
        display.progress(jlength*i + j, ilength*jlength)

        for a, aperture in enumerate(regions['blue']):

            x0, y0, r = aperture

            if ( j - x0 )**2 + ( i - y0 )**2 < r**2 :
                mu, sigma = background['blue'][a]
                data['blue'][i,j] = np.random.normal(mu, sigma)

display.complete()

# define edges of cropped image
xmin = {
        'ir'  : cropped['ir'][0]   - cropped['ir'][2]   / 2,
        'red' : cropped['red'][0]  - cropped['red'][2]  / 2,
        'blue': cropped['blue'][0] - cropped['blue'][2] / 2
    }

xmax = {
        'ir'  : cropped['ir'][0]   + cropped['ir'][2]   / 2,
        'red' : cropped['red'][0]  + cropped['red'][2]  / 2,
        'blue': cropped['blue'][0] + cropped['blue'][2] / 2
    }

ymin = {
        'ir'  : cropped['ir'][1]   - cropped['ir'][3]   / 2,
        'red' : cropped['red'][1]  - cropped['red'][3]  / 2,
        'blue': cropped['blue'][1] - cropped['blue'][3] / 2
    }

ymax = {
        'ir'  : cropped['ir'][1]   + cropped['ir'][3]   / 2,
        'red' : cropped['red'][1]  + cropped['red'][3]  / 2,
        'blue': cropped['blue'][1] + cropped['blue'][3] / 2
    }

# crop the images
image = {

        'ir': data['ir'][
                np.floor(xmin['ir']):np.floor(xmax['ir']),
                np.floor(ymin['ir']):np.floor(ymax['ir'])
            ],

        'red': data['red'][
                np.floor(xmin['red']):np.floor(xmax['red']),
                np.floor(ymin['red']):np.floor(ymax['red'])
            ],

        'blue': data['blue'][
                np.floor(xmin['blue']):np.floor(xmax['blue']),
                np.floor(ymin['blue']):np.floor(ymax['blue'])
            ]
    }

# new index vectors
x = {
        'ir'   : np.arange( np.shape(image['ir'])[1]   ),
        'red'  : np.arange( np.shape(image['red'])[1]  ),
        'blue' : np.arange( np.shape(image['blue'])[1] )
    }

y = {
        'ir'   : np.arange( np.shape(image['ir'])[0]   ),
        'red'  : np.arange( np.shape(image['red'])[0]  ),
        'blue' : np.arange( np.shape(image['blue'])[0] )
    }

# resample (higher density pixels)
newx, newy, xx, yy, zz, f = {}, {}, {}, {}, {}, {}

print('\n Resampling models to higher density ... `ir`, ', end='')
f['ir']    = interpolate.interp2d(x['ir'], y['ir'], image['ir'])
newx['ir'] = np.linspace(0, len(x['ir']) - 1, 1000)
newy['ir'] = np.linspace(0, len(y['ir']) - 1, 1000)
#xx['ir'], yy['ir'] = np.meshgrid(newx['ir'], newy['ir'])
zz['ir'] = f['ir'](newx['ir'], newy['ir'])

print('`red`, ', end='')
f['red']    = interpolate.interp2d(x['red'], y['red'], image['red'])
newx['red'] = np.linspace(0, len(x['red']) - 1, 1000)
newy['red'] = np.linspace(0, len(y['red']) - 1, 1000)
#xx['red'], yy['red'] = np.meshgrid(newx['red'], newy['red'])
zz['red'] = f['red'](newx['red'], newy['red'])

print('`blue`, ', end='')
f['blue']    = interpolate.interp2d(x['blue'], y['blue'], image['blue'])
newx['blue'] = np.linspace(0, len(x['blue']) - 1, 1000)
newy['blue'] = np.linspace(0, len(y['blue']) - 1, 1000)
#xx['blue'], yy['blue'] = np.meshgrid(newx['blue'], newy['blue'])
zz['blue'] = f['blue'](newx['blue'], newy['blue'])

print('done')

# physical coordinates
d = 33726 #pcs

x = {
        'ir'   : np.linspace( -d/2, d/2, len(newx['ir'])),
        'red'  : np.linspace( -d/2, d/2, len(newx['red'])),
        'blue' : np.linspace( -d/2, d/2, len(newx['blue']))
    }

y = {
        'ir'   : np.linspace( -d/2, d/2, len(newy['ir'])),
        'red'  : np.linspace( -d/2, d/2, len(newy['red'])),
        'blue' : np.linspace( -d/2, d/2, len(newy['blue']))
    }

xx['ir'],   yy['ir']   = np.meshgrid( x['ir'],   y['ir']   )
xx['red'],  yy['red']  = np.meshgrid( x['red'],  y['red']  )
xx['blue'], yy['blue'] = np.meshgrid( x['blue'], y['blue'] )

# normalize/enhance data
zz["ir"][ np.where( zz["ir"] < np.median(zz["ir"]) ) ] = 0.0
zz["ir"] = zz["ir"] ** 3
zz["ir"] /= 2 * zz["ir"].max()

zz["red"][ np.where( zz["red"] < np.median(zz["red"]) ) ] = 0.0
zz["red"] = zz["red"] ** 3
zz["red"] /= 2 * zz["red"].max()

zz["blue"][ np.where( zz["blue"] < np.median(zz["blue"]) ) ] = 0.0
zz["blue"] = zz["blue"] ** 3
zz["blue"] /= 2 * zz["blue"].max()

print('\n Saving `ir` data to csv file `ncg1300-ir.csv` ...')
output = list(zz['ir'])
output = [ [ str(value) for value in row ] for row in output ]
output = [ ', '.join(row) for row in output ]
output = [ row + '\n' for row in output ]

with open('ngc1300-ir.csv', 'w') as outfile:
    outfile.writelines(output)

print('\n Saving `red` data to csv file `ncg1300-red.csv` ...')
output = list(zz['red'])
output = [ [ str(value) for value in row ] for row in output ]
output = [ ', '.join(row) for row in output ]
output = [ row + '\n' for row in output ]

with open('ngc1300-red.csv', 'w') as outfile:
    outfile.writelines(output)

print('\n Saving `blue` data to csv file `ncg1300-blue.csv` ...')
output = list(zz['blue'])
output = [ [ str(value) for value in row ] for row in output ]
output = [ ', '.join(row) for row in output ]
output = [ row + '\n' for row in output ]

with open('ngc1300-blue.csv', 'w') as outfile:
    outfile.writelines(output)
