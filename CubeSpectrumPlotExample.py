# Script to make publishable-quality spectrum plot of input image cube

import astropy.units as u
import numpy as np
import regions
from astropy.utils import data
from spectral_cube import SpectralCube
import pyspeckit
from astropy.io.misc import yaml
# The input yaml file contains all of the input variables needed by the script
infile = 'CubeSpectrumPlotExampleInput.yaml'

with open(infile) as fh:
    params = yaml.load(fh)

cubefile = params['cubefile'] # Cube from which the spectrum is to be extracted
regfile = params['regfile'] # Regions file containing positions through which spectra can be extracted
target = params['target'] # Target name for plot annotation
figfile = params['figfile'] # Output figure file name
velconvention = params['velconvention'] # Velocity axis convention used for spectrum
regplot = int(params['regplot']) # Which region number to produce spectrum towards
smoothfact = int(params['smoothfact']) # Smooth spectrum to smoothfact resolution in km/s
yminval = float(params['yminval']) # Intensity axis minimum value
ymaxval = float(params['ymaxval']) # Intensity axis maximum value
# Allow for empty line marker variables
if params['linenames'] != None:
    linenames = list(params['linenames'].split(", ")) # Line names to use for line ID markers
    linexvals = u.Quantity(list(map(float, params['linexvals'].split(", "))), u.MHz) # X-axis values for linenames
    linenames_sizes = np.array(list(map(int,params['linenames_sizes'].split(", ")))) # Size in points for each linename
    arrow_tips = list(map(int,params['arrow_tips'].split(", "))) # Size of arrow tips used in annotate
    # Define spectral line annotation box locations
    box_locs = [x+25 for x in arrow_tips]
vregion = u.Quantity(float(params['vregion']), u.km/u.s) # Nominal velocity of region used to derive spectrum to properly mark the line center velocity positions

cube = SpectralCube.read(cubefile) # Read cube
print(cube)

regionlist = regions.read_ds9(regfile) # Read regions file
cubespec = cube.subcube_from_regions([regionlist[regplot]]).mean(axis=(1,2)) # Define spectrum of Region regplot so that we can define the spectral axis
cubespec_mJyperbeam = cubespec.to(u.mJy) # Convert intensity units of spectrum to mJy/beam

cubespec_hdu = cubespec_mJyperbeam.hdu

restf = cubespec.wcs.wcs.restfrq*u.Hz
if velconvention == 'optical':
    offset = restf-(vregion).to(u.GHz,u.doppler_optical(restf))
elif velconvention == 'radio':
    offset = restf-(vregion).to(u.GHz,u.doppler_radio(restf))
print(offset)
# The following will define the spectral axis units using the specified velocity convention
spectral_axis = cubespec_mJyperbeam.with_spectral_unit(u.GHz, velocity_convention=velconvention, rest_value=restf).spectral_axis + offset
sp = pyspeckit.Spectrum(xarr=spectral_axis, data=cubespec_mJyperbeam.value, header={}) # Extract spectrum from cubespec_mJyperbeam
# The following will convert sp to the specified velocity convention read from the input yaml file
if velconvention == 'optical':
    sp.xarr.convert_to_unit('km/s', refX=restf, equivalencies=u.doppler_optical(restf)) # Do this first so I can smooth to smoothfact km/s spectral resolution
elif velconvention == 'radio':
    sp.xarr.convert_to_unit('km/s', refX=restf, equivalencies=u.doppler_radio(restf))
sp.smooth(smoothfact/np.abs(sp.xarr.cdelt(approx=True).value)) # Smooth to smoothfact km/s, which is smoothfact/cdelt channels
sp.xarr.convert_to_unit('GHz') # This will give me a frequency axis
sp.plotter(ymin=yminval,ymax=ymaxval)
sp.plotter.label(xlabel='Rest Frequency (GHz)',ylabel='Flux Density (mJy beam$^{-1}$)') # X- and Y-axis labeling
sp.plotter.axis.plot([sp.xarr.min().value,sp.xarr.max().value],[0,0],color='k',linestyle='-',zorder=-5) # Draw a line at y=0
if params['linenames'] != None:
    sp.plotter.line_ids(linenames,linexvals,auto_yloc=False,auto_yloc_fraction=0.83,label1_size=9,arrow_tip=arrow_tips,box_loc=box_locs) # Plot line IDs read from yaml file
# NOTE: annotate must be called *after* line_ids
sp.plotter.axis.annotate(s=target,xy=(0.05,0.9),xycoords='axes fraction') # Annotate with source name
sp.plotter.savefig(figfile) # Save to output file
