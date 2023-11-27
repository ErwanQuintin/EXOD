import numpy as np
from astropy import wcs
from astropy.table import Table, vstack
from scipy.stats import binned_statistic_dd
from scipy.stats import poisson
from scipy.ndimage import convolve
from astropy.io import fits
from astropy.convolution import convolve
from astropy import wcs
from astropy.table import Table
from scipy.sparse import csr_matrix
from matplotlib.colors import LogNorm
#from scipy.ndimage import convolve
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import median_filter
import cmasher as cmr
import os



size_arcsec = 10 #Size of a end pixel in arcsec
pixel_size = size_arcsec/0.05 #Size of a end pixel in DetX DetY values
extent=70000
nb_pixels = int(extent/pixel_size)

data_path=os.getcwd()+'/Testing/0831790701/'

def create_fake_burst(cubeshape, time_interval, time_peak_fraction, position, width_time, amplitude):
    time = np.arange(0, cubeshape[-1])
    peak_data = np.zeros(cubeshape, dtype=int)
    time_peak = time_peak_fraction*len(time)
    width_bins = width_time / time_interval

    #Poissonian PSF
    sigma_2d = (6.6 / size_arcsec)# * 2.355  # 6.6" FWHM, 4.1" per pixel, so (6.6/4.1) pixel FWHM, 2.355 to convert FWHM in sigma
    for x in range(int(position[0]-10*sigma_2d), int(position[0]+10*sigma_2d)):
        for y in range(int(position[1]-10*sigma_2d), int(position[1]+10*sigma_2d)):
            sqdist = (x-position[0])**2+(y-position[1])**2
            psf = (1/(2*np.pi*np.sqrt(sigma_2d)))*np.exp(-(sqdist)/(2*(sigma_2d**2)))
            peak_data[x,y]+=np.random.poisson(psf*(amplitude*time_interval)*np.exp(-(time-time_peak)**2/(2*(width_bins**2))))
    #peak_data = convolve(peak_data, np.ones((3,3,1), dtype=np.int64),mode='constant', cval=0.0)
    return peak_data

def test_eef():
    position=[1000,1000]
    peak_data = create_fake_burst(cube_EPIC.shape, time_interval, time_peak_fraction=0.6, position=position,
                                  width_time=1200, amplitude=1e-1)
    all_photons = np.sum(peak_data)
    photonimage = np.sum(peak_data, axis=2)
    radii = np.linspace(0, 140, 100)
    tab_radii=[]
    for x in range(peak_data.shape[0]):
        for y in range(peak_data.shape[1]):
            for photon in range(photonimage[x][y]):
                tab_radii.append(size_arcsec*np.sqrt((x-position[0])**2+(y-position[1])**2))
    #plt.plot(np.cumsum(tab_radii, bins=radii, cumulative=True)[0]/all_photons)#Sum all photons in a given radius
    fig, (ax1,ax2,ax3)=plt.subplots(1,3)
    ax1.hist(tab_radii, bins=radii, density=True,cumulative=True)
    ax2.imshow(np.sum(cube_EPIC+peak_data,axis=2),norm=LogNorm())
    ax3.imshow(np.sum(cube_EPIC,axis=2),norm=LogNorm())

def test_image():
    # Imaging test of combined EPIC data
    raw_images = np.zeros((nb_pixels,nb_pixels))
    for inst in ('PN','M1','M2'):
        data = fits.open('/home/erwan/Documents/GitHub/EXOD/Testing/0831790701/'+inst+"_patternclean.fits")[1].data
        raw_images+=np.histogram2d(data['X'],data['Y'],bins=(np.linspace(0,extent, nb_pixels+1), np.linspace(0,extent, nb_pixels+1)))[0] #np.array(grouped_events)
    plt.figure()
    plt.imshow(raw_images, norm=LogNorm())


#3D cube test of combined EPIC data
data_pn = Table(fits.open(f'{data_path}PN_cleanpattern.fits')[1].data)['X','Y','TIME','RAWX','RAWY','CCDNR']
data_pn = data_pn[~((data_pn['CCDNR']==4)&(data_pn['RAWX']==12))&
                  ~((data_pn['CCDNR']==5)&(data_pn['RAWX']==11))&
                  ~((data_pn['CCDNR']==10)&(data_pn['RAWX']==28))] #Bad rows in Struder et al. 2001b
data_pn = data_pn[~(data_pn['RAWX']==0)&~(data_pn['RAWX']==64)&
                  ~(data_pn['RAWY']==0)&~(data_pn['RAWY']==200)] #Eject borders. Might need to adapt this to observing modes
data_M1 = Table(fits.open(f'{data_path}M1_cleanpattern.fits')[1].data)['X','Y','TIME']
data_M2 = Table(fits.open(f'{data_path}M2_cleanpattern.fits')[1].data)['X','Y','TIME']
data_EPIC = vstack((data_pn,data_M1,data_M2))

#Choose the time binning, and peak properties
time_interval=5000
peak_width=250
amplitude=1e1

#Create the data cube
n_bins = int(np.ceil((np.max(data_EPIC['TIME']) - np.min(data_EPIC['TIME'])) / time_interval))
stop_time = np.min(data_EPIC['TIME']) + n_bins * time_interval
start_time = np.min(data_EPIC['TIME'])
time_windows = np.arange(start_time, stop_time+1, time_interval)
cube_EPIC = binned_statistic_dd((data_EPIC['X'],data_EPIC['Y'],data_EPIC['TIME']), values=None, statistic='count',
                                    bins=(np.linspace(0,extent, nb_pixels+1),
                                          np.linspace(0,extent, nb_pixels+1),
                                          time_windows))[0]
#Crop the cube
size_mean_box=5
indices_image = np.where(np.sum(cube_EPIC, axis=2) > 1)
cube_EPIC=cube_EPIC[np.min(indices_image[0])-size_mean_box:np.max(indices_image[0])+1+size_mean_box,
                    np.min(indices_image[1])-size_mean_box:np.max(indices_image[1])+1+size_mean_box]
cube_EPIC[np.where(np.sum(cube_EPIC, axis=2) < 1)] = np.full(len(time_windows) - 1, np.nan)

#Add the bursting source at coordinates X,Y and T, expressed as fraction of total length
x=int(0.6*cube_EPIC.shape[0])
y=int(0.6*cube_EPIC.shape[1])
time_peak_fraction=0.6

before = np.sum(cube_EPIC[x,y])
cube_EPIC+=create_fake_burst(cube_EPIC.shape, time_interval, time_peak_fraction=time_peak_fraction, position=[x,y],
                             width_time=peak_width, amplitude=amplitude)

#Compute the background desourced by convolving with a constant buffer
k = np.ones((size_mean_box,size_mean_box,1))
k[int((size_mean_box - 1) / 2), int((size_mean_box - 1) / 2),:]=0
image = np.sum(cube_EPIC, axis=2)
threshold=np.nanpercentile(image.flatten(), 95)
space_mean_nosource = convolve(np.where(np.repeat((image<threshold)[:,:,np.newaxis],28, axis=2),
                                        cube_EPIC, np.nan), k)#, boundary='extend')

#Compare the total and desourced images
fig, (ax1,ax2) = plt.subplots(1,2)
ax1.imshow(np.sum(cube_EPIC,axis=2),norm=LogNorm(vmin=10,vmax=1000))
ax2.imshow(np.sum(space_mean_nosource,axis=2),norm=LogNorm(vmin=10,vmax=1000))
ax1.set_title("Full image")
ax2.set_title("Desourced background image")


#Now transient detections: likelihoods compared between stacked image and max over time dimension
det_likelihood_max = np.nanmax(-poisson.logpmf(cube_EPIC, space_mean_nosource), axis=2)
det_likelihood_stack = -poisson.logpmf(np.sum(cube_EPIC,axis=2), np.sum(space_mean_nosource,axis=2))

fig, (ax1,ax2) = plt.subplots(1,2)
ax1.imshow(np.nanmax(-poisson.logpmf(cube_EPIC, space_mean_nosource), axis=2), norm=LogNorm())
ax1.set_title("Likelihood of arising from the mean-estimated background")
# map_likelihoodvariable = np.nanmax(-poisson.logpmf(cube_EPIC, space_mean_nosource+np.nanmedian(cube_EPIC-space_mean_nosource, axis=2)[...,None]), axis=2)
# map_likelihoodvariable = np.where((~np.isnan(np.sum(cube_EPIC, axis=2)))&np.isinf(map_likelihoodvariable),0.1,map_likelihoodvariable)
ax2.imshow(np.maximum(det_likelihood_stack-det_likelihood_max,0.1), norm=LogNorm(vmin=0.1))
ax2.set_title("Likelihood difference between best frame and stacked image")


