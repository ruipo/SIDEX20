############################### Imports ##################################################
import os
curdir = '/Users/Rui/Documents/Graduate/Research/SIDEX/SIDEX20/feature_clustering/'
os.chdir(curdir)

import numpy as np
from sidex_load import sidex_readin
import matplotlib.pyplot as plt
import matplotlib.dates as mds
import librosa
import librosa.display

path = '/Users/Rui/Documents/Graduate/Research/SIDEX/SIDEX20/sidex20_calibration/'
first_file = 0
last_file = 90
FS = 1000
NUM_CHANNELS = 16

n_fft = 1024
hop_length = int(n_fft/2)
fmax = 250 # max frequency to examine
n_mels = 128 #

data,time = sidex_readin(path,FS,NUM_CHANNELS,first_file,last_file)

t_ep_start = 1579980698.215000
time = t_ep_start+time
tcoords = mds.epoch2num(time)

tcalib = np.loadtxt(curdir+'calib_event_eptimes.txt')
tcalib_coords = mds.epoch2num(tcalib)


plt.plot(tcoords,data[:,0])
plt.plot(tcalib_coords,np.zeros(tcalib_coords.shape[0]),'r*')
ax = plt.gca()
ax.xaxis.set_major_formatter(formatter=mds.DateFormatter('%H:%M:%S'))
plt.grid()
plt.show()


num_nfft_tot = int(np.ceil(np.shape(data)[0]/hop_length)) # calculate number of total nfft bins for stats data
flist = librosa.core.mel_frequencies(n_mels=n_mels, fmin=0.0, fmax=fmax, htk=False)
tlist = (1/FS)*np.linspace(0,np.shape(data)[0],num_nfft_tot) # get list of times for spectrogram x-axis
tcoords = mds.epoch2num(tlist+t_ep_start)
S = np.zeros((n_mels,num_nfft_tot,NUM_CHANNELS)) # initialize stats data spectrogram matrix

for c in range(12): # for each channel, determine stats data mel-spectrogram
  S[:,:,c] = librosa.feature.melspectrogram(y=np.array(data[:,c]), sr=FS, n_fft=n_fft, hop_length=hop_length,fmax=fmax,n_mels=n_mels)
S = np.mean(S,axis=2)

fig = plt.figure(figsize=(20,8))

ax1 = plt.subplot(1,1,1,autoscale_on=True)
librosa.display.specshow(50*np.log10(S),x_coords=tcoords, y_coords=flist,x_axis='time',y_axis='mel', sr=FS, fmax=fmax)
plt.plot(tcalib_coords,10*np.zeros(tcalib_coords.shape[0]),'r*')
cax = plt.colorbar()
cax.set_label('dB',labelpad=-30, y=1.05, rotation=0,fontsize=15)
ax = plt.gca()
ax.xaxis.set_major_formatter(formatter=mds.DateFormatter('%H:%M:%S'))
plt.xticks(fontsize=15)
plt.yticks(np.arange(0, fmax, 192),fontsize=15)
#plt.clim(110,160)
# plt.xlabel('Time')
# plt.ylabel('Frequency (Hz)')
# plt.title('Conventional')
ax1.set_xlabel('Time',fontsize=15)
ax1.set_ylabel('Frequency (Hz)',fontsize=15)
ax1.set_title('Conventional',fontsize=15)
plt.show()

