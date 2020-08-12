function [X_HiV,Y_HiV,X_est,Y_est] = MPD(data_samp, f1, f2, FS, k)
bandfilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',f1,'CutoffFrequency2',f2,'SampleRate',FS);

data_filt = filtfilt(bandfilt,data_samp);

vz = imag(hilbert(data_filt(:,1)));
vx = data_filt(:,2);
vy = data_filt(:,3);

Y_HiV = movmean(vz.*vy,k);
X_HiV = movmean(vz.*vx,k);


B = Y_HiV\-X_HiV;
X_est = [-200:1:200];
Y_est = X_est*B;