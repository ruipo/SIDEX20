
function [vg] = groupVelocityEst(data_samp, time_samp, f1, f2, FS,source_time, dist)
bandfilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',f1,'CutoffFrequency2',f2,'SampleRate',FS);

data_filt = filtfilt(bandfilt,data_samp);

hilt = hilbert(data_filt);

[~,ind] = max(abs(hilt));
time_event = time_samp(ind);

tdiff = time_event-source_time;
vg = dist/tdiff;