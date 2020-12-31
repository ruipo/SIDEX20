function [X_HiV,Y_HiV,X_est,Y_est] = MPD(data_samp, k)
%bandfilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',f1,'CutoffFrequency2',f2,'SampleRate',FS);

%data_filt = filtfilt(bandfilt,data_samp);
data_filt = data_samp;
%vz = gradient(data_filt(:,1));
vz = imag(hilbert(data_filt(:,1)));
vx = data_filt(:,2);
vy = data_filt(:,3);

wins = floor(length(vz)/k);
Y_HiV = zeros(wins,1);
X_HiV = Y_HiV;
vzy = vz.*vy;
vzx = vz.*vx;

Y_HiV = movmean(vzy,k);
X_HiV = movmean(vzx,k);

% for w = 1:wins
%     Y_HiV(w) = mean(vzy(w*k-k+1:w*k));
%     X_HiV(w) = mean(vzx(w*k-k+1:w*k));
% end

xcal = X_HiV;%*cos(deg2rad(theta))-Y_HiV*sin(deg2rad(theta));
ycal = Y_HiV;%X_HiV*sin(deg2rad(theta))+Y_HiV*cos(deg2rad(theta));

% B = xcal\ycal;
% X_est = [-1000:4:999];
% Y_est = X_est*B;

X_est = linspace(-1500,1500,500);
p = polyfit(xcal,ycal,1);
Y_est = polyval(p,X_est);
