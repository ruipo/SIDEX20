%-----------------------------%
% Written by: Rui Chen
% Email: ruic@mit.edu
% Last updated: 02/28/2020
%-----------------------------%

function[tdoa_mat,xps_tdoa,yps_tdoa] = tdoa_sidex(zdata,xdata,ydata,start_sample,end_sample,FS,method)

%---------------------------------------------------------------------------------------------------------------------------------%
% Calculates the time difference of arrival between different channels of the input data in the z axis direction. 
% Also calculates the time difference of arrival between data in the x/y channels and the z channels. 

% Inputs: zdata, xdata, ydata - [timeseries, chn]; 
%         start_sample, end_sample - event start and end sample number (float value);
%         FS - data sampling frequency (in Hz);
%         method - which method to use for tdoa calculation ('max' or 'envelope');
%           'max' selects the max value of the event time series and compares the location of the max peak to compute tdoa;
%           'envelope' fits an envelope to the event time series and compares the location of the envelope peak to compute tdoa;
%
% Outputs: tdoa_mat - tdoa between different channels in z-axis data (in secs);
%          xps_tdoa - tdoa between event on x-axis data and the z-axis data for each channel [tdoa, chn];
%          yps_tdoa - tdoa between event on y-axis data and the z-axis data for each channel [tdoa, chn];
%---------------------------------------------------------------------------------------------------------------------------------%

% configure x,y,z events
num_chn = size(zdata,1);

tdoa_mat = zeros(num_chn);
xps_tdoa = zeros(num_chn,1);
yps_tdoa = zeros(num_chn,1);

if end_sample > size(zdata,2)
    zevent = zdata(:,start_sample:end); 
%     xevent = xdata(:,start_sample:end);
%     yevent = ydata(:,start_sample:end);
else
    zevent = zdata(:,start_sample:end_sample); 
%     xevent = xdata(:,start_sample:end_sample);
%     yevent = ydata(:,start_sample:end_sample);
end

% hilbert method
if strcmp(method,'hilbert')

    for i = 1:num_chn
        for ii = i:num_chn

            [corr,lags] = xcorr(abs(hilbert(zevent(i,:))),abs(hilbert(zevent(ii,:))));
% 
%             if i == ii
%                 [xcorrs,xlags] = xcorr(abs(hilbert(xevent(i,:))),abs(hilbert(zevent(ii,:))));
%                 [ycorrs,ylags] = xcorr(abs(hilbert(yevent(i,:))),abs(hilbert(zevent(ii,:))));
% 
%                 xcorrs = xcorrs(xlags <= 1);
%                 xlags = xlags(xlags <= 1);
%                 ycorrs = ycorrs(xlags <= 1);
%                 ylags = ylags(xlags <= 1);
% 
%                 [pks,locs] = findpeaks(xcorrs);
%                 [~,locpx] = max(pks);
% 
%                 xps_tdoa(i) = xlags(locs(locpx));
% 
%                 [pks,locs] = findpeaks(ycorrs);
%                 [~,locpy] = max(pks);
% 
%                 yps_tdoa(i) = ylags(locs(locpy));
% 
%             end

            [~,loc] = max(corr);
            tdoa_mat(i,ii) = lags(loc);
            tdoa_mat(ii,i) = -lags(loc);
        end
    end
end

% max method
if strcmp(method,'hilbert_max')
   
    for i = 1:num_chn
        for ii = i:num_chn

            [~,lag1] = max(abs(hilbert(zevent(i,:))));
            [~,lag2] = max(abs(hilbert(zevent(ii,:))));
            
            tdoa_mat(i,ii) = (lag1-lag2);
            tdoa_mat(ii,i) = -(lag1-lag2);

        end
    end
end


% max method
if strcmp(method,'max')
   
    for i = 1:num_chn
        for ii = i:num_chn

            [~,lag1] = max(abs(zevent(i,:)));
            [~,lag2] = max(abs(zevent(ii,:)));

%             if i == ii
% 
%               [~,xlag] = max(abs(xevent(i,:)));
%               [~,ylag] = max(abs(yevent(i,:)));
% 
%               if xlag-lag1 < 0
%                   xps_tdoa(i) = xlag-lag1;
%               else
%                   continue
%               end
% 
%               if ylag-lag1 < 0
%                   yps_tdoa(i) = ylag-lag1;
%               else
%                   continue
%               end
% 
%             end
            
            tdoa_mat(i,ii) = (lag1-lag2);
            tdoa_mat(ii,i) = -(lag1-lag2);

        end
    end
end

%envelope method
if strcmp(method,'envelope')

    for i = 1:num_chn
        for ii = i:num_chn

            [corr,lags] = xcorr(envelope(abs(zevent(i,:)),55,'peaks'),envelope(abs(zevent(ii,:)),55,'peaks'));

%             if i == ii
%                 [xcorrs,xlags] = xcorr(envelope(abs(xevent(i,:)),55,'peaks'),envelope(abs(zevent(ii,:)),55,'peaks'));
%                 [ycorrs,ylags] = xcorr(envelope(abs(yevent(i,:)),55,'peaks'),envelope(abs(zevent(ii,:)),55,'peaks'));
% 
%                 xcorrs = xcorrs(xlags <= 1);
%                 xlags = xlags(xlags <= 1);
%                 ycorrs = ycorrs(xlags <= 1);
%                 ylags = ylags(xlags <= 1);
% 
%                 [pks,locs] = findpeaks(xcorrs);
%                 [~,locpx] = max(pks);
% 
%                 xps_tdoa(i) = xlags(locs(locpx));
% 
%                 [pks,locs] = findpeaks(ycorrs);
%                 [~,locpy] = max(pks);
% 
%                 yps_tdoa(i) = ylags(locs(locpy));
% 
%             end

            [~,loc] = max(corr);
            tdoa_mat(i,ii) = lags(loc);
            tdoa_mat(ii,i) = -lags(loc);
        end
    end
end

%xcorr method
if strcmp(method,'xcorr')

    for i = 1:num_chn
        for ii = i:num_chn

            [corr,lags] = xcorr(zevent(i,:),zevent(ii,:));

            [~,loc] = max(corr);
            tdoa_mat(i,ii) = lags(loc);
            tdoa_mat(ii,i) = -lags(loc);
        end
    end
end

%xcorr method
if strcmp(method,'finddelay')

    for i = 1:num_chn
        for ii = i:num_chn

            lags = finddelay(zevent(i,:),zevent(ii,:));

            tdoa_mat(i,ii) = -lags;
            tdoa_mat(ii,i) = lags;
        end
    end
end

tdoa_mat = tdoa_mat./FS;
xps_tdoa = xps_tdoa./FS;
yps_tdoa = yps_tdoa./FS;

end
