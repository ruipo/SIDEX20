%-----------------------------%
% Written by: Rui Chen
% Email: ruic@mit.edu
% Last updated: 02/28/2020
%-----------------------------%

function[start_sample,end_sample,chn_triggered_mat,lta_mat,sta_mat,step_size] = sl_eSelect(data,FS,lta_size_s,sta_size_s,pem_s,pet_s,ratio_thres,weight_thres)
%---------------------------------------------------------------------------------------------------------------------------------%
% Detects large amplitude transient events in time series data based on the ratio of short-time and long-time window averages.

% INPUTS 
% data - [timeseries,chn] array
% FS - sampling frequency in hz (int)
% lta_size_s - long time window size in seconds (double)
% sta_size_s - short time window size in seconds (double)
% pem_s - pre-event buffer window in seconds (double)
% pet_s - post-event buffer window in seconds (double)
% ration_thres - ratio between sta/lta threshold for detection of event (double)
% weight_thres - minimum number of channel to be triggered to event detection - 1;

% OUTPUTS:
% start_sample - start_time of events [start sample number by number of events]
% end_time of times - end_time of events [end sample number by number of events]
% chn_triggered_mat - chns triggered for each event 
% lta_mat - long time window average for each window 
% sta_mat - short time window average for each window
% step_size - window size for each short time window
%---------------------------------------------------------------------------------------------------------------------------------%

%initialized output matrices
chn_triggered_mat = []; %initialized chn-triggered result mat
start_sample = []; %initialize start-sample result mat
end_sample = []; %initialize end-sample result mat
sta_mat = zeros(size(data)); %initialized sta result mat
lta_mat = zeros(size(data)); %initialized lta result mat

%define parameters
num_chn = size(data,1); %data(channels, samples);
lta_size = round(lta_size_s*FS); %LTA window in samples
sta_size = round(sta_size_s*FS); %STA window in samples
pem = round(pem_s*FS); %pre-event window
pet = round(pet_s*FS); %post-event window

start_delay = pem_s*FS; %pem_s seconds of delay before calculation of LTA
step_size = sta_size; %step size after each iteration

lta_start_sample = start_delay+1; %start first lta 1 sample after the start delay
lta_end_sample = lta_start_sample+lta_size; %end first lta after the lta size
sta_start_sample = start_delay+round(lta_size/2)-round(sta_size/2); %start sta at the center of the first lta
sta_end_sample = sta_start_sample+step_size; %end sta after the sta size

status = 'ntriggered'; %set initial detection status to not triggered

%start event detection
while lta_end_sample <= size(data,2)

    voting_weights = zeros(1,num_chn); %initially set the voting weight of each chn to 0
    
    % define lta, sta, find means of both
    lta_win = abs(data(:,lta_start_sample:lta_end_sample)); %define the lta window
    sta_win = abs(data(:,sta_start_sample:sta_end_sample)); %define the sta window

    lta = nanmean(lta_win,2); %find mean of lta win for each chn %size(number of chn, 1)
    lta_mat(:,sta_start_sample:sta_end_sample) = repmat(lta,1,length(sta_start_sample:sta_end_sample)); %record lta result for this window in results matrix
    
    sta = nanmean(sta_win,2);%find mean of sta window for each chn %size(number of chn, 1)
    sta_mat(:,sta_start_sample:sta_end_sample) = repmat(sta,1,length(sta_start_sample:sta_end_sample)); %record sta result for this window in results matrix
    
    % test for tigger of each chn
    ratio = sta./lta; %calculate the ratio of sta to lta
    chn_triggered = ratio > ratio_thres; %if the ratio is greater than the threshold ratio for any channel, mark that chn as triggered.
    voting_weights(chn_triggered) = 1; %set the voting weights of the triggered chn to 1.
    
    % if more than half of geophones in any direction is triggered, record the start time if the system is not triggered already.
    if sum(voting_weights)>weight_thres  
        if strcmp(status,'ntriggered')
            start_sample = [start_sample sta_start_sample-pem];
            chn_triggered_mat = [chn_triggered_mat chn_triggered]; %size(num_chn, timestep)
            
            status = 'triggered'; 
        end
    
    % if not more than half of geophones in any direction is triggered and the system is currently triggered, check to see if over half of the geophones have become untriggered, 
    % if so, untrigger the system.
    else
        
        if strcmp(status,'triggered')
            chn_untriggered = ratio < 0.85*ratio_thres; % untrigger system if the ratio drops below 75% of threshold value.
            voting_weights(chn_untriggered) = 1;
            
            if sum(voting_weights)>weight_thres % check to see if enough channels are triggered
                end_sample = [end_sample sta_end_sample+pet];
                status = 'ntriggered'; 
            end
        end   
    end
    
    % update start and end of next window.
    sta_start_sample = sta_start_sample + step_size;
    sta_end_sample = sta_end_sample + step_size;
    lta_start_sample = lta_start_sample + step_size;
    lta_end_sample = lta_end_sample + step_size;
        
end
    

end

      
    
    
    