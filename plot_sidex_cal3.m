clear all
close all

%function[] = plot_sidex(prefix)
%prefix = '/Volumes/SIDEX/data_folder/test2/';
FS = 1000;
cal_file = '~/data/sidex/sidex_2020_01_25/GPSLOG00.TXT'
cal_fs=5;
plot_all=0;
plot_corr=0;
plotall=0;
time_file = '~/data/sidex/sidex_2020_01_25/GPSLOG_2020-01-25.txt'
midch=4;
array_gps_file ='~/data/sidex/2020_geophone_GPS.txt';
%xpos = [30; -31; 0; 31];
%ypos = [-35; 15; 0; 27];
%geophone_GPS= [71.332947,-156.407381; 71.33338,-156.409125;71.333275,-156.408282; 71.333504,-156.407474]; % Old coordinates!

% get hydrophone delays
gethdelays=0;
getallcov=1;


% correct for element mapping:
geophone_GPS=[71.332947,-156.407381;71.33338,-156.409125;71.333504,-156.407474;71.333275,-156.408282];%Rui placement

%geophone_GPS=[71.332947,-156.407381;71.333275,-156.408282;71.333504,-156.407474;71.33338,-156.409125];%Erin guess?

standalone_GPS= [71.333270, -156.408087;71.335713, -156.398203;71.334492,-156.416451;71.329938,-156.401615];

% set run date string:
rundate = '20200125';

% file to export time_file to:
time_mat = [time_file(1:end-4) '.mat'];

prefix='~/data/sidex/sidex_2020_01_25/';
directory = dir([prefix 'Sidex_' rundate '*.txt']);
ch_plot=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
[b,a]=butter(6,30/(1000/2),'low');
[d,c]=butter(6,1/(1000/2),'high');
% % create a duration vector for each file:
% for i=1:length(directory)
%     % get the timestamp, convert to duration
%     daq_time_vec(i)=directory(i).name
%     
%     filename = [prefix directory(i).name];
% end

% Read in Cal file:

%%
% Read in GPS file:

if ~isfile(time_mat)
fid = fopen(time_file);
tline = fgetl(fid);
%  Read GPGGA (code based on
%  https://github.com/balzer82/Matlab-NMEA-File-Reader/blob/master/readGPS.m)
ii=1;
while ischar(tline)
    if strncmp(tline(17:end),'$GPGGA',6)
        data = textscan(tline(17:end),'%s%f%f%c%f%c%f%f%f%f',1,'delimiter',','); 
        % compute UTC(HHMMSS.SSS), Universal Time Coordinated 
        hour = floor(data{2}/10000); 
        minute = floor((data{2}-hour*10000)/100); 
        second = round(data{2}-floor(data{2}/100)*100); 
        UTC(ii) = hours(hour)+minutes(minute)+seconds(second); 
        
        %compute latitude(DDMM.MMM) and longitude(DDDMM.MMMM) 
        lat_degree = floor(data{3}/100); 
        lat_decimal = round((data{3}-lat_degree*100)/60*10^6)/10^6; 
        lat_A = strcat(num2str(lat_degree+lat_decimal),data{4}); 

        long_degree= floor(data{5}/100); 
        long_decimal = round((data{5}-long_degree*100)/60*10^6)/10^6; 
        long_A = strcat(num2str(long_degree+long_decimal),data{6});

        numsats = data{8};
        quality = data{7};
        
        % now that it is parsed, write vector of UTC v. daqtime:
        daqtime{ii}=tline(1:15); % string that you will want to match
        ii=ii+1;
        GPS_data(ii,:)=[hour, minute, second, numsats, lat_decimal, long_decimal]; 
        
    else
        %continue
    end
    tline = fgetl(fid);
end
UTC.Format='hh:mm:ss';

% write .mat:
save(time_mat,'GPS_data','UTC','daqtime');
fclose(fid);
else
    load(time_mat);
end
duration_plot=ceil(FS);
duration_corr=ceil(FS);%seconds
%%
% locate features in cal data, for each feature, identify the string we are going for in the file name:
A=readtable(cal_file);
cal_time=table2array(A(:,1));
cal_date = table2array((A(:,2)));
cal_lat=table2array(A(:,3));
cal_long=table2array(A(:,4));

a_x = table2array(A(:,6));
a_y=table2array(A(:,7));
a_z=table2array(A(:,8));
abs_accel = a_x.^2+a_y.^2+a_z.^2;

% next, find peaks: 
[pks,locs]=findpeaks(-a_z,cal_fs,'MinPeakHeight',8,'MinPeakDistance',2);


for kk=10:length(pks)
    
    %close all
    % identify, read in and plot the files associated w/ each peak:
    [val,ind]=min(abs(UTC-cal_time(locs(kk)*cal_fs)))
    
    % find associated filenames for file of that minute and the minute
    % before to ensure you capture the event:
    daq_timestr = daqtime{ind}; 
    daq_hour = str2num(daq_timestr(end-5:end-4));
    daq_minute = str2num(daq_timestr(end-3:end-2));
    file_minute = dir([prefix 'Sidex_' daq_timestr(1:end-2) '*.txt']);
    
    if daq_minute ~= 0
        timestr_before=[daq_timestr(1:end-4) num2str(daq_minute-1,'%02i')];   
    else
        timestr_before=[daq_timestr(1:end-6) num2str(daq_hour-1) num2str(59,'%02i')];
    end
    file_before = dir([prefix 'Sidex_' timestr_before  '*.txt']);
    % find the duration from GPS associated w/ file_before:
        
    timestr_before = file_before.name(7:21);
    % next load file_minute and file_before:
    M = [dlmread([file_before.folder '/' file_before.name], ',', 2, 0);dlmread([file_before.folder '/' file_minute.name], ',', 2, 0)]';
    
    % normalize all channels:
    for ii=ch_plot
        M(ii,:)=M(ii,:)./max(M(ii,:));
    end
    
    % filter:
    Mf=filter(b,a,M');
    Mf=filter(d,c,Mf)';
    M(1:12,:)=Mf(1:12,:);
    
    
   % M(15,:)=Mf(15,:)';
    
    % create time vector:
    file_ind_time=find(contains(daqtime,timestr_before),1);
    gps_file_start = UTC(file_ind_time);
    daq_gps_file_time=UTC(file_ind_time)+seconds((0:(length(M)-1))/FS); % gps time vector for daq data!
   
    % do tdoa on the data: 
    
    [~,start_ind] = min(abs(daq_gps_file_time-cal_time(locs(kk)*cal_fs)));
    data_trunk = M(:,start_ind:(start_ind+duration_plot));
    [pp,ll]=findpeaks(data_trunk(15,:),'MinPeakHeight',0.005,'MinPeakDistance',100);
    
    start_ind=start_ind+ll(1)-100;
    
    % calc range vector for this event:
    lat_event = cal_lat(locs(kk)*cal_fs);
    long_event = cal_long(locs(kk)*cal_fs);
    
    % calculate range matrix: difference in distance sound travels to each geophone:
    ExpectedRangeMat = getLocalDist(lat_event,long_event,geophone_GPS);
    hRangeEvent(kk) = getLocalRange(lat_event,long_event,geophone_GPS(3,1),geophone_GPS(3,2));
    ExpectedRangeMat(:,5)=0;
    ExpectedRangeMat(5,:)=0;
    
    % plot the full sequence first: 
    
    
    eindd = duration_plot+start_ind;
    eind2 = duration_corr+start_ind;
    
    % try lp filter:
    
    
     figure(6)
     subplot(2,2,4)
      hold on
      
       plot(seconds(cal_time),a_z*max(abs(M(7,:)))/(2*max(a_z)))
       
   for ii=ch_plot
        plot(seconds(daq_gps_file_time),M(ii,:)+(ii))
       
        xlim([seconds(daq_gps_file_time(start_ind)),seconds(daq_gps_file_time(start_ind+duration_plot))])       
   
   end
   if plotall
    % plot data in M v. accelerometer data, starting with event:
    %calibration time:
    % next, find indices associated with the event: must occur after jump,
    % so run TDOA on the 2 s after the jump occured.
    
    
     % plot the file v. accelerometer info w/ time referenced correctly
     
     figure(3)
     clf
    subplot(2,1,1)
      hold on
   % index all times to ch0:
   z_delays=[];%4 geophones + 1 hydrophone
   x_delays=[];
   y_delays=[];
   
      
   % instead of xcorr, try something else: try for peak find.
   offset_z=0
   for ii=ch_plot
        
        if ii~=14 && ii~=15 && ii~=13
            plot(seconds(daq_gps_file_time),M(ii,:)+(ii-1)*max(abs(M(7,start_ind:eindd)))/2)
        elseif ii==13
            plot(seconds(daq_gps_file_time),M(ii,:)*max(abs(M(7,start_ind:eindd)))/max(M(ii,start_ind:eindd))+(ii-1)*max(abs(M(7,start_ind:eindd)))/2)
        else
            plot(seconds(daq_gps_file_time),M(ii,:)/100+(ii-1)*max(abs(M(7,start_ind:eindd)))/2)
        end
        xlim([seconds(daq_gps_file_time(start_ind)),seconds(daq_gps_file_time(start_ind+duration_plot))])       
        % peakfind in current channel:
        
        [~,locs2]=findpeaks(-M(ii,start_ind:eind2),FS,'NPeaks',1,'SortStr','descend','MinPeakProminence',10e-4);
        index = ceil((ii+1)/4);
        
        
        if ii==1 || ii==4 || ii==7 || ii==10 || ii==15
            if ~isempty(locs2)
            if offset_z==0 
                offset_z = locs2;
            end
        
            % z axis
            
            z_delays(index)=locs2-offset_z;
            end
        elseif ii==2 || ii == 5 || ii==8 || ii==11
            x_delays(index)=locs2-offset_x
        elseif ii==3 || ii==6 || ii==9 || ii==12
            y_delays(index)=locs2-offset_y
        end
        
            
   end
   
   subplot(2,1,2)
   plot(seconds(cal_time),a_z)
   
   xlim([seconds(daq_gps_file_time(start_ind)),seconds(daq_gps_file_time(start_ind+duration_plot))])
   %%
   figure
   % plot range vectors:
   rangeMat = getRangeVec_from0(lat_event,long_event, geophone_GPS)
   
   plot(x_delays)
   plot(y_delays)
   plot(z_delays)
   figure(10)
   hold on
   plot(rangeMat./z_delays)
    end
   
   if seconds(val) > 1
        disp('something wrong here- the files are mis-aligned!')
    else
    %M = M(:,start_ind:(start_ind+duration_plot));
    % Next, calculate delays between the channels: 
    TDOA_Z_mat=zeros(5,5);
    TDOA_X_mat=zeros(5,5);
    TDOA_Y_mat=zeros(5,5);
    
    if gethdelays
    % estimate delay for geophone 3 / hydrophone:
    x_i = (midch-1)*3+2;
    y_i=(midch-1)*3+3;
    z_i=(midch-1)*3+1;
    x_j = 15;
    y_j=15;
    z_j=15;
    [corr_x,lags_x] = xcorr(M(x_i,start_ind:eind2),M(x_j,start_ind:eind2));
    [corr_y,lags_y] = xcorr(M(y_i,start_ind:eind2),M(y_j,start_ind:eind2));
    [corr_z,lags_z] = xcorr(M(z_i,start_ind:eind2),M(z_j,start_ind:eind2));
    
    [~,loc_x] = max(abs(corr_x));
    [~,loc_y] = max(abs(corr_y));
    [~,loc_z] = max(abs(corr_z));
    TDOA_hydrophone_geophone(kk,:)= [lags_x(loc_x)/FS,lags_y(loc_y)/FS,lags_z(loc_z)/FS]
    
    end
    
    if getallcov
    
    % create distance mat:
    dist_mat = zeros(5,5)
    figure(2)
    
    for ii=1:4
        for jj=1:4
            x_i = (ii-1)*3+2;
            y_i=(ii-1)*3+3;
            z_i=(ii-1)*3+1;
            x_j = (jj-1)*3+2;
            y_j=(jj-1)*3+3;
            z_j=(jj-1)*3+1;
            
            [corr_x,lags_x] = xcorr(M(x_i,start_ind:eind2),M(x_j,start_ind:eind2));
            [corr_y,lags_y] = xcorr(M(y_i,start_ind:eind2),M(y_j,start_ind:eind2));
            [corr_z,lags_z] = xcorr(M(z_i,start_ind:eind2),M(z_j,start_ind:eind2));
            [corr_xy,lags_xy] = xcorr(M(x_i,start_ind:eind2).^2+M(y_i,start_ind:eind2).^2,M(x_j,start_ind:eind2).^2+M(y_j,start_ind:eind2).^2);
            subplot(4,4,jj+(ii-1)*4)
            hold on
            plot(lags_x/FS,corr_x./max(corr_x))
            plot(lags_y/FS,corr_y./max(corr_y))
            plot(lags_z/FS,corr_z./max(corr_z))
            plot(lags_xy/FS,corr_xy./max(corr_xy))
            [~,loc_x] = max(abs(corr_x));
            [~,loc_y] = max(abs(corr_y));
            [~,loc_z] = max(abs(corr_z));
            
            % Find TDOA for  s around peak location
            TDOA_X_mat(ii,jj)= lags_x(loc_x)/FS;
            TDOA_Y_mat(ii,jj)= lags_y(loc_y)/FS;
            TDOA_Z_mat(ii,jj)= lags_z(loc_z)/FS;
            
            indmat(ii,jj)=ii;
            indmat2(ii,jj)=jj;
        end
    end  
    %pause
   
    
    figure(1)
    
    subplot(2,2,1)    
    pcolor(TDOA_X_mat)
    title(['X delay, peak index = ' num2str(kk)])
    subplot(2,2,2)
    pcolor(TDOA_Y_mat)
    title('Y delay')
    subplot(2,2,3) 
    pcolor(TDOA_Z_mat)
    title('Z delay')
    subplot(2,2,4) 
    pcolor(ExpectedRangeMat)
    title('expected range diff')
    
    saveas(gcf,['Delays' num2str(kk) '.png'])
    
     % with built-up TDOA mat, calculate cp, cs:
    c_z = abs(ExpectedRangeMat./TDOA_Z_mat);
    c_y = abs(ExpectedRangeMat./TDOA_Y_mat);
    c_x = abs(ExpectedRangeMat./TDOA_X_mat);
    figure(6)
    subplot(2,2,1)
    pcolor(c_z)
    caxis([200,1000])
    colorbar
     title('soundspeed z (m/s)')
    subplot(2,2,2)
    pcolor(c_x)
    caxis([200,1000])
    colorbar
     title('soundspeed x (m/s)')
    
    subplot(2,2,3)
    pcolor(c_y)
    caxis([200,1000])
    colorbar
    title('soundspeed y (m/s)')
    saveas(gcf,['soundspeeds_' num2str(kk) '.png'])
    
    
    
    end
   end
end
%%
figure(1)
plot(hRangeEvent,TDOA_hydrophone_geophone,'.')


% find c for the ice:
c_water=1500;
c_x = hRangeEvent./(-hRangeEvent./c_water + TDOA_hydrophone_geophone(:,1)');
c_y= hRangeEvent./(-hRangeEvent./c_water + TDOA_hydrophone_geophone(:,2)');
c_z=hRangeEvent./(-hRangeEvent./c_water + TDOA_hydrophone_geophone(:,3)');

figure(2)
hold on
plot(c_x,'x')
plot(c_y,'.')
plot(c_z,'o')


pause

%%

% Pull the file from google drive: 


% Plot the file w/ 



N=size(directory)
pause

data = [];

stind=1;%length(directory)-10;
eind=length(directory);
startna=split(directory(stind).name,'Sidex_');
date_start=startna{2};

for i = stind:eind
    filename = [prefix directory(i).name]
    err=0;
    try 
        M = dlmread(filename, ',', 2, 0);
    catch
        disp('file error')
        M=[];
    end
    size(M)
   % pause
    if size(M,2)==16
    
    data = [data; M(1:(length(M)-1),:)];
    end
    %figure(2)
    %clf
    %hold on
    ch_plot=[1,2,3,4,5,6,7,8,9,10,11,12,14,15]%1:16;%[1,2,3,4,5,6,10,11,12,13,15,16]
   % for ii=ch_plot
    %    
     %   plot(M(:,ii)+(ii-1)*max(abs(M(:,13))))
   % end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
t = 0:1/FS:size(data,1)/FS - 1/FS;
figure(2)
hold on
[b,a]=butter(6,40/(50/2),'low');
figure(3)
hold on
for ii=ch_plot
    figure(2)
    if ii~=14 && ii~=15
   plot(t/60,data(:,ii)+(ii-1)*max(abs(data(:,1))/2))
    else
        plot(t/60,data(:,ii)/10+(ii-1)*max(abs(data(:,1))/2))
    end
   xlabel(['raw data: minutes after date: ' date_start])
   figure(3)
   
   data_filt=filter(b,a,data(:,ii));
   data_filt_scale=(data_filt-mean(data_filt))/max(abs(data_filt));
   max(data_filt)
   max(data_filt_scale)
   %pause
   plot(t/60,data_filt_scale+(ii-1)*0.5)
   xlabel(['filtered data: minutes after date: ' date_start])
end
% for j = 1:size(data,2)
% >> plot_event
%    figure
%     plot(t./60,data(:,j))
%     xlabel('time (min)')
%     ylabel('amplitude')
%     title(['chn = ', num2str(j)]);
%     set(gca,'fontsize',25);
%     grid on
% end
% %end
