
%------------------------------------------------%
% Written by: Rui Chen & Erin Fischell
% Email: ruic@mit.edu
% Last updated: 02/28/2020
% Calibration data analysis code for SIDEX 2020
%------------------------------------------------%

%% Set Input Files
clear all
close all

FS = 1000;
cal_file = './GPSLOG00.TXT';
cal_fs=5;
time_file = './GPSLOG_2020-01-25.txt';

array_gps_file ='./2020_geophone_GPS.txt';
geophone_GPS = [71.332947,-156.407381; 71.33338,-156.409125; 71.333504,-156.407474; 71.333275,-156.408282];

% set run date string:
rundate = '20200125';

% file to export time_file to:
time_mat = [time_file(1:end-4) '.mat'];

prefix='./sidex20_calibration/';
directory = dir([prefix 'Sidex_' rundate '*.txt']);
ch_plot=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];

%% Read in GPS file:
if ~isfile(time_mat)
fid = fopen(time_file);
tline = fgetl(fid);
% Read GPGGA (code based on https://github.com/balzer82/Matlab-NMEA-File-Reader/blob/master/readGPS.m)
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
        continue
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

%% Get acceleration features and location of each calibration event (z-axis)

A=readtable(cal_file);
cal_time=table2array(A(:,1));
cal_date = table2array((A(:,2)));
cal_lat=table2array(A(:,3));
cal_long=table2array(A(:,4));

a_x = table2array(A(:,6)); %x_accel
a_y=table2array(A(:,7));   %y_accel
a_z=table2array(A(:,8));   %z_accel
abs_accel = a_x.^2+a_y.^2+a_z.^2;

a_z_f = a_z;
a_z_f(abs(a_z)<5)=0;

a_z_f(1:200) = 0;
a_z_f(650:700) = 0;
a_z_f(2550:2600) = 0;
a_z_f(3600:3750) = 0;
a_z_f(4150:4200) = 0;
a_z_f(4650:4800) = 0;
a_z_f(5450:5500) = 0;
a_z_f(6025:6100) = 0;
a_z_f(6400:6700) = 0;
a_z_f(6850:7000) = 0;
a_z_f(7250:7300) = 0;
a_z_f(7450:7600) = 0;
a_z_f(9700:9850) = 0;
a_z_f(10000:10050) = 0;
a_z_f(10350:10450) = 0;
a_z_f(10700:10800) = 0;
a_z_f(10980:11010) = 0;
a_z_f(11900:12000) = 0;
a_z_f(12500:12800) = 0;
a_z_f(12950:13000) = 0;
a_z_f(14100:17500) = 0;
a_z_f(17650:17700) = 0;
a_z_f(18000:18200) = 0;
a_z_f(18980:19000) = 0;
a_z_f(19400:19900) = 0;
a_z_f(20300:20400) = 0;
a_z_f(21450:21500) = 0;
a_z_f(21650:21900) = 0;
a_z_f(22250:22600) = 0;
a_z_f(22800:23000) = 0;
a_z_f(23150:23200) = 0;
a_z_f(23850:23900) = 0;
a_z_f(24450:24500) = 0;
a_z_f(24700:25500) = 0;

% next, find peaks in z_accel: 
[pks,locs] = findpeaks(a_z_f,'MinPeakDistance',6);
locs_start = zeros(length(locs),1);
for ll = 1:length(locs)
    val = pks(ll);
    loc_val = locs(ll);
    val0 = a_z_f(loc_val);
    while (val ~= 0) || (val0 ~=0)
        loc_val = loc_val-1;
        val = a_z_f(loc_val);
        val0 = a_z_f(loc_val-1);
    end
    
    locs_start(ll) = loc_val;
end

locs_start(49) = [];
locs_start(149) = [];
locs_start(149) = [];
locs_start(149) = [];

lat_event = zeros(length(locs_start),1);
long_event = zeros(length(locs_start),1);
UTC_event = NaT(length(locs_start),1);

for ll = 1:length(locs_start)
        
    lat_event(ll) = cal_lat(locs_start(ll)); %event lat
    long_event(ll) = cal_long(locs_start(ll)); %event long
    UTC_event(ll) = datetime([rundate ' ' char(cal_time(locs_start(ll)))],'InputFormat', 'yyyyMMdd HH:mm:ss.SSS','Format','d-MMM-y HH:mm:ss.SSS'); %event UTC time
end

epochtime_event = posixtime(UTC_event); %event epoch time

figure; plot(long_event,lat_event,'r*')
hold on
plot(geophone_GPS(:,2),geophone_GPS(:,1),'ko')

%% Calculate distance of calibration events to each geophone

calib_dist = zeros(size(geophone_GPS,1), length(lat_event));
for gg = 1:size(geophone_GPS,1)
    for ee = 1:length(lat_event)
        calib_dist(gg,ee) = 1000*lldistkm(geophone_GPS(gg,:),[lat_event(ee) long_event(ee)]);
    end
end

%% Get calibration data from directory and plot calibration events and the geophone time series

N = size(directory)

data = [];

stind=1;
eind=length(directory);
start_time = datetime([directory(stind).name(7:14) num2str(str2num(directory(stind).name(16:17))+5) num2str(str2num(directory(stind).name(18:21))+4) '.' directory(stind).name(23:25)],'InputFormat','yyyyMMddHHmmss.SSS');
start_time = posixtime(start_time);
startna=split(directory(stind).name,'Sidex_');
date_start=startna{2};

for i = stind:eind
    i
    filename = [prefix directory(i).name];
    err=0;
    try 
        M = dlmread(filename, ',', 2, 0);
    catch
        disp('file error')
        M=[];
    end

    if size(M,2)==16
    
    data = [data; M(1:(length(M)-1),:)];
    end

    ch_plot=[1,2,3,4,5,6,7,8,9,10,11,12,14,15];

end

t = 0:1/FS:size(data,1)/FS - 1/FS;
t = t+start_time; % Calibration Data epoch time

figure %(for geophone 1)
plot(t,data(:,1))
hold on
plot(t,data(:,4))
plot(t,data(:,7))
plot(t,data(:,10))
plot(epochtime_event,zeros(length(epochtime_event),1),'*')

%% Event detection in z-axis

mult_vec = zeros(length(t),1);

for ee = 1:4:length(epochtime_event)
    [~,ind1] = min(abs(t-epochtime_event(ee)));
    [~,ind2] = min(abs(t-epochtime_event(ee+3)));
    
   mult_vec(ind1-FS:ind2+4*FS) = 1;
end
    
mult_vec_mat = repmat(mult_vec,1,4);  

zdata = (mult_vec_mat.*data(:,[1 4 7 10])).';
xdata = (mult_vec_mat.*data(:,[2 5 8 11])).';
ydata = (mult_vec_mat.*data(:,[3 6 9 12])).';

lta_size_s = 5; %LTA window in seconds
sta_size_s = 0.2; %STA window in seconds
pem_s = 0.5;
pet_s = 0.5;
ratio_thres = 2;
weight_thres = 0;

znum_chn = size(zdata,1);
[zstart_sample,zend_sample,zchn_triggered_mat,~,~,~] = sl_eSelect(zdata,FS,lta_size_s,sta_size_s,pem_s,pet_s,ratio_thres,weight_thres);

found_events_lat = zeros(length(zstart_sample),1);
found_events_long = zeros(length(zstart_sample),1);
found_events_t = zeros(length(zstart_sample),1);

for en = 1:length(zstart_sample)
   tevent = t(zstart_sample(en));
   
   [~,loct] = min(abs(tevent-epochtime_event));
   found_events_lat(en) = lat_event(loct);
   found_events_long(en) = long_event(loct);
   found_events_t(en) = epochtime_event(loct);
    
end

[found_events_x,found_events_y] = ll2xy(found_events_lat,found_events_long,geophone_GPS(1,1),geophone_GPS(1,2));

figure
plot(xpos,ypos,'ro')
hold on
plot(found_events_x,found_events_y,'k*')

    
%% Plot detected events if desired

n = 1; % <------ SET DETECTED EVENT NUMBER

t_u_chn = zchn_triggered_mat(:,n);

figure
for chn = 1:znum_chn
   subplot(znum_chn,1,znum_chn-chn+1)
   title(date)

   if t_u_chn(chn) == 1 %triggered and used
       color = 'r';
   elseif t_u_chn(chn) == 0 %not triggered and not used
       color = 'k';
   end

   plot(t(zstart_sample(n):zend_sample(n)),10E3*zdata(chn,zstart_sample(n):zend_sample(n)),[color,'-.']);
   hold on
   plot(t(zstart_sample(n):zend_sample(n)),10E3*xdata(chn,zstart_sample(n):zend_sample(n)),['b-.']);
   plot(t(zstart_sample(n):zend_sample(n)),10E3*ydata(chn,zstart_sample(n):zend_sample(n)),['b-.']);
   plot(t(zstart_sample(n):zend_sample(n)),1*ratio_thres*ones(1,length(zstart_sample(n):zend_sample(n))),'k-.','linewidth',1.5)
   xlabel('Time')
   ylabel('Amplitude')
   xlim([t(zstart_sample(n)) t(zend_sample(n))])
   set(gca,'fontsize',15);
   grid on
end


%% Event localization

xpos = [0;-62.068;-3.310;-32.067];
ypos = [0;48.147;61.936;36.472];

zdatatest = [];
xdatatest = [];
ydatatest = [];

bpFilt = designfilt('bandpassfir','FilterOrder',500, ...
         'CutoffFrequency1',8,'CutoffFrequency2',32, ...
         'SampleRate',FS);
%fvtool(bpFilt)

for i = 1:4
    zdatafilt(i,:) = filtfilt(bpFilt,zdata(i,:));
    xdatafilt(i,:) = filtfilt(bpFilt,xdata(i,:));
    ydatafilt(i,:) = filtfilt(bpFilt,ydata(i,:));
end

num_event = min(length(zstart_sample),length(zend_sample));
FS = 1000;
c_range = 50:1:1000;
N = 1000;
plotting = 0;

figure(1)
plot(xpos,ypos,'ro');
xlabel('X position (m)')
xlabel('Y position (m)')
xlim([-350 350]);
ylim([-150 350]);
xL = xlim;
yL = ylim;
line([0 0], yL,'color','black');
line(xL, [0 0],'color','black');
grid on
hold on

pos_est_mat = zeros(2,length(zstart_sample));
c_est_mat = zeros(length(zstart_sample),1);
error_mat = zeros(length(zstart_sample),1);

parfor i = [1:1:162]
    i 
    
    ss = zstart_sample(i);
    es = zend_sample(i);
    
    calib_act = [found_events_x(i) found_events_y(i)];
    
    try
        [loc_est,c_est,err,tdoa_mat] = loc_est_calib(zdatatest,xdatatest,ydatatest,xpos,ypos,ss,es,FS,c_range,N,plotting,calib_act);
    catch
        disp('Skipping this')
    end
    
    pos_est_mat(:,i) = loc_est;
    c_est_mat(i) = c_est;
    error_mat(i) = err;

    disp(['x estimate = ', num2str(loc_est(1)), '; y estimate = ',num2str(loc_est(2)) '.']);
    disp(['propagation speed estimate = ', num2str(c_est), 'm/s.']);
    disp(['estimate error = ', num2str(err), '.']);
    
    if c_est < 300
        color = 'r';
    elseif (300 <= c_est) && (c_est < 600)
        color = 'b';
    elseif (600 <= c_est) && (c_est < 1000)
        color = 'g';
    else
        color = 'm';
    end
    
    figure(1);
    plot(loc_est(1),loc_est(2),[color '*']);
   
end

%% Location Estimation plot
path = '/Users/Rui/Documents/Graduate/Research/SIDEX/sidex20/'; % set path to where you want to save figures; create folder "calib_results" at this path location;

figure('units','normalized','outerposition',[0 0 1 1])
plot(xpos,ypos,'ro', 'MarkerFaceColor', 'r');
xlabel('X position (m)')
ylabel('Y position (m)')
xlim([-350 350]);
ylim([-150 350]);
xL = xlim;
yL = ylim;
line([0 0], yL,'color','black');
line(xL, [0 0],'color','black');
grid on
hold on
set(gca,'fontsize',20);

h(1) = plot(NaN,NaN,'ro','MarkerFaceColor', 'r');
h(2) = plot(NaN,NaN,'ks','MarkerSize',8);
h(3) = plot(NaN,NaN,'r*','MarkerSize',6);
h(4) = plot(NaN,NaN,'b*','MarkerSize',6);
h(5) = plot(NaN,NaN,'m*','MarkerSize',6);
h(6) = plot(NaN,NaN,'g*','MarkerSize',6);
legend(h, 'Geophone Locations','True Location','Estimated Location (v < 250m/s)','Estimated Location (250m/s <= v < 500m/s)','Estimated Location (500m/s <= v < 750m/s)','Estimated Location (v >= 750m/s)');
set(0,'DefaultLegendAutoUpdate','off')

for i = 1:length(found_events_t)
    i
    c_est = c_est_mat(i);
    if c_est < 250
        color = 'r';
    elseif (250 <= c_est) && (c_est < 500)
        color = 'b';
    elseif (500 <= c_est) && (c_est < 750)
        color = 'm';
    else
        color = 'g';
    end
    
    title([datestr(datetime(found_events_t(i), 'convertfrom','posixtime')) ' UTC'])

    plot(found_events_x(i),found_events_y(i),'ks','MarkerSize',8);
    plot(pos_est_mat(1,i),pos_est_mat(2,i),[color '*'],'MarkerSize',6);
    quiver(pos_est_mat(1,i),pos_est_mat(2,i),found_events_x(i)-pos_est_mat(1,i),found_events_y(i)-pos_est_mat(2,i),0,'color',[1 0.5 1]);
    saveas(gcf,[path,'calib_results/' num2str(i) '.png']);
    pause(0.1)
end

%% Location Estimation Error plot

distfc = sqrt((found_events_x-xpos(4)).^2 + (found_events_y-ypos(4)).^2);

figure
semilogy(distfc,error_mat,'b.','MarkerSize',10)
grid on
set(gca,'fontsize',20)
hold on
f=fit(distfc,error_mat,'poly3','Normalize','on','Robust','on');
plot(f)
xlabel('Distance from Center Geophone (m)')
ylabel('Error in Location Estimation (m)')
title('Estimation Error vs. Distance from Center Geophone')

%% Propagation Speed Estimation plot

figure
cf = fit([pos_est_mat(1,:).' pos_est_mat(2,:).'],c_est_mat,'lowess');
figure
plot(cf,[pos_est_mat(1,:).' pos_est_mat(2,:).'],c_est_mat)

[cx,cy] = meshgrid(-350:1:350, -150:1:350);
c_fit = zeros(size(cx));
for i = 1:size(cx,1)
    c_fit(i,:) = cf([cx(i,:).' cy(i,:).']);
end

c_fit(c_fit>1000) = NaN;
c_fit(c_fit<50) = NaN;
figure
h = pcolor(cx(1,:),cy(:,1),c_fit);
caxis([50 1000])
set(h,'Edgecolor','None')
hold on
plot(xpos,ypos,'ro', 'MarkerFaceColor', 'r');
xlabel('X position (m)')
ylabel('Y position (m)')
xlim([-350 350]);
ylim([-150 350]);
xL = xlim;
yL = ylim;
line([0 0], yL,'color','black');
line(xL, [0 0],'color','black');
plot(pos_est_mat(1,:),pos_est_mat(2,:),'k*')
title('Propagation Speed Estimation')
colormap(flipud(jet))
colorbar
set(gca,'fontsize',20)

%% Location Estimation Results movie

testf = dir([path 'calib_results/*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([path 'calib_results/calib_results.avi']);
  writerObj.FrameRate = 8;
  open(writerObj);
  
  for frame = 1:length(testf)
    disp([num2str(frame) '/' num2str(length(testf))])
    % Construct an output image file name.
    outputBaseFileName = sprintf(testf(frame).name);
    outputFullFileName = fullfile(testf(frame).folder, outputBaseFileName);
    % Read the image in from disk.
    thisFrame = imread(outputFullFileName);
    % Convert the image into a "movie frame" structure.
    recalledMovie(frame) = im2frame(thisFrame);
    % Write this frame out to a new video file.
    writeVideo(writerObj, thisFrame);
  end
  close(writerObj);

%% Calibration Event plot
path = '/Users/Rui/Documents/Graduate/Research/SIDEX/sidex20/'; % set path to where you want to save figures; create folder "calib_figs" at this path location;
figure('units','normalized','outerposition',[0 0 1 1])
plot(geophone_GPS(:,2),geophone_GPS(:,1),'ko','linewidth',2)
hold on

for ll = 1:length(locs_start)
ll
plot(long_event(ll),lat_event(ll),'r*')

set(gca,'fontsize',25)
grid on
xlabel('Longitude')
xlim([-156.42 -156.39])
ylabel('Latitude')
ylim([71.332 71.337])
title(['Time: ' , char(cal_time(locs_start(ll))), ' UTC'])
saveas(gcf,[path,'calib_figs/' num2str(ll) '.png']);
pause(0.1)
end

%% Calibration Events Video

testf = dir([path 'calib_figs/*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([path 'calib_figs/calib_events.avi']);
  writerObj.FrameRate = 8;
  open(writerObj);
  
  for frame = 1:length(testf)
    disp([num2str(frame) '/' num2str(length(testf))])
    % Construct an output image file name.
    outputBaseFileName = sprintf(testf(frame).name);
    outputFullFileName = fullfile(testf(frame).folder, outputBaseFileName);
    % Read the image in from disk.
    thisFrame = imread(outputFullFileName);
    % Convert the image into a "movie frame" structure.
    recalledMovie(frame) = im2frame(thisFrame);
    % Write this frame out to a new video file.
    writeVideo(writerObj, thisFrame);
  end
  close(writerObj);
