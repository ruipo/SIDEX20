%% Read in sigicod node data
clear;clc
close all
directory = dir('/Volumes/RUIC_Backup/SIDEX20_data/SIDEx_node_data/transients_with_4_nodes_TD/*.csv');
%directory = dir('/Volumes/RUIC_Backup/SIDEX20_data/SIDEx_node_data/bad_data/*.csv');

%% Read in stand-alone geophone data
for dd = 11
Table = readtable(['/Volumes/RUIC_Backup/SIDEX20_data/SIDEx_node_data/transients_with_4_nodes_TD/' directory(dd).name], 'HeaderLines',1,'Format','%D%f%f%f%f%f%f%f%f%f%f%f%f');
%Table = readtable(['/Volumes/RUIC_Backup/SIDEX20_data/SIDEx_node_data/bad_data/' directory(dd).name], 'HeaderLines',1,'Format','%D%f%f%f%f%f%f%f%f%f%f%f%f');

%node_start = 1; %dd
%node_trim = 1; %dd

%node_start = 1; %dd1
%node_trim = 25*4096; %dd1

%node_start = 1; %dd2
%node_trim = 5*4095; %dd2

%node_start = 1; %dd3
%node_trim = 25*4096; %dd3

%node_start = 1; %dd4
%node_trim = 1; %dd4

%node_start = 20*4096; %dd5
%node_trim = 10*4096; %dd5

%node_start = 4096*13; %dd6
%node_trim = 4096*10; %dd6

%node_start = 5*4096; %dd7
%node_trim = 30*4096; %dd7

%node_start = 20*4096; %dd7
%node_trim = 15*4096; %dd7

%node_start = 0.5*4096; %dd8
%node_trim = 35*4096; %dd8

%node_start = 1; %dd9
%node_trim = 31*4096; %dd9

%node_start = 1; %dd10
%node_trim = 31*4096; %dd10

%node_start = 1; %dd11
%node_trim = 35*4096; %dd11

node_start = 6*4096; %dd11
node_trim = 25*4096; %dd11

FS_node = 4096;

tn = table2array(Table(:,1));

z_node = table2array(Table(:,[2 5 8 11]));
tempmat = ~isnan(z_node);
for i = 1:4
    tempvect1(i) = find(tempmat(:,i),1,'first');
    tempvect2(i) = find(tempmat(:,i),1,'last');
end
z_node = z_node(max(tempvect1):min(tempvect2),:);
tn = tn(max(tempvect1):min(tempvect2));

x_node = table2array(Table(:,[3 6 9 12]));
tempmat = ~isnan(x_node);
for i = 1:4
    tempvect1(i) = find(tempmat(:,i),1,'first');
    tempvect2(i) = find(tempmat(:,i),1,'last');
end
x_node = x_node(max(tempvect1):min(tempvect2),:);

y_node = table2array(Table(:,[4 7 10 13]));
tempmat = ~isnan(y_node);
for i = 1:4
    tempvect1(i) = find(tempmat(:,i),1,'first');
    tempvect2(i) = find(tempmat(:,i),1,'last');
end
y_node = y_node(max(tempvect1):min(tempvect2),:);

z_noise = z_node;
x_noise = x_node;
y_noise = y_node;

if node_trim > 0
    tn = tn(node_start:end-node_trim);
    z_node = z_node(node_start:end-node_trim,:);
    x_node = x_node(node_start:end-node_trim,:);
    y_node = y_node(node_start:end-node_trim,:);
end

figure('units','normalized','outerposition',[0 0 1 1])
for j = 1:4
    subplot(4,1,j)
    plot(tn,z_node(:,j))
    xlim([tn(1) tn(end)]);
    set(gca,'fontsize',12);
    grid on
    %ylim([-0.01 0.01])
    ylabel(['N' num2str(j)])
    set(gca,'fontsize',15)
end
%disp(directory(dd).name)
%pause
%close all
end

%% tdoa localization estimate set up
FS_node = 4096;
xpos = [];
ypos = [];
geophone_GPS = [71.332947,-156.407381; 71.33338,-156.409125; 71.333504,-156.407474; 71.333275,-156.408282];
standalone_GPS =[71.3357 -156.3982; 71.3299 -156.4016; 71.3345 -156.4165; 71.3333 -156.4081];

for g = 1:4
[xtemp, ytemp] = ll2xy(standalone_GPS(g,1),standalone_GPS(g,2),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
xpos(g) = xtemp;
ypos(g) = ytemp;
end

zdata = [];
xdata = [];
ydata = [];

for g = 1:4
    zdata(g,:) = z_node(:,g).';
    xdata(g,:) = y_node(:,g).';
    ydata(g,:) = x_node(:,g).';
end

bpFiltnode = designfilt('bandpassfir','FilterOrder',500, ...
         'CutoffFrequency1',12,'CutoffFrequency2',20, ...
         'SampleRate',FS_node);
%fvtool(bpFilt)

zdatafilt = [];
xdatafilt = [];
ydatafilt = [];
znoisefilt = [];
xnoisefilt = [];
ynoisefilt = [];

for i = 1:4
    zdatafilt(i,:) = filtfilt(bpFiltnode,z_node(:,i).');
    xdatafilt(i,:) = filtfilt(bpFiltnode,y_node(:,i).');
    ydatafilt(i,:) = filtfilt(bpFiltnode,x_node(:,i).');
    znoisefilt(i,:) = filtfilt(bpFiltnode,z_noise(:,i).');
    xnoisefilt(i,:) = filtfilt(bpFiltnode,y_noise(:,i).');
    ynoisefilt(i,:) = filtfilt(bpFiltnode,x_noise(:,i).');
end

xdatafilt_p = xdatafilt;
ydatafilt_p = ydatafilt;
xdatafilt = xdatafilt_p.*cos(3*pi/2) - ydatafilt_p.*sin(3*pi/2);
ydatafilt = xdatafilt_p.*sin(3*pi/2) + ydatafilt_p.*cos(3*pi/2);

xnoisefilt_p = xnoisefilt;
ynoisefilt_p = ynoisefilt;
xnoisefilt = xnoisefilt_p.*cos(3*pi/2) - ynoisefilt_p.*sin(3*pi/2);
ynoisefilt = xnoisefilt_p.*sin(3*pi/2) + ynoisefilt_p.*cos(3*pi/2);

znoise = zeros(1,size(znoisefilt,1));
for zz = 1:size(znoisefilt,1)
    znoise(zz) = mean(abs(znoisefilt(zz,:)));
end
xnoise = zeros(1,size(xnoisefilt,1));
for xx = 1:size(xnoisefilt,1)
    xnoise(xx) = mean(abs(xnoisefilt(xx,:)));
end
ynoise = zeros(1,size(ynoisefilt,1));
for yy = 1:size(ynoisefilt,1)
    ynoise(yy) = mean(abs(ynoisefilt(yy,:)));
end

noise_mat = [znoise; xnoise; ynoise];

figure
for j = 1:4
subplot(4,1,j)
plot(tn,zdatafilt(j,:))
hold on
plot(tn,ydatafilt(j,:));
plot(tn,xdatafilt(j,:));
title(['Geophone ' num2str(j)]);
legend('z-dir','T-dir','L-dir')
%xlabel('Time')
%ylabel('Amplitude')
%datetick('x','HH:MM:SS:FFF','keepticks');
xlim([tn(1) tn(end)]);
%title(date);
set(gca,'fontsize',12);
grid on
%k = k+1;
end

%% tdoa localization estimate ambimax
c_range = 336;
N = 1000;
plotting = 0;
plotmap = 1;
calib_act = [NaN NaN];

figure(1)
plot(xpos,ypos,'ro');
xlabel('X position (m)')
xlabel('Y position (m)')
xlim([-1050 1050]);
ylim([-1050 1050]);
xL = xlim;
yL = ylim;
line([0 0], yL,'color','black');
line(xL, [0 0],'color','black');
grid on
hold on

try
    [loc_est,c_est,err,tdoa_mat] = loc_est_ambi(zdatafilt,xdatafilt,ydatafilt,xpos,ypos,1,length(zdatafilt),FS_node,c_range,N,plotting,calib_act,noise_mat,plotmap,dd);
catch
    disp('Skipping this...')
end

figure(1)
hold on
plot(loc_est(1),loc_est(2),'k*','MarkerSize',8)
ylim([-1050 1050]);
xlim([-1050 1050]);
set(gca,'fontsize',20)
xlabel('x-position (m)')
ylabel('y-position (m)')
title(['Propagation Speed estimate: ' num2str(c_est) ' m/s']);
disp(loc_est)

%% loc_est_todaMatch

FS = 4096;
c_range = 336;
xrange = [-1000 1000];
yrange = [-1000 1000];
ds = 1;
plotting = 1;
calib_act = [NaN NaN];

%try
    [loc_est,c_est,err,rms_err_mat] = loc_est_tdoaMatch(zdatafilt,xpos,ypos,1,length(zdatafilt),FS,xrange,yrange,ds,c_range,plotting,calib_act,dd);
%catch
%    disp('Skipping this...')
%end

disp(['x estimate = ', num2str(loc_est(1)), '; y estimate = ',num2str(loc_est(2)) '.']);
disp(['propagation speed estimate = ', num2str(c_est), 'm/s.']);
disp(['estimate error = ', num2str(err), '.']);


%% Plot results
load('2020_locations.mat')
xpos = [];
ypos = [];
crackxpos = zeros(length(crack_long),1);
crackypos = crackxpos;
geophone_GPS = [71.332947,-156.407381; 71.33338,-156.409125; 71.333504,-156.407474; 71.333275,-156.408282];
standalone_GPS =[71.3357 -156.3982; 71.3299 -156.4016; 71.3345 -156.4165; 71.3333 -156.4081];

for g = 1:4
[xtemp, ytemp] = ll2xy(standalone_GPS(g,1),standalone_GPS(g,2),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
xpos(g) = xtemp;
ypos(g) = ytemp;
end

for cc = 1:length(crack_lat)
[xtemp, ytemp] = ll2xy(crack_lat(cc),crack_long(cc),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
crackxpos(cc) = xtemp;
crackypos(cc) = ytemp;
end

[crackxpos,inds] = sort(crackxpos);
crackypos = crackypos(inds);

%x_est = [129 -195 305 138 49 -159 241 128 365 227 1238 356 354];
%y_est = [-193 -327 -74 391 241 374 -43 -192 -62 168 -275 -53 -54];

%ambi_max method
%x_est = [402 -194 305 -922 49 73 31 241 75 238 -360 313 355 351];
%y_est = [-246 -325 -74 -225 239 217 -38 -43 -23 -39 -476 -72 -53 -55];

%tdoaMatch method
x_est = [217 -181 1000 1000 29 1000 523 354 296 296 611 744 736];
y_est = [-210 -217 -189 1000 116 1000 -84 -45 -47 113 -122 -99 -101];

figure
plot(xpos,ypos,'k.', 'MarkerSize',30);
xlabel('X position (m)')
xlabel('Y position (m)')
xlim([-1050 1050]);
ylim([-1050 1050]);
xL = xlim;
yL = ylim;
grid on
hold on
plot(x_est(1:4), y_est(1:4), 'gd','MarkerSize',10,'MarkerFaceColor','g');
plot(x_est(5:end), y_est(5:end), 'bp','MarkerSize',10,'MarkerFaceColor','b');
plot(crackxpos,crackypos,'r','linewidth',1.5);
set(gca,'fontsize',20)
xlabel('X-Position (m)')
ylabel('Y-Position (m)')
line([0 0], yL,'color','black');
line(xL, [0 0],'color','black');
legend('Node Geophones','01/26-02/05 Events','03/17 Events','Crack in Ice Cover')
axis square
