%% Plot Event time series

%t1 = 58169;%1/26
%t2 = 65651;%1/26

%t1 = 21000; %1/27_1 
%t2 = 25000; %1/27_1

t1 = 14000; %1/27_2 
t2 = 16000; %1/27_2

%t1 = 32000;%1/28_1
%t2 = 33500;%1/28_1

%t1 = 8500;%1/28_2
%t2 = 10000;%1/28_2

%t1 = 15600;%1/31
%t2 = 30400;%1/31

xdata = data(t1:t2,[2 5 8 11]).*cos(3*pi/2) - data(t1:t2,[3 6 9 12]).*sin(3*pi/2);
ydata = data(t1:t2,[2 5 8 11]).*sin(3*pi/2) + data(t1:t2,[3 6 9 12]).*cos(3*pi/2);
zdata = data(t1:t2,[1 4 7 10]);
t = t(t1:t2);

figure
k=1;
for j = [1 2 3 4]
    subplot(4,1,k)
    plot(t,zdata(:,j))
    hold on
    ylabel(['G' num2str(k)])
    set(gca,'fontsize',12);
    grid on
    set(gca,'fontsize',15)
    k = k+1;
end

%% tdoa localization estimate set up
load('/Users/Rui/Documents/Graduate/Research/SIDEX/SIDEX20/new_analysis/2020_locations.mat')
FS = 1000;
xpos = [0;-62.068;-3.310;-32.067];
ypos = [0;48.147;61.936;36.472];
xcrack = [];
ycrack = [];

geophone_GPS = [71.332947,-156.407381; 71.33338,-156.409125; 71.333504,-156.407474; 71.333275,-156.408282];
standalone_GPS =[71.3357 -156.3982; 71.3299 -156.4016; 71.3345 -156.4165; 71.3333 -156.4081];

for g = 1:length(crack_lat)
[xtemp, ytemp] = ll2xy(crack_lat(g),crack_long(g),geophone_GPS(1,1),geophone_GPS(1,2));
xcrack(g) = xtemp;
ycrack(g) = ytemp;
end
[xcrack,indc] = sort(xcrack);
ycrack = ycrack(indc);

bpFiltnode = designfilt('bandpassfir','FilterOrder',500, ...
         'CutoffFrequency1',12,'CutoffFrequency2',20, ...
         'SampleRate',FS);
%fvtool(bpFilt)

zdatafilt = [];
xdatafilt = [];
ydatafilt = [];

for i = 1:4
    zdatafilt(i,:) = filtfilt(bpFiltnode,zdata(:,i));
    xdatafilt(i,:) = filtfilt(bpFiltnode,xdata(:,i));
    ydatafilt(i,:) = filtfilt(bpFiltnode,ydata(:,i));
end

znoise = zeros(1,size(zdatafilt,1));
for zz = 1:size(zdatafilt,1)
    znoise(zz) = mean(abs(zdatafilt(zz,:)));
end
xnoise = zeros(1,size(xdatafilt,1));
for xx = 1:size(xdatafilt,1)
    xnoise(xx) = mean(abs(xdatafilt(xx,:)));
end
ynoise = zeros(1,size(ydatafilt,1));
for yy = 1:size(ydatafilt,1)
    ynoise(yy) = mean(abs(ydatafilt(yy,:)));
end

noise_mat = [znoise; xnoise; ynoise];

figure
for j = 1:4
subplot(4,1,j)
plot(t,zdatafilt(j,:))
hold on
plot(t,ydatafilt(j,:));
plot(t,xdatafilt(j,:));
title(['Geophone ' num2str(j)]);
legend('z-dir','x-dir','y-dir')
%xlabel('Time')
%ylabel('Amplitude')
%datetick('x','HH:MM:SS:FFF','keepticks');
%title(date);
set(gca,'fontsize',12);
grid on
%k = k+1;
end

%% tdoa localization estimate
c_range = 336;
N = 1000;
plotting = 1;
plotmap = 1;
calib_act = [0 0];

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

%try
    [loc_est,c_est,err,tdoa_mat] = loc_est_ambi(zdatafilt,xdatafilt,ydatafilt,xpos,ypos,1,length(zdatafilt),FS,c_range,N,plotting,calib_act,noise_mat,plotmap,1);
%catch
%    disp('Skipping this...')
%end

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
disp(err)

