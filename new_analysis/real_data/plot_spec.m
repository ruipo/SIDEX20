%% Plot Event time series
t1 = 58169;%1/26
t2 = 65651;%1/26
%t1 = 21000; %1/27_1 
%t2 = 25000; %1/27_1
%t1 = 14000; %1/27_2 
%t2 = 16000; %1/27_2
%t1 = 19600;%1/31
%t2 = 20400;%1/31
%t1 = 32000;%1/28_1
%t2 = 33500;%1/28_1
%t1 = 8500;%1/28_2
%t2 = 10000;%1/28_2

% Read in stand-alone geophone data
toffset = 6.2500e-05; % define time offset between wired and stand-alone geophone data (derived imperically)
t_trim = 33; % trim off data at the end of stand-alone geophone data (in seconds)
FS_node = 4096;
Table = readtable('/Volumes/RUIC_Backup/SIDEX20_data/SIDEx_node_data/transients_with_4_nodes_TD/2020-01-26_14_00_05.csv', 'HeaderLines',1,'Format','%D%f%f%f%f%f%f%f%f%f%f%f%f');
tn = table2array(Table(:,1));
tn = tn - toffset;
z_node = table2array(Table(:,[2 5 8 11]));
tempmat = ~isnan(z_node);
for i = 1:4
    tempvect1(i) = find(tempmat(:,i),1,'first');
    tempvect2(i) = find(tempmat(:,i),1,'last');
end
z_node = z_node(max(tempvect1):min(tempvect2)-t_trim*FS_node,:);
tn = tn(max(tempvect1):min(tempvect2)-t_trim*FS_node);
x_node = table2array(Table(:,[3 6 9 12]));
tempmat = ~isnan(x_node);
for i = 1:4
    tempvect1(i) = find(tempmat(:,i),1,'first');
    tempvect2(i) = find(tempmat(:,i),1,'last');
end
x_node = x_node(max(tempvect1):min(tempvect2)-t_trim*FS_node,:);
y_node = table2array(Table(:,[4 7 10 13]));
tempmat = ~isnan(y_node);
for i = 1:4
    tempvect1(i) = find(tempmat(:,i),1,'first');
    tempvect2(i) = find(tempmat(:,i),1,'last');
end
y_node = y_node(max(tempvect1):min(tempvect2)-t_trim*FS_node,:);

% interp wired geophone data to same sampling frequency as stand-alone geophone data
data_interp = interp1(t(t1:t2),data(t1:t2,:),tn,'spline');

figure
k=1;
for j = [1 4 7 10]
    subplot(8,1,k)
    plot(t(t1:t2),data(t1:t2,j))
    hold on
    plot(tn,data_interp(:,j))
    %xlabel('Time')
    %ylabel('Amplitude')
    %datetick('x','HH:MM:SS:FFF','keepticks');
    xlim([tn(1) tn(end)]);
    %title(date);
    set(gca,'fontsize',12);
    grid on
    k = k+1;
end

for j = 1:4
    subplot(8,1,k)
    plot(tn,z_node(:,j))
    xlim([tn(1) tn(end)]);
    set(gca,'fontsize',12);
    grid on
    k = k+1;
end
    
%% X-Z, Y-Z plot
figure
j = 1;
for i = [1 4 7 10]
subplot(2,2,j)
plot(data(t1:t2,i+1)./max(abs(data(t1:t2,i+1))),data(t1:t2,i)./max(abs(data(t1:t2,i))))
hold on
plot(data(t1:t2,i+2)./max(abs(data(t1:t2,i+2))),data(t1:t2,i)./max(abs(data(t1:t2,i))))
xlabel('X Axis')
ylabel('Z Axis')
grid on
set(gca,'fontsize',20)
xlim([-1 1])
ylim([-1 1])
title(['Geophone ' num2str(j)])
legend('X-Z','Y-Z')
j = j+1;
end
%% plot channels on same plot
figure
plot(t(t1:t2),data(t1:t2,[1 4 7 10]))
hold on
xlabel('Time')
ylabel('Amplitude')
title(date);
set(gca,'fontsize',20);
grid on

%% Spectrogram plot
chn = 1;
FS = 1000;
L = 128; % Window length (1ms)
R = L/2; % Overlap percentage
NFFT = L; 
time = 0:1/FS:size(data,1)/FS - 1/FS;

[S1,f,ts] = spectrogram(data(t1:t2,chn),hanning(L),R,NFFT,FS);
[S2,f,ts] = spectrogram(data(t1:t2,chn+3),hanning(L),R,NFFT,FS);
[S3,f,ts] = spectrogram(data(t1:t2,chn+6),hanning(L),R,NFFT,FS);
[S4,f,ts] = spectrogram(data(t1:t2,chn+9),hanning(L),R,NFFT,FS);

figure
subplot(4,1,1)
imagesc(ts,f,20*log10(abs(S1)))
set(gca,'YDir','normal')
colormap jet
caxis([-100 -30])
colorbar
%xlim([19 23])
ylim([0 100])
subplot(4,1,2)
imagesc(ts,f,20*log10(abs(S2)))
set(gca,'YDir','normal')
colormap jet
caxis([-100 -30])
colorbar
%xlim([19 23])
ylim([0 100])
subplot(4,1,3)
imagesc(ts,f,20*log10(abs(S3)))
set(gca,'YDir','normal')
colormap jet
caxis([-100 -30])
colorbar
%xlim([19 23])
ylim([0 100])
subplot(4,1,4)
imagesc(ts,f,20*log10(abs(S4)))
set(gca,'YDir','normal')
colormap jet
caxis([-100 -30])
colorbar
%xlim([19 23])
ylim([0 100])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% tdoa localization estimate set up
FS = 1000;
FS_node = 4096;
xpos = [0;-62.068;-3.310;-32.067];
ypos = [0;48.147;61.936;36.472];
xpos_node = [];
ypos_node = [];

geophone_GPS = [71.332947,-156.407381; 71.33338,-156.409125; 71.333504,-156.407474; 71.333275,-156.408282];
standalone_GPS =[71.3357 -156.3982; 71.3299 -156.4016; 71.3345 -156.4165; 71.3333 -156.4081];

for g = 1:4
[xtemp, ytemp] = ll2xy(standalone_GPS(g,1),standalone_GPS(g,2),geophone_GPS(1,1),geophone_GPS(1,2));
xpos(g+4) = xtemp;
ypos(g+4) = ytemp;
xpos_node(g) = xtemp;
ypos_node(g) = ytemp;
end

zdata = data_interp(:,[1 4 7 10]).';
xdata = data_interp(:,[2 5 8 11]).';
ydata = data_interp(:,[3 6 9 12]).';

bpFilt = designfilt('bandpassfir','FilterOrder',500, ...
         'CutoffFrequency1',8,'CutoffFrequency2',32, ...
         'SampleRate',FS);
bpFiltnode = designfilt('bandpassfir','FilterOrder',500, ...
         'CutoffFrequency1',4,'CutoffFrequency2',490, ...
         'SampleRate',FS_node);
%fvtool(bpFilt)

zdatafilt = [];
xdatafilt = [];
ydatafilt = [];

for i = 1:4
    zdatafilt(i,:) = filtfilt(bpFiltnode,zdata(i,:));
    xdatafilt(i,:) = filtfilt(bpFiltnode,xdata(i,:));
    ydatafilt(i,:) = filtfilt(bpFiltnode,ydata(i,:));
    zdatafilt(i+4,:) = filtfilt(bpFiltnode,z_node(:,i).');
    xdatafilt(i+4,:) = filtfilt(bpFiltnode,y_node(:,i).');
    ydatafilt(i+4,:) = filtfilt(bpFiltnode,x_node(:,i).');
end

figure
k=[1 3 5 7 2 4 6 8];
for j = 1:8
subplot(4,2,k(j))
plot(tn,zdatafilt(j,:))
hold on
plot(tn,ydatafilt(j,:));
plot(tn,xdatafilt(j,:));
title(['Geophone ' num2str(j)]);
if j<5
legend('z-dir','x-dir','y-dir')
else  
legend('z-dir','T-dir','L-dir')
end
%xlabel('Time')
%ylabel('Amplitude')
%datetick('x','HH:MM:SS:FFF','keepticks');
xlim([tn(1) tn(end)]);
%title(date);
set(gca,'fontsize',12);
grid on
%k = k+1;
end

% figure
% k=[1 3 5 7];
% for j = 1:4
% subplot(4,2,k(j))
% plot(tn,znodefilt(j,:))
% hold on
% plot(tn,ynodefilt(j,:));
% plot(tn,xnodefilt(j,:));
% title(['Geophone ' num2str(j)]);
% legend('z-dir','T-dir','L-dir')
% %xlabel('Time')
% %ylabel('Amplitude')
% %datetick('x','HH:MM:SS:FFF','keepticks');
% xlim([tn(1) tn(end)]);
% %title(date);
% set(gca,'fontsize',12);
% grid on
% %k = k+1;
% end

% xpos(5) = []; xpos(6) = [];
% ypos(5) = []; ypos(6) = [];
% zdatafilt(5,:) = [];
% zdatafilt(6,:) = [];

% Use phase to estimate direction of arrival
psivec_filt = [];
figure
for j = 1:8
thetavec = 0:1:360;
xz = (xdatafilt(j,:).*zdatafilt(j,:));
%xz(xz<0) = 0;
yz = (ydatafilt(j,:).*zdatafilt(j,:));
%yz(yz<0) = 0;
xzp = (xdatafilt(j,:).*gradient(zdatafilt(j,:)));
%xzp(xzp>0) = 0;
yzp = (ydatafilt(j,:).*gradient(zdatafilt(j,:)));
%yzp(yzp>0) = 0;

psivec_z = cos(deg2rad(thetavec)).^2*sum(xz.^2) + sin(deg2rad(thetavec)).^2*sum(yz.^2) + 2*cos(deg2rad(thetavec)).*sin(deg2rad(thetavec))*sum(xz.*yz);
psivec_zp = cos(deg2rad(thetavec)).^2*sum(xzp.^2) + sin(deg2rad(thetavec)).^2*sum(yzp.^2) + 2*cos(deg2rad(thetavec)).*sin(deg2rad(thetavec))*sum(xzp.*yzp);
psivec_filt(j,:) = psivec_z/min(psivec_z) + psivec_zp/min(psivec_zp);
[~, locs] = findpeaks(psivec_filt(j,:));
if length(locs)<2
    locs = [locs 1];
end
r1 = xdatafilt(j,:)*cos(deg2rad(thetavec(locs(1)))) + ydatafilt(j,:)*sin(deg2rad(thetavec(locs(1))));
r2 = xdatafilt(j,:)*cos(deg2rad(thetavec(locs(2)))) + ydatafilt(j,:)*sin(deg2rad(thetavec(locs(2))));

sigma1 = 0;
for i = 2:length(r1)-1
    sigma1 = sigma1 + zdatafilt(j,i)*(r1(i+1)-r1(i-1));
end
    
sigma2 = 0;
for i = 2:length(r1)-1
    sigma2 = sigma2 + zdatafilt(j,i)*(r2(i+1)-r2(i-1));
end
disp([sigma1 sigma2])
[~,ind] = max([sigma1 sigma2]);

subplot(4,2,j);
polarplot(deg2rad(thetavec), 10*log10(psivec_filt(j,:)))
hold on
polarplot(deg2rad(thetavec(locs(ind)))*ones(10,1), linspace(0,15,10))
end




%% tdoa localization estimate
c_range = 384;
N = 1000;
plotting = 1;

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
    [loc_est,c_est,err,tdoa_mat] = loc_est_hyp(zdatafilt,xpos,ypos,1,length(zdatafilt),FS_node,c_range,N,plotting);
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
disp(err)
