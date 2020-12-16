%% Data Band-pass filter
f1 = 12;
f2 = 20;
FS = 1000;

bandfilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',f1,'CutoffFrequency2',f2,'SampleRate',FS);
%data_filt = filtfilt(bandfilt,data);
data_cor = data;
data_cor(:,[2 5 8 11]) = data(:,[2 5 8 11]).*cos(3*pi/2) - data(:,[3 6 9 12]).*sin(3*pi/2);
data_cor(:,[3 6 9 12]) = data(:,[2 5 8 11]).*sin(3*pi/2) + data(:,[3 6 9 12]).*cos(3*pi/2);
data_filt = filtfilt(bandfilt,data_cor);

znoise = zeros(1,4);
index = 1;
for zz = [1 4 7 10]
    znoise(index) = mean(abs(data_filt(:,zz)));
    index = index + 1;
end
xnoise = zeros(1,4);
index = 1;
for xx = [2 5 8 11]
    xnoise(index) = mean(abs(data_filt(:,xx)));
    index = index + 1;
end
ynoise = zeros(1,4);
index = 1;
for yy = [3 6 9 12]
    ynoise(index) = mean(abs(data_filt(:,yy)));
    index = index + 1;
end

event_inds = xlsread('event_inds.xlsx',1,'B2:C189');

figure
plot(t,data_filt(:,1))
hold on
plot(t,data_filt(:,4))
plot(t,data_filt(:,7))
plot(t,data_filt(:,10))
plot(epochtime_event,zeros(length(epochtime_event),1),'*')

%% plot event location and moveout
i = 9;

event = data_cor(event_inds(i,1):event_inds(i,2),:);
zevent = event(:,[1 4 7 10]);
xevent = event(:,[2 5 8 11]);
yevent = event(:,[3 6 9 12]);

zSNR = 10*log10(max(abs(zevent))./znoise);
xSNR = 10*log10(max(abs(xevent))./xnoise);
ySNR = 10*log10(max(abs(yevent))./ynoise);
SNR_mat = [zSNR; xSNR; ySNR];

geophone_loc = [24.36 -36.64;-37.71 11.51;21.05 25.30;-7.71 -0.167];
[eventx, eventy]= ll2xy(lat_event(i),long_event(i),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));

deg = rad2deg(atan2(eventy,eventx));
if deg<0
    ang_act = 180+deg;
else
    ang_act = deg;
end
revent = xevent.*cos(deg2rad(ang_act)) - yevent.*sin(deg2rad(ang_act));

taxis = 0:1/FS:size(zevent,1)/FS-(1/FS);

figure
plot(geophone_loc(:,1),geophone_loc(:,2),'ro');
hold on
plot(eventx,eventy ,'k*')
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
set(gca,'fontsize',20)

figure
subplot(4,2,1)
plot(taxis, zevent(:,1));
title('Geophone 1 Vertical Axis')
set(gca,'fontsize',20)
ylabel('R = 49.76 m')
ylim([-0.0025 0.0025])
grid on
subplot(4,2,3)
plot(taxis, zevent(:,3));
title('Geophone 3 Vertical Axis')
set(gca,'fontsize',20)
ylabel('R = 48.76 m')
ylim([-0.0025 0.0025])
grid on
subplot(4,2,5)
plot(taxis, zevent(:,2));
title('Geophone 2 Vertical Axis')
set(gca,'fontsize',20)
ylabel('R = 29.92 m')
ylim([-0.0025 0.0025])
grid on
subplot(4,2,7)
plot(taxis, zevent(:,4));
title('Geophone 4 Vertical Axis')
set(gca,'fontsize',20)
ylabel('R = 11.39 m')
ylim([-0.0025 0.0025])
grid on
subplot(4,2,2)
plot(taxis, revent(:,1));
title('Geophone 1 Radial Axis')
set(gca,'fontsize',20)
ylabel('R = 49.76 m')
ylim([-0.0025 0.0025])
grid on
subplot(4,2,4)
plot(taxis, revent(:,3));
title('Geophone 3 Radial Axis')
set(gca,'fontsize',20)
ylabel('R = 48.76 m')
ylim([-0.0025 0.0025])
grid on
subplot(4,2,6)
plot(taxis, revent(:,2));
title('Geophone 2 Radial Axis')
set(gca,'fontsize',20)
ylabel('R = 29.92 m')
ylim([-0.0025 0.0025])
grid on
subplot(4,2,8)
plot(taxis, revent(:,4));
title('Geophone 4 Radial Axis')
set(gca,'fontsize',20)
ylabel('R = 11.39 m')
ylim([-0.0025 0.0025])
grid on

figure
for i = 1:4
    subplot(2,2,i)
    plot(movmean(zevent(500:end,i),50)./max(abs(movmean(zevent(500:end,i),50))),movmean(revent(500:end,i),50)./max(abs(movmean(revent(500:end,i),50))))
    xlabel('H-axis')
    ylabel('Z-axis')
    xlim([-1 1])
    ylim([-1 1])
    grid on
    axis square
    title(['Geophone ' num2str(i)])
    set(gca,'fontsize',20)
end


%% Spectrogram View
figure
win_len = 128;

plotting = 0;
Sz_vec_mat = zeros(length(event_inds),65);
Sx_vec_mat = zeros(length(event_inds),65);
Sy_vec_mat = zeros(length(event_inds),65);

for ee = 1:length(event_inds)
    ee
    
    if ~isnan(event_inds(ee,1))
    
    d1 = lldistkm([lat_event(ee) long_event(ee)],[geophone_GPS(1,1) geophone_GPS(1,2)]);
    d2 = lldistkm([lat_event(ee) long_event(ee)],[geophone_GPS(2,1) geophone_GPS(2,2)]);
    d3 = lldistkm([lat_event(ee) long_event(ee)],[geophone_GPS(3,1) geophone_GPS(3,2)]);
    d4 = lldistkm([lat_event(ee) long_event(ee)],[geophone_GPS(4,1) geophone_GPS(4,2)]);
    
    event = data(event_inds(ee,1):event_inds(ee,2)+500,:);
    event_time = t(event_inds(ee,1):event_inds(ee,2)+500);
    source_time = epochtime_event(ee);
    t_start = event_time(1)-source_time;
    
    eventz = event(:,[1 4 7 10]);
    eventx = event(:,[2 5 8 11]);
    eventy = event(:,[3 6 9 12]);
    
    [Sz,~,~] = spectroView(eventz,FS, win_len);
    [Sx,~,~] = spectroView(eventx,FS, win_len);
    [Sy,fvec,tvec] = spectroView(eventy,FS, win_len);
    
    if plotting
        specSubplot(Sz(:,:,1),tvec+t_start,fvec,4,3,1)
        title('Z-axis Spectrogram')
        ylabel({['R = ' num2str(round(d1*1000,2)) ' m'],'Frequeny (Hz)'})
        specSubplot(Sz(:,:,2),tvec+t_start,fvec,4,3,4)
        ylabel({['R = ' num2str(round(d2*1000,2)) ' m'],'Frequeny (Hz)'})
        specSubplot(Sz(:,:,3),tvec+t_start,fvec,4,3,7)
        ylabel({['R = ' num2str(round(d3*1000,2)) ' m'],'Frequeny (Hz)'})
        specSubplot(Sz(:,:,4),tvec+t_start,fvec,4,3,10)
        ylabel({['R = ' num2str(round(d4*1000,2)) ' m'],'Frequeny (Hz)'})
        specSubplot(Sx(:,:,1),tvec+t_start,fvec,4,3,2)
        title('X-axis Spectrogram')
        specSubplot(Sx(:,:,2),tvec+t_start,fvec,4,3,5)
        specSubplot(Sx(:,:,3),tvec+t_start,fvec,4,3,8)
        specSubplot(Sx(:,:,4),tvec+t_start,fvec,4,3,11)
        xlabel('Time after event (sec)')
        specSubplot(Sy(:,:,1),tvec+t_start,fvec,4,3,3)
        title('Y-axis Spectrogram')
        specSubplot(Sy(:,:,2),tvec+t_start,fvec,4,3,6)
        specSubplot(Sy(:,:,3),tvec+t_start,fvec,4,3,9)
        specSubplot(Sy(:,:,4),tvec+t_start,fvec,4,3,12)
        
        pause
        clf
    end
    
    Sz_vec_mat(ee,:) = mean(mean(Sz,3),2);
    Sx_vec_mat(ee,:) = mean(mean(Sx,3),2);
    Sy_vec_mat(ee,:) = mean(mean(Sy,3),2);
    
    else
    Sz_vec_mat(ee,:) = NaN(65,1);
    Sx_vec_mat(ee,:) = NaN(65,1);
    Sy_vec_mat(ee,:) = NaN(65,1);
    
    end
    
end

figure
subplot(3,1,1)
plot(fvec, 20*log10(abs(Sz_vec_mat)./1E-6),'.');
hold on
plot(fvec, nanmean(20*log10(abs(Sz_vec_mat)./1E-6),1),'k','linewidth',2);
title('Z-axis')
set(gca,'fontsize',20)
ylim([0 60])
subplot(3,1,2)
plot(fvec, 20*log10(abs(Sx_vec_mat)./1E-6),'.');
hold on
plot(fvec, nanmean(20*log10(abs(Sx_vec_mat)./1E-6),1),'k','linewidth',2);
ylabel('dB re 1 uPa')
title('X-axis')
set(gca,'fontsize',20)
ylim([0 60])
subplot(3,1,3)
plot(fvec, 20*log10(abs(Sy_vec_mat)./1E-6),'.');
hold on
plot(fvec, nanmean(20*log10(abs(Sy_vec_mat)./1E-6),1),'k','linewidth',2);
xlabel('Frequency (Hz)')
title('Y-axis')
set(gca,'fontsize',20)
ylim([0 60])

%% Group velocity estimate with toda and multiple filter method

fcs = 2:2:50;
vg_mat = zeros(length(event_inds),3,length(fcs),6);
bw = 10;
d1 = lldistkm([lat_event(1) long_event(1)],[geophone_GPS(1,1) geophone_GPS(1,2)]);
d2 = lldistkm([lat_event(2) long_event(2)],[geophone_GPS(2,1) geophone_GPS(2,2)]);
d3 = lldistkm([lat_event(3) long_event(3)],[geophone_GPS(3,1) geophone_GPS(3,2)]);
d4 = lldistkm([lat_event(4) long_event(4)],[geophone_GPS(4,1) geophone_GPS(4,2)]);
range_rec = [d1 d2 d3 d4]*1000;

for ee = 1:length(event_inds)
    ee
    if ~isnan(event_inds(ee,1))
        event = data(event_inds(ee,1):event_inds(ee,2)+500,:);
        event_time = t(event_inds(ee,1):event_inds(ee,2)+500);

        eventz = event(:,[1 4 7 10]);
        eventx = event(:,[2 5 8 11]);
        eventy = event(:,[3 6 9 12]);

        try   
            [vg_estz] = tdoaGroupVelocityEst(eventz, event_time, fcs, bw, range_rec, FS);
            [vg_estx] = tdoaGroupVelocityEst(eventx, event_time, fcs, bw, range_rec, FS);
            [vg_esty] = tdoaGroupVelocityEst(eventy, event_time, fcs, bw, range_rec, FS);

            vg_mat(ee,1,:,:) = vg_estz;
            vg_mat(ee,2,:,:) = vg_estx;
            vg_mat(ee,3,:,:) = vg_esty;

        catch
            disp(['Skipping ' num2str(ee)])
        end
    
    end
    
end

% Z-axis
vgz_mat = squeeze(vg_mat(:,1,:,:));
vgz_mat(vgz_mat <= 10) = NaN;
vgz_mat(vgz_mat == Inf) = NaN;
vgz_mat(vgz_mat >= 3000) = NaN;
figure
hold on
for ee = 1:length(event_inds)
    %plot(fcs,squeeze(vgz_mat(ee,:,:)),'k*')
    plot(fcs,nanmean(squeeze(vgz_mat(ee,:,:)),2),'k*')
    grid on
    
end

plot(fcs,nanmean(nanmean(vgz_mat,3),1))

% X-axis
vgx_mat = squeeze(vg_mat(:,2,:,:));
vgx_mat(vgx_mat <= 10) = NaN;
vgx_mat(vgx_mat == Inf) = NaN;
vgx_mat(vgx_mat >= 3000) = NaN;
figure
hold on
for ee = 1:length(event_inds)
    %plot(fcs,squeeze(vgz_mat(ee,:,:)),'k*')
    plot(fcs,nanmean(squeeze(vgx_mat(ee,:,:)),2),'k*')
    grid on
    
end

plot(fcs,nanmean(nanmean(vgx_mat,3),1))

% Y-axis
vgy_mat = squeeze(vg_mat(:,2,:,:));
vgy_mat(vgy_mat <= 10) = NaN;
vgy_mat(vgy_mat == Inf) = NaN;
vgy_mat(vgy_mat >= 3000) = NaN;
figure
hold on
for ee = 1:length(event_inds)
    %plot(fcs,squeeze(vgz_mat(ee,:,:)),'k*')
    plot(fcs,nanmean(squeeze(vgy_mat(ee,:,:)),2),'k*')
    grid on
    
end

plot(fcs,nanmean(nanmean(vgy_mat,3),1))

% All axes
vg_est_tot = [nanmean(nanmean(vgz_mat,3)); nanmean(nanmean(vgx_mat,3)); nanmean(nanmean(vgy_mat,3))];
figure;
plot(fcs,mean(vg_est_tot,1))
grid on
xlabel('Frequency (Hz)')
ylabel('Group Velocity (m/s)')

%% Motion Product Detector

%thetas = [275 275 300 311];
%thetas = [270 270 270 270];

figure

for ee = 1:188
    if ~isnan(event_inds(ee,1))
        [event1x, event1y]= ll2xy(lat_event(ee),long_event(ee),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
        geophone_loc = [24.36 -36.64;-37.71 11.51;21.05 25.30;-7.71 -0.167];
        event = data_filt(event_inds(ee,1):event_inds(ee,2)+1000,:);
        plot(geophone_loc(:,1),geophone_loc(:,2),'ro')
        hold on
        plot(event1x, event1y, 'b*')
        
        zSNR = 10*log10(max(abs(event(:,[1 4 7 10])))./znoise);
        xSNR = 10*log10(max(abs(event(:,[2 5 8 11])))./xnoise);
        ySNR = 10*log10(max(abs(event(:,[3 6 9 12])))./ynoise);
        SNR_mat = [zSNR; xSNR; ySNR];

        for chn = 1:4
            if min(SNR_mat(:,chn)) > 10
                chn
                datain = event(1:end,chn+(chn-1)*2:chn+(chn-1)*2+2);
                [X_HiV, Y_HiV, X_est, Y_est] = MPD(datain,200);
%                 plot(movmean(datain,50))
%                 pause
%                 clf
                plot(X_est+geophone_loc(chn,1),Y_est+geophone_loc(chn,2),'k');
            else
                datain = event(1:end,chn+(chn-1)*2:chn+(chn-1)*2+2);
                [X_HiV, Y_HiV, X_est, Y_est] = MPD(datain,200);
%                 plot(movmean(datain,50))
%                 pause
%                 clf
                plot(X_est+geophone_loc(chn,1),Y_est+geophone_loc(chn,2),'r');
            end
        end
        grid on
        title(['Calibration Event ' num2str(ee)]);
        xlim([-150 150])
        ylim([-150 150])
        xlabel('X (m)')
        ylabel('Y (m)')
        set(gca, 'fontsize',20)
        axis square
        pause
        clf;
    end
end

%% plot V-H

tdoa_xy = zeros(4,188);
figure

for ee = 1:188
    if ~isnan(event_inds(ee,1))
        event = data(event_inds(ee,1):event_inds(ee,2),:);
        
        p = 1;
        for chn = 1:4
            datain = event(:,chn+(chn-1)*2:chn+(chn-1)*2+2);
            [corr,lags] = xcorr(abs(hilbert(datain(:,2))),abs(hilbert(datain(:,3))));
            [~,loc] = max(corr);
            tdoa_xy(chn,ee) = lags(loc)/FS;
        
            subplot(2,4,p)
            plot(datain(:,2)./max(abs(datain(:,2))),datain(:,1)./max(abs(datain(:,1))))
            xlabel('X-axis')
            ylabel('Z-axis')
            xlim([-1 1])
            ylim([-1 1])
            grid on
            axis square
            title(['Geophone ' num2str(chn)])
            set(gca,'fontsize',20)
            p = p+1;
            subplot(2,4,p)
            plot(datain(:,3)./max(abs(datain(:,3))),datain(:,1)./max(abs(datain(:,1))))
            xlabel('Y-axis')
            ylabel('Z-axis')
            xlim([-1 1])
            ylim([-1 1])
            grid on
            axis square
            title(['Geophone ' num2str(chn)])
            set(gca,'fontsize',20)
            p = p+1;

        end
        pause
        clf
    end
end
%% Event Beamform

zevent = data(:,[1 4 7 10]);
geophone_loc = [24.36 -36.64 0;-37.71 11.51 0;21.05 25.30 0;-7.71 -0.167 0];

FS = 1000;
elev = 0;
az = 0:1:360;
c = 350;
f_range = [16 32];
NFFT = 16;
overlap = 0;
weighting = 'uniform';


ang_est_mat = zeros(1,length(event_inds));

for i = 1:188
    i 
    for c= 350
    ss = event_inds(i,1);
    es = event_inds(i,2);
    if isnan(ss) || isnan(es)
        ang_est_mat(i) = NaN;
    else
        figure(1)
        plot(geophone_loc(:,1),geophone_loc(:,2),'ro');
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
        [eventx, eventy]= ll2xy(lat_event(i),long_event(i),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
        plot(eventx,eventy,'*')
        deg = rad2deg(atan2(eventy,eventx));
        if deg<0
            ang_act = 360+deg;
        else
            ang_act = deg;
        end

        %try
            window = hanning(es-ss+1);
            [beamform_output,~,~] = beamform_3D(zevent(ss:es,:),geophone_loc,FS,elev,az,c,f_range,NFFT,window,overlap,weighting);
            ang_list = squeeze(mean(mean(beamform_output,1),4));
            [~,ind] = max(abs(ang_list));
            ang_est = az(ind);
            if (ang_est>90 && ang_est<=270)
                xlist = -500:0;
                ylist = xlist*tan(deg2rad(ang_est));
                plot(xlist,ylist,'b')
            else
                xlist = 0:500;
                ylist = xlist*tan(deg2rad(ang_est));
                plot(xlist,ylist,'b')
            end
        %catch
            %disp('Skipping this...')
        %end

        ang_est_mat(:,i) = ang_est;

        disp(['ang estimate = ', num2str(ang_est), '.']);
        disp(['ang actual = ', num2str(ang_act), '.']);

    end
    
    pause;
    clf;
    end
   
end

%% phase event direction of arrival estimate

zdatafilt = data_filt(:,[1 4 7 10]).';
xdatafilt = data_filt(:,[2 5 8 11]).';
ydatafilt = data_filt(:,[3 6 9 12]).'; 

for i = 1:188
    i 
    ss = event_inds(i,1);
    es = event_inds(i,2);
    if isnan(ss) || isnan(es)
        ang_est_mat(i) = NaN;
    else
        figure(1)
        plot(geophone_loc(:,1),geophone_loc(:,2),'ro');
        xlabel('X position (m)')
        xlabel('Y position (m)')
        xlim([-350 350]);
        ylim([-350 350]);
        xL = xlim;
        yL = ylim;
        line([0 0], yL,'color','black');
        line(xL, [0 0],'color','black');
        grid on
        hold on
        [eventx, eventy]= ll2xy(lat_event(i),long_event(i),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
        plot(eventx,eventy,'*')
        
        psivec_filt = [];
        theta_est = [];
%         figure(2)
        for j = 1:4
            thetavec = 0:1:360;
            
            xz = (xdatafilt(j,ss:es).*zdatafilt(j,ss:es));
            xz(xz<0) = 0;
            yz = (ydatafilt(j,ss:es).*zdatafilt(j,ss:es));
            yz(yz<0) = 0;
            xzp = (xdatafilt(j,ss:es).*gradient(zdatafilt(j,ss:es)));
            xzp(xzp>0) = 0;
            yzp = (ydatafilt(j,ss:es).*gradient(zdatafilt(j,ss:es)));
            yzp(yzp>0) = 0;
            
            psivec_z = cos(deg2rad(thetavec)).^2*sum(xz.^2) + sin(deg2rad(thetavec)).^2*sum(yz.^2) + 2*cos(deg2rad(thetavec)).*sin(deg2rad(thetavec))*sum(xz.*yz);
            psivec_zp = cos(deg2rad(thetavec)).^2*sum(xzp.^2) + sin(deg2rad(thetavec)).^2*sum(yzp.^2) + 2*cos(deg2rad(thetavec)).*sin(deg2rad(thetavec))*sum(xzp.*yzp);
            psivec_filt(j,:) = psivec_z/min(psivec_z) + psivec_zp/min(psivec_zp);
            [~, locs] = findpeaks(psivec_filt(j,:));
            if length(locs)>1
                length(locs)
                r1 = xdatafilt(j,:)*cos(deg2rad(thetavec(locs(1)))) + ydatafilt(j,:)*sin(deg2rad(thetavec(locs(1))));
                r2 = xdatafilt(j,:)*cos(deg2rad(thetavec(locs(2)))) + ydatafilt(j,:)*sin(deg2rad(thetavec(locs(2))));

                sigma1 = 0;
                for k = 2:length(r1)-1
                    sigma1 = sigma1 + zdatafilt(j,k)*(r1(k+1)-r1(k-1));
                end

                sigma2 = 0;
                for k = 2:length(r1)-1
                    sigma2 = sigma2 + zdatafilt(j,k)*(r2(k+1)-r2(k-1));
                end

                [~,ind] = max([sigma1 sigma2]);
                disp([sigma1 sigma2])
                theta_est(j) = rad2deg(deg2rad(thetavec(locs(ind)))+3*pi/2);
            else
                theta_est(j) = rad2deg(deg2rad(thetavec(locs(1)))+3*pi/2);
            end
                
            
            figure(1)
            if (theta_est(j)>90 && theta_est(j)<=270)
                xlist = -500:0;
                ylist = xlist*tan(deg2rad(theta_est(j)));
                plot(xlist+geophone_loc(j,1),ylist+geophone_loc(j,2),'b')
            else
                xlist = 0:500;
                ylist = xlist*tan(deg2rad(theta_est(j)));
                plot(xlist+geophone_loc(j,1),ylist+geophone_loc(j,2),'b')
            end
            
%             figure(2)
%             subplot(2,2,j);
%             polarplot(deg2rad(thetavec), 10*log10(psivec_filt(j,:)))
%             hold on
%             polarplot(deg2rad(thetavec(locs(ind)))*ones(10,1), linspace(0,15,10))
        end
        
    end
    pause;
    clf;
end


%% Event Locatization

zevent = data_filt(:,[1 4 7 10]).';
xevent = data_filt(:,[2 5 8 11]).';
yevent = data_filt(:,[3 6 9 12]).'; 
zeventuf = data_cor(:,[1 4 7 10]).';
xeventuf = data_cor(:,[2 5 8 11]).';
yeventuf = data_cor(:,[3 6 9 12]).'; 

%plot xyz time series for event i
i = 1;
[eventx, eventy]= ll2xy(lat_event(i),long_event(i),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
t1 = event_inds(i,1);
t2 = event_inds(i,2);

figure
k=[1 3 5 7];
for j = 1:4
subplot(4,1,j)
plot(t(t1:t2),zevent(j,t1:t2))
hold on
plot(t(t1:t2),yevent(j,t1:t2));
plot(t(t1:t2),xevent(j,t1:t2));
title(['Geophone ' num2str(j) '; Event Location: ' num2str(eventx) ', ' num2str(eventy)]);
legend('z-dir','x-dir','y-dir')
%xlabel('Time')
%ylabel('Amplitude')
%datetick('x','HH:MM:SS:FFF','keepticks');
xlim([t(t1) t(t2)]);
ylim([-0.002 0.002])
%title(date);
set(gca,'fontsize',12);
grid on
%k = k+1;
end
%%
geophone_loc = [24.36 -36.64;-37.71 11.51;21.05 25.30;-7.71 -0.167];

FS = 1000;
c_range = 242;
N = 5000;
plotting = 1;
noise_mat = [znoise; xnoise; ynoise];

figure
plot(geophone_loc(:,1),geophone_loc(:,2),'ro');
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

pos_est_mat = zeros(2,length(event_inds));
c_est_mat = zeros(length(event_inds),1);
error_mat = zeros(length(event_inds),1);

for i = 166
    i 
    
    ss = event_inds(i,1);
    es = event_inds(i,2);
    if isnan(ss) || isnan(es)
        pos_est_mat(:,i) = [NaN; NaN];
        c_est_mat(i) = NaN;
        error_mat(i) = NaN;
    else
        [eventx, eventy]= ll2xy(lat_event(i),long_event(i),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
        calib_act = [eventx eventy];

        %try
            [loc_est,c_est,err,~] = loc_est_calib_testfn(zevent,xevent, yevent,zeventuf,xeventuf, yeventuf,geophone_loc(:,1),geophone_loc(:,2),ss,es,FS,c_range,N,plotting,calib_act,noise_mat);
        %catch
        %    disp('Skipping this...')
        %end

        pos_est_mat(:,i) = loc_est;
        c_est_mat(i) = c_est;
        error_mat(i) = err;

        disp(['x estimate = ', num2str(loc_est(1)), '; y estimate = ',num2str(loc_est(2)) '.']);
        disp(['propagation speed estimate = ', num2str(c_est), 'm/s.']);
        disp(['estimate error = ', num2str(err), '.']);
    %     
    %     if c_est < 250
    %         color = 'r';
    %     elseif (250 <= c_est) && (c_est < 500)
    %         color = 'b';
    %     elseif (500 <= c_est) && (c_est < 750)
    %         color = 'g';
    %     else
    %         color = 'm';
    %     end
    %     
    %     figure(1);
    %     plot(loc_est(1),loc_est(2),[color '*']);
    end
   
end
%%
xlist = -350:1:350;
ylist = -350:1:350;
clist = 350;
err_mat = tdoa_cmp(xlist, ylist, clist, geophone_loc(:,1), geophone_loc(:,2), tdoa_mat);
[xind,yind,cind] = ind2sub(size(err_mat),find(err_mat == min(min(min(err_mat)))));

disp([xlist(xind) ylist(yind) clist(cind)])

figure
h = pcolor(xlist,ylist,err_mat(:,:,cind));
set(h,'Edgecolor','None')

%% Location Estimation plot
path = '/Users/Rui/Documents/Graduate/Research/SIDEX/SIDEX20/new_analysis/'; % set path to where you want to save figures; create folder "calib_results" at this path location;

figure('units','normalized','outerposition',[0 0 1 1])
plot(geophone_loc(:,1),geophone_loc(:,2),'ro', 'MarkerFaceColor', 'r');
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
h(3) = plot(NaN,NaN,'b*','MarkerSize',6);
%h(4) = plot(NaN,NaN,'b*','MarkerSize',6);
%h(5) = plot(NaN,NaN,'m*','MarkerSize',6);
%h(6) = plot(NaN,NaN,'g*','MarkerSize',6);
legend(h, 'Geophone Locations','True Location','Estimated Location');% (v < 250m/s)','Estimated Location (250m/s <= v < 500m/s)','Estimated Location (500m/s <= v < 750m/s)','Estimated Location (v >= 750m/s)');
set(0,'DefaultLegendAutoUpdate','off')

for i = 1:188
    i
    [eventx, eventy]= ll2xy(lat_event(i),long_event(i),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
    calib_act = [eventx eventy];
    
    c_est = c_est_mat(i);
    if c_est < 200
        color = 'r';
    elseif (200 <= c_est) && (c_est < 500)
        color = 'b';
    elseif (500 <= c_est) && (c_est < 750)
        color = 'm';
    else
        color = 'g';
    end
    
    title([datestr(datetime(epochtime_event(i), 'convertfrom','posixtime')) ' UTC'])

    plot(eventx,eventy,'ks','MarkerSize',8,'linewidth',1.25);
    plot(pos_est_mat(1,i),pos_est_mat(2,i),[color '*'],'MarkerSize',6,'linewidth',1.25);
    quiver(pos_est_mat(1,i),pos_est_mat(2,i),eventx-pos_est_mat(1,i),eventy-pos_est_mat(2,i),0,'color',[1 0.5 1]);
    saveas(gcf,[path,'loc_est_figs2/' num2str(i) '.png']);
    pause(0.1)
end

%% Location Estimation Results movie

testf = dir([path 'loc_est_figs2/*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([path 'loc_est_figs2/loc_est.avi']);
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
  
  
%% Location Estimation Error plot
distfc = zeros(188,1);
for i = 1:188
    [eventx, eventy]= ll2xy(lat_event(i),long_event(i),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
    distfc(i) = sqrt((eventx).^2 + (eventy).^2);
end
figure
semilogy(distfc,error_mat,'b.','MarkerSize',10)
grid on
set(gca,'fontsize',20)
hold on
% f=fit(distfc(~isnan(error_mat)),error_mat(~isnan(error_mat)),'exp2');
% plot(f)
xlabel('Distance from Center Geophone (m)')
ylabel('Error in Location Estimation (m)')
title('Estimation Error vs. Distance from Center Geophone')

%% Group Velocity Estimation

FS = 1000;
f1 = 16;
f2 = 32;
vgs_est = zeros(length(event_inds),1);

for ee = 1:length(event_inds)
    ee
    if ~isnan(event_inds(ee,1))
        vg_est = zeros(12,1);
        event = data(event_inds(ee,1):event_inds(ee,2),:);
        event_time = t(event_inds(ee,1):event_inds(ee,2));
        for chn = 1:12
            datain = event(:,chn);
            source_time = epochtime_event(ee);

            if chn <= 3
                dist = calib_dist(1,1);
            elseif 3 < chn <= 6
                dist = calib_dist(2,1);
            elseif 6 < chn <= 9
                dist = calib_dist(3,1);
            else
                dist = calib_dist(4,1);
            end

            vg_est(chn) = groupVelocityEst(datain, event_time, f1, f2, FS,source_time, dist);
        end
        vgs_est(ee) = mean(vg_est);
    else
        vgs_est(ee) = NaN;
    end
end
% vgs_est(159) = NaN;
% vgs_est(vgs_est>100) = NaN;
n = 4; % average every n values
vgs_est_mean = arrayfun(@(i) mean(vgs_est(i:i+n-1)),1:n:length(vgs_est)-n+1)'; 

figure
plot3(long_event(1:4:end),lat_event(1:4:end),vgs_est_mean,'*')
grid on
xlabel('Longitude')
ylabel('Latitude')
zlabel('Group Speed (m/s)')
xlim([-156.41 -156.398]);
ylim([71.332 71.337]);
zlim([0 80]);
title('Group Velocity Estimation (Freq = 16-32Hz)')
set(gca,'fontsize',25);

% cf = fit([long_event(~isnan(vgs_est_mean)) lat_event(~isnan(vgs_est_mean))],vgs_est(~isnan(vgs_est_mean)),'lowess');
% % figure
% % plot([long_event(~isnan(vgs_est)) lat_event(~isnan(vgs_est))],vgs_est(~isnan(vgs_est)))
% 
% [cx,cy] = meshgrid(-156.4096:0.0001:-156.3998, 71.3327:0.00001:71.3360);
% c_fit = zeros(size(cx));
% for i = 1:size(cx,1)
%     c_fit(i,:) = cf([cx(i,:).' cy(i,:).']);
% end
% % 
% %c_fit(c_fit>50) = NaN;
% c_fit(c_fit<0) = NaN;
% figure
% h = pcolor(cx(1,:),cy(:,1),c_fit);
% % caxis([50 1000])
% set(h,'Edgecolor','None')
% hold on
% plot(geophone_GPS(:,2),geophone_GPS(:,1),'ro', 'MarkerFaceColor', 'r');
% xlabel('Longitude')
% ylabel('Latitude')
% xlim([-156.4096 -156.3998]);
% ylim([71.3327 71.3360]);
% xL = xlim;
% yL = ylim;
% %plot(long_event(~isnan(vgs_est)),lat_event(~isnan(vgs_est)),'k*')
% title('Group Velocity Estimation (Freq = 16-32Hz)')
% colormap(flipud(jet))
% colorbar
% set(gca,'fontsize',20)


%% 
diff = event_inds(:,2) - event_inds(:,1);
pad = (3000-diff)/2;
zevent = data(:,[1 4 7 10]).';
xevent = data(:,[2 5 8 11]).';
yevent = data(:,[3 6 9 12]).'; 
zSTFT_mat = [];
xSTFT_mat = [];
ySTFT_mat = [];

nn = 1;
%figure
for ee = 1:length(event_inds)
    ee
    if ~isnan(event_inds(ee,1))
        for chn = 1:4
            segment = zevent(chn,event_inds(ee,1)-pad(ee):event_inds(ee,2)+pad(ee));
            [STFT, f, t] = stft(segment,hanning(64),32,64,1000);
            [~,yind] = find(STFT==max(max(STFT)));
            STFT = circshift(STFT,size(STFT,2)/2-yind,2);
            zSTFT_mat(:,:,nn) = STFT;
            nn = nn+1;
%             pcolor(t,f,20*log10(abs(STFT)))
%             title(['event #' num2str(ee) ', channel = ' num2str(chn)])
%             colorbar
%             caxis([-80 -40])
%             pause
%             clf;
        end
    end
end

nn = 1;
%figure
for ee = 1:length(event_inds)
    if ~isnan(event_inds(ee,1))
        for chn = 1:4
            segment = xevent(chn,event_inds(ee,1)-pad(ee):event_inds(ee,2)+pad(ee));
            [STFT, f, t] = stft(segment,hanning(64),32,64,1000);
            [~,yind] = find(STFT==max(max(STFT)));
            STFT = circshift(STFT,size(STFT,2)/2-yind,2);
            xSTFT_mat(:,:,nn) = STFT;
            nn = nn+1;
%             pcolor(t,f,20*log10(abs(STFT)))
%             title(['event #' num2str(ee) ', channel = ' num2str(chn)])
%             colorbar
%             caxis([-80 -40])
%             pause
%             clf;
        end
    end
end

nn = 1;
%figure
for ee = 1:length(event_inds)
    if ~isnan(event_inds(ee,1))
        for chn = 1:4
            segment = yevent(chn,event_inds(ee,1)-pad(ee):event_inds(ee,2)+pad(ee));
            [STFT, f, t] = stft(segment,hanning(64),32,64,1000);
            [~,yind] = find(STFT==max(max(STFT)));
            STFT = circshift(STFT,size(STFT,2)/2-yind,2);
            ySTFT_mat(:,:,nn) = STFT;
            nn = nn+1;
%             pcolor(t,f,20*log10(abs(STFT)))
%             title(['event #' num2str(ee) ', channel = ' num2str(chn)])
%             colorbar
%             caxis([-80 -40])
%             pause
%             clf;
        end
    end
end

[sinv, tinv] = istft(mean(zSTFT_mat,3), hanning(64), hanning(64),32,64, 1000);
ztemplate_filt = filtfilt(bandfilt,sinv);
[sinv, tinv] = istft(mean(xSTFT_mat,3), hanning(64), hanning(64),32,64, 1000);
xtemplate_filt = filtfilt(bandfilt,sinv);
[sinv, tinv] = istft(mean(ySTFT_mat,3), hanning(64), hanning(64),32,64, 1000);
ytemplate_filt = filtfilt(bandfilt,sinv);
%% Matched Filtering

zevent = data_filt(:,[1 4 7 10]).';
xevent = data_filt(:,[2 5 8 11]).';
yevent = data_filt(:,[3 6 9 12]).'; 
zmatched_data = zeros(size(zevent));
xmatched_data = zeros(size(xevent));
ymatched_data = zeros(size(yevent));

for i = 1:4
full_nfft = 5399910;
Fdata = fft(zevent(i,:), full_nfft);
Freplica = conj(fft(ztemplate_filt, full_nfft));
full_matched_filter = repmat(Freplica,1,1);
zmatched_data(i,:) = bsxfun(@rdivide, ifft(Fdata .* full_matched_filter), (sqrt(length(ztemplate_filt)*length(zevent(i,:)))*sqrt(var(zevent(i,:)))*sqrt(var(ztemplate_filt))));
%figure;
%plot(matched_data(i,:));
end

for i = 1:4
full_nfft = 5399910;
Fdata = fft(xevent(i,:), full_nfft);
Freplica = conj(fft(xtemplate_filt, full_nfft));
full_matched_filter = repmat(Freplica,1,1);
xmatched_data(i,:) = bsxfun(@rdivide, ifft(Fdata .* full_matched_filter), (sqrt(length(xtemplate_filt)*length(xevent(i,:)))*sqrt(var(xevent(i,:)))*sqrt(var(xtemplate_filt))));
%figure;
%plot(matched_data(i,:));
end

for i = 1:4
full_nfft = 5399910;
Fdata = fft(yevent(i,:), full_nfft);
Freplica = conj(fft(ytemplate_filt, full_nfft));
full_matched_filter = repmat(Freplica,1,1);
ymatched_data(i,:) = bsxfun(@rdivide, ifft(Fdata .* full_matched_filter), (sqrt(length(ytemplate_filt)*length(yevent(i,:)))*sqrt(var(yevent(i,:)))*sqrt(var(ytemplate_filt))));
%figure;
%plot(matched_data(i,:));
end