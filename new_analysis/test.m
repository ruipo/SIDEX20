%% Data Band-pass filter
f1 = 16;
f2 = 32;
FS = 1000;

bandfilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',f1,'CutoffFrequency2',f2,'SampleRate',FS);
data_filt = filtfilt(bandfilt,data);
event_inds = xlsread('event_inds.xlsx',1,'B2:C189');

figure
plot(t,data_filt(:,1))
hold on
plot(t,data_filt(:,4))
plot(t,data_filt(:,7))
plot(t,data_filt(:,10))
plot(epochtime_event,zeros(length(epochtime_event),1),'*')

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

%% Motion Product Detector

FS = 1000;
f1 = 16;
f2 = 32;
thetas = [275 275 300 311];
%thetas = [100 100 100 100];

figure

for ee = 1:188
    if ~isnan(event_inds(ee,1))
        [event1x, event1y]= ll2xy(lat_event(ee),long_event(ee),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
        geophone_loc = [24.36 -36.64;-37.71 11.51;21.05 25.30;-7.71 -0.167];
        event = data(event_inds(ee,1):event_inds(ee,2),:);
        plot(geophone_loc(:,1),geophone_loc(:,2),'ro')
        hold on
        plot(event1x, event1y, 'b*')

        for chn = 1:4
            datain = event(:,chn+(chn-1)*2:chn+(chn-1)*2+2);
            [X_HiV, Y_HiV, X_est, Y_est] = MPD(datain,f1,f2,FS,10,thetas(chn));
            plot(X_est+geophone_loc(chn,1),Y_est+geophone_loc(chn,2),'k');
        end
        grid on
        title(['event number: ' num2str(ee)]);
        xlim([-250 250])
        ylim([-250 250])
        pause
        clf;
    end
end

%% plot V-H

FS = 1000;
f1 = 16;
f2 = 32;
tdoa_xy = zeros(4,188);
figure

for ee = 1:188
    if ~isnan(event_inds(ee,1))
        event = data_filt(event_inds(ee,1):event_inds(ee,2),:);

        for chn = 1:4
            datain = event(:,chn+(chn-1)*2:chn+(chn-1)*2+2);
            [corr,lags] = xcorr(abs(hilbert(datain(:,2))),abs(hilbert(datain(:,3))));
            [~,loc] = max(corr);
            tdoa_xy(chn,ee) = lags(loc)/FS;
        
            subplot(1,2,1)
            plot(datain(:,2),datain(:,1))
            xlabel('X-dir')
            ylabel('Z-xir')
            grid on
            axis equal
            title(['event ' num2str(ee) 'Geophone ' num2str(chn)])
            subplot(1,2,2)
            plot(datain(:,3),datain(:,1))
            xlabel('Y-dir')
            ylabel('Z-xir')
            grid on
            axis equal
            title(['Geophone ' num2str(chn)])
            pause

            clf
        end
        

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
%% Event Locatization

zevent = data_filt(:,[1 4 7 10]).';
xevent = data_filt(:,[2 5 8 11]).';
yevent = data_filt(:,[3 6 9 12]).'; 
geophone_loc = [24.36 -36.64;-37.71 11.51;21.05 25.30;-7.71 -0.167];

FS = 1000;
c_range = 350;
N = 1000;
plotting = 1;

figure
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

pos_est_mat = zeros(2,length(event_inds));
c_est_mat = zeros(length(event_inds),1);
error_mat = zeros(length(event_inds),1);

for i = 185
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

        try
            [loc_est,c_est,err,tdoa_mat] = loc_est_calib_temp(zevent,xevent,yevent,geophone_loc(:,1),geophone_loc(:,2),ss,es,FS,c_range,N,plotting,calib_act);
        catch
            disp('Skipping this...')
        end

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
h(3) = plot(NaN,NaN,'r*','MarkerSize',6);
h(4) = plot(NaN,NaN,'b*','MarkerSize',6);
h(5) = plot(NaN,NaN,'m*','MarkerSize',6);
h(6) = plot(NaN,NaN,'g*','MarkerSize',6);
legend(h, 'Geophone Locations','True Location','Estimated Location (v < 250m/s)','Estimated Location (250m/s <= v < 500m/s)','Estimated Location (500m/s <= v < 750m/s)','Estimated Location (v >= 750m/s)');
set(0,'DefaultLegendAutoUpdate','off')

for i = 1:188
    i
    [eventx, eventy]= ll2xy(lat_event(i),long_event(i),mean(geophone_GPS(:,1)),mean(geophone_GPS(:,2)));
    calib_act = [eventx eventy];
    
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
    
    title([datestr(datetime(epochtime_event(i), 'convertfrom','posixtime')) ' UTC'])

    plot(eventx,eventy,'ks','MarkerSize',8,'linewidth',1.25);
    plot(pos_est_mat(1,i),pos_est_mat(2,i),[color '*'],'MarkerSize',6,'linewidth',1.25);
    quiver(pos_est_mat(1,i),pos_est_mat(2,i),eventx-pos_est_mat(1,i),eventy-pos_est_mat(2,i),0,'color',[1 0.5 1]);
    saveas(gcf,[path,'loc_est_figs/' num2str(i) '.png']);
    pause(0.1)
end

%% Location Estimation Results movie

testf = dir([path 'loc_est_figs/*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([path 'loc_est_figs/loc_est.avi']);
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


