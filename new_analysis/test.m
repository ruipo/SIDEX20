%% Data Band-pass filter
f1 = 16;
f2 = 32;
FS = 1000;

bandfilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',f1,'CutoffFrequency2',f2,'SampleRate',FS);
data_filt = filtfilt(bandfilt,data);


figure
plot(t,data_filt(:,1))
hold on
plot(t,data_filt(:,4))
plot(t,data_filt(:,7))
plot(t,data_filt(:,10))
plot(epochtime_event,zeros(length(epochtime_event),1),'*')

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
thetas = [300 275 300 311];
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