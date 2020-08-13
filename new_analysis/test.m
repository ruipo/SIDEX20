%% Data Band-pass filter
f1 = 16;
f2 = 32;
FS = 1000;

bandfilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',f1,'CutoffFrequency2',f2,'SampleRate',FS);
data_filt = filtfilt(bandfilt,data);


figure
plot(data_filt(:,1))
hold on
plot(data_filt(:,4))
plot(data_filt(:,7))
plot(data_filt(:,10))
plot(epochtime_event,zeros(length(epochtime_event),1),'*')

event1 = data(364900:366400,:);
event1_time = t(364900:366400);
event2 = data(367700:369000,:);
event2_time = t(367700:369000);
event3 = data(370600:371700,:);
event3_time = t(370600:371700);
event4 = data(373000:374200,:);
event4_time = t(373000:374200);

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
vgs_est(159) = NaN;
vgs_est(vgs_est>100) = NaN;

cf = fit([long_event(~isnan(vgs_est)) lat_event(~isnan(vgs_est))],vgs_est(~isnan(vgs_est)),'lowess');
% figure
% plot([long_event(~isnan(vgs_est)) lat_event(~isnan(vgs_est))],vgs_est(~isnan(vgs_est)))

[cx,cy] = meshgrid(-156.4096:0.0001:-156.3998, 71.3327:0.00001:71.3360);
c_fit = zeros(size(cx));
for i = 1:size(cx,1)
    c_fit(i,:) = cf([cx(i,:).' cy(i,:).']);
end
% 
c_fit(c_fit>50) = NaN;
c_fit(c_fit<10) = NaN;
figure
h = pcolor(cx(1,:),cy(:,1),c_fit);
% caxis([50 1000])
set(h,'Edgecolor','None')
hold on
plot(geophone_GPS(:,2),geophone_GPS(:,1),'ro', 'MarkerFaceColor', 'r');
xlabel('Longitude')
ylabel('Latitude')
xlim([-156.4096 -156.3998]);
ylim([71.3327 71.3360]);
xL = xlim;
yL = ylim;
plot(long_event(~isnan(vgs_est)),lat_event(~isnan(vgs_est)),'k*')
title('Propagation Speed Estimation')
colormap(flipud(jet))
colorbar
set(gca,'fontsize',20)


%% Motion Product Detector

FS = 1000;
f1 = 16;
f2 = 32;

geophone_loc = [24.36 -36.64;-37.71 11.51;21.05 25.30;-7.71 0.167];

figure
plot(geophone_loc(:,1),geophone_loc(:,2),'ro')
hold on
for chn = 1:4
    datain = event1(:,chn+(chn-1)*2:chn+(chn-1)*2+2);
    [X_HiV, Y_HiV, X_est, Y_est] = MPD(datain,f1,f2,FS,100);
    plot(X_est+geophone_loc(chn,1),Y_est+geophone_loc(chn,2),'k')
end
xlim([-50 50])
ylim([-50 50])
% figure
% plot(X_HiV,Y_HiV,'.')
% xlim([-2.5E-8 2.5E-8])
% ylim([-2.5E-8 2.5E-8])
% hold on
% plot(X_est,Y_est)

