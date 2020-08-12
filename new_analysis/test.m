%% Group Velocity Estimation

FS = 1000;
f1 = 16;
f2 = 32;
vg_est = [];

for chn = 1:12
    datain = data_samp(:,chn);
    source_time = epochtime_event(1);
    
    if chn <= 3
        dist = calib_dist(1,1);
    elseif 3 < chn <= 6
        dist = calib_dist(2,1);
    elseif 6 < chn <= 9
        dist = calib_dist(3,1);
    else
        dist = calib_dist(4,1);
    end
        
    vg_est(chn) = groupVelocityEst(datain, time_samp, f1, f2, FS,source_time, dist);
end


%% Motion Product Detector

FS = 1000;
f1 = 16;
f2 = 32;

geophone_loc = [24.36 -36.64;-37.71 11.51;21.05 25.30;-7.71 0.167];

figure
plot(geophone_loc(:,1),geophone_loc(:,2),'ro')
hold on
for chn = 1:4
    datain = data_samp(:,chn+(chn-1)*2:chn+(chn-1)*2+2);
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

