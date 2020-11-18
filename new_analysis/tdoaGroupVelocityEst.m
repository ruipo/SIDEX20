
function [vg_est] = tdoaGroupVelocityEst(data_samp, time_samp, fcs, bw, range, FS)

vg_est = zeros(length(fcs),factorial(size(data_samp,2)-1));

for ff = 1:length(fcs)
    fc = fcs(ff);
    f1 = fc-bw/2;
    f2  = fc+bw/2;
    
    if f1 <= 0
        f1 = 1;
    end
    
    if f2 >= FS/2
        f2 = FS/2-1;
    end
    
    bandfilt = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',f1,'CutoffFrequency2',f2,'SampleRate',FS);

    data_filt = filtfilt(bandfilt,data_samp);

    hilt = hilbert(data_filt);
    
    vg = [];
    for chn1 = 1:size(data_filt,2)
        for chn2 = chn1+1 :size(data_filt,2)
        
            [~,ind1] = max(abs(hilt(:,chn1)));
            te1 = time_samp(ind1);
            r1 = range(chn1);
            
            [~,ind2] = max(abs(hilt(:,chn2)));
            te2 = time_samp(ind2);
            r2 = range(chn2);
            
            tdiff = abs(te2-te1);
            rdiff = abs(r2-r1);
            
            vg = [vg rdiff/tdiff];
        end   
    end
    
    vg_est(ff,:) = vg;
  
end