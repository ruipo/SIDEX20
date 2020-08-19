
% Plots the power spectral density of the input data.
% data must be [# of data points] by [# of channels].
% Can input just one channel or multiple channels in a list.
% Need to specify window size, type of window
% (hanning(),hamming(),rect(),etc...), amount of overlap (0-1),NFFT
% size (atleast 2*windowsize-1), and FS.

function [psd,psd_avg,psd_var,freq] = psd(data,window_size,window,overlap,NFFT,FS,channels,data_name)
    psd_var = 0;
    % One channel input
    if length(channels) == 1
        
        chn = channels;
        ts = data(:,chn);
        
        % Zero pad if ts is not the length of a power of 2
        data_len = 2^nextpow2(length(data));
        zero_len = data_len - length(data);
        ts = [ts;zeros(zero_len,1)];
        
        window_start = round(window_size-window_size*overlap);
        num_window = round(length(ts)/window_start)-1;
        ts_f = zeros(num_window,NFFT);
        ts_f_onesided = zeros(num_window,NFFT/2+1);
        psd = zeros(num_window,NFFT/2+1);
        
        for l = 1:num_window
            %disp([num2str(l), '/', num2str(num_window)])
            ts_f(l,:) = fft(window.*ts(l*window_start-window_start+1:l*window_start-window_start+window_size),NFFT);
            ts_f_onesided(l,:) = ts_f(l,1:NFFT/2+1);
            psd(l,:) = (1/(FS*norm(window,2)^2)) * abs(ts_f_onesided(l,:)).^2;
            psd(l,2:end-1) = 2*psd(l,2:end-1);
        end
        
        psd(~any(psd,2), :) = [];
        psd_avg = median(psd,1);
        psd_var = var(psd,1);
        freq = 0:FS/NFFT:FS/2;
        
        figure
        semilogx(freq,10*log10(psd_avg/(1E-6)^2),'linewidth',1.5)
%         hold on
%         semilogx(freq,10*log10(psd_var/(1E-6)^2),'linewidth',1.5)
%         grid on
    end

    % Multiple channle inputs
    if length(channels) > 1
%         figure
        psd_avg_chn = zeros(NFFT/2+1,length(channels));
        
        for i = 1:length(channels)
           chn = channels(i);
           ts = data(:,chn);

            % Zero pad if ts is not the length of a power of 2
            data_len = 2^nextpow2(length(data));
            zero_len = data_len - length(data);
            ts = [ts;zeros(zero_len,1)];

            window_start = round(window_size-window_size*overlap);
            num_window = round(length(ts)/window_start)-1;
            ts_f = zeros(num_window,NFFT);
            ts_f_onesided = zeros(num_window,NFFT/2+1);
            psd = zeros(num_window,NFFT/2+1);

            for l = 1:num_window
                ts_f(l,:) = fft(window.*ts(l*window_start-window_start+1:l*window_start-window_start+window_size),NFFT);
                ts_f_onesided(l,:) = ts_f(l,1:NFFT/2+1);
                psd(l,:) = (1/(FS*norm(window,2)^2)) * abs(ts_f_onesided(l,:)).^2;
                psd(l,2:end-1) = 2*psd(l,2:end-1);
            end

            psd_avg = mean(psd,1);
            psd_var = var(psd,1);
            psd_avg_chn(:,chn) = psd_avg;
            freq = 0:FS/NFFT:FS/2;
            
            semilogx(freq,10*log10(psd_avg/(1E-6)^2),'linewidth',1.5)
            grid on
            hold on
            
        end
        psd = psd_avg_chn;
    end

set(gca,'Fontsize',20);
title(['Power Spectral Density; Time = ',data_name,'; Window size = ', num2str(window_size)]);
xlabel('Frequency (Hz)');
xlim([1 FS/2]);
ylabel('Power/Frequency (dB/Hz)');

end