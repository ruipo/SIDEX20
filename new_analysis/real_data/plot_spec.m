%% Plotting
% t1 = 59000;%1/26
% t2 = 61000;%1/26
t1 = 19600;%1/31
t2 = 20400;%1/31
figure
for j = 1:12
plot(t(t1:t2),data(t1:t2,j)+(j-1)*0.75*max(max(abs(data(t1:t2,1:12)))))
hold on
xlabel('Time')
ylabel('Amplitude')
%datetick('x','HH:MM:SS:FFF','keepticks');
xlim([t(t1) t(t2)]);
title(date);
set(gca,'fontsize',20);
grid on
end
%%
figure
plot(t(t1:t2),data(t1:t2,[1 4 7 10]))
hold on
xlabel('Time')
ylabel('Amplitude')
title(date);
set(gca,'fontsize',20);
grid on

%%
chn = 1;
FS = 1000;
L = 128; % Window length (1ms)
R = L/2; % Overlap percentage
NFFT = L; 
time = 0:1/FS:size(data,1)/FS - 1/FS;

[S1,f,t] = spectrogram(data(:,chn),hanning(L),R,NFFT,FS);
[S2,f,t] = spectrogram(data(:,chn+3),hanning(L),R,NFFT,FS);
[S3,f,t] = spectrogram(data(:,chn+6),hanning(L),R,NFFT,FS);
[S4,f,t] = spectrogram(data(:,chn+9),hanning(L),R,NFFT,FS);

figure
subplot(4,1,1)
imagesc(t,f,20*log10(abs(S1)))
set(gca,'YDir','normal')
colormap jet
caxis([-100 -0])
colorbar
xlim([19 23])
ylim([0 100])
subplot(4,1,2)
imagesc(t,f,20*log10(abs(S2)))
set(gca,'YDir','normal')
colormap jet
caxis([-100 -0])
colorbar
xlim([19 23])
ylim([0 100])
subplot(4,1,3)
imagesc(t,f,20*log10(abs(S3)))
set(gca,'YDir','normal')
colormap jet
caxis([-100 -0])
colorbar
xlim([19 23])
ylim([0 100])
subplot(4,1,4)
imagesc(t,f,20*log10(abs(S4)))
set(gca,'YDir','normal')
colormap jet
caxis([-100 -0])
colorbar
xlim([19 23])
ylim([0 100])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

