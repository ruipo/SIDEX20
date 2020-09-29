
%% SP1-V

filename = 'sidex_exp_sp5_vertical.txt';
fout = read_moveout_asc(filename);
TT = array2timetable(fout,'SampleRate',1024);
for i = 1:length(TT.Properties.VariableNames)
    TT.Properties.VariableNames(i) = {[num2str((i-1)*10) 'm']};
end
smallTTV = TT(1:1024,fliplr([2 4 6 8 10 12 14 16 18 20 22]));

filename = 'sidex_exp_sp5_horiz.txt';
fout = read_moveout_asc(filename);
TT = array2timetable(fout,'SampleRate',1024);
for i = 1:length(TT.Properties.VariableNames)
    TT.Properties.VariableNames(i) = {[num2str((i-1)*10) 'm']};
end
smallTTH = TT(1:1024,fliplr([2 4 6 8 10 12 14 16 18 20 22]));

figure
subplot(1,2,1)
stackedplot(smallTTV)
set(gca,'fontsize',12)
title('Explosive Source, SP5-V')
grid on
subplot(1,2,2)
stackedplot(smallTTH)
set(gca,'fontsize',12)
title('Explosive Source, SP5-H')
grid on

%% spec plot

chn = 1;
FS = 1024;
L = 32; % Window length (1ms)
R = L/2; % Overlap percentage
NFFT = L; 

[S1,f,t] = spectrogram(smallTTV(:,end-10).Variables,hanning(L),R,NFFT,FS);
[S2,f,t] = spectrogram(smallTTV(:,end-8).Variables,hanning(L),R,NFFT,FS);
[S3,f,t] = spectrogram(smallTTV(:,end-6).Variables,hanning(L),R,NFFT,FS);
[S4,f,t] = spectrogram(smallTTV(:,end-4).Variables,hanning(L),R,NFFT,FS);
[S5,f,t] = spectrogram(smallTTV(:,end-2).Variables,hanning(L),R,NFFT,FS);
[S6,f,t] = spectrogram(smallTTV(:,end).Variables,hanning(L),R,NFFT,FS);

figure
subplot(6,1,1)
imagesc(t,f,20*log10(abs(S1)))
set(gca,'YDir','normal')
colormap jet
caxis([-70 -20])
colorbar
%xlim([19 23])
%ylim([0 200])
ylabel({'range = 210m', 'Frequency (Hz)'})
title('Explosive Source, SP5-V')
set(gca,'fontsize',12)
subplot(6,1,2)
imagesc(t,f,20*log10(abs(S2)))
set(gca,'YDir','normal')
colormap jet
caxis([-70 -20])
colorbar
%xlim([19 23])
%ylim([0 200])
ylabel({'range = 170m', 'Frequency (Hz)'})
set(gca,'fontsize',12)
subplot(6,1,3)
imagesc(t,f,20*log10(abs(S3)))
set(gca,'YDir','normal')
colormap jet
caxis([-70 -10])
colorbar
%xlim([19 23])
%ylim([0 200])
ylabel({'range = 130m', 'Frequency (Hz)'})
set(gca,'fontsize',12)
subplot(6,1,4)
imagesc(t,f,20*log10(abs(S4)))
set(gca,'YDir','normal')
colormap jet
caxis([-80 -0])
colorbar
%xlim([19 23])
%ylim([0 200])
ylabel({'range = 90m', 'Frequency (Hz)'})
set(gca,'fontsize',12)
subplot(6,1,5)
imagesc(t,f,20*log10(abs(S5)))
set(gca,'YDir','normal')
colormap jet
caxis([-85 10])
colorbar
%xlim([19 23])
%ylim([0 200])
ylabel({'range = 50m', 'Frequency (Hz)'})
set(gca,'fontsize',12)
subplot(6,1,6)
imagesc(t,f,20*log10(abs(S6)))
set(gca,'YDir','normal')
colormap jet
caxis([-90 20])
colorbar
%xlim([19 23])
%ylim([0 200])
xlabel('Time (s)')
ylabel({'range = 10m', 'Frequency (Hz)'})
set(gca,'fontsize',12)

%% vg_est
fvals = f(2:5);
ranges = [0 10 50 90 130 170 210]; 
[~,inds] = max(abs(S6(2:5,:)).');
t6 = t(inds);
[~,inds] = max(abs(S5(2:5,:)).');
t5 = t(inds);
[~,inds] = max(abs(S4(2:5,:)).');
t4 = t(inds);
[~,inds] = max(abs(S3(2:5,:)).');
t3 = t(inds);
[~,inds] = max(abs(S2(2:5,:)).');
t2 = t(inds);
[~,inds] = max(abs(S1(2:5,:)).');
t1 = t(inds);
tmat = [t6;t5;t4;t3;t2;t1];

figure
plot([0; tmat(:,1)],ranges,'b*')
hold on
plot([0; tmat(:,1)],([0; tmat(:,1)]\ranges.')*[0; tmat(:,1)],'-b')
plot([0; tmat(:,2)],ranges,'r*')
plot([0; tmat(:,2)],([0; tmat(:,2)]\ranges.')*[0; tmat(:,2)],'-r')
plot([0; tmat(:,3)],ranges,'g*')
plot([0; tmat(:,3)],([0; tmat(:,3)]\ranges.')*[0; tmat(:,3)],'-g')
plot([0; tmat(:,4)],ranges,'k*')
plot([0; tmat(:,4)],([0; tmat(:,4)]\ranges.')*[0; tmat(:,4)],'-k')
grid on
title('Explosive Source, SP5-V')
set(gca,'Fontsize',20)
xlabel('Arrival Time (sec)')
ylabel('Range (m)')
legend('32Hz',strcat('v_{est} = ',num2str(([0; tmat(:,1)]\ranges.')),'m/s'),'64Hz',strcat('v_{est} = ',num2str(([0; tmat(:,2)]\ranges.')),'m/s'),'96Hz',strcat('v_{est} = ',num2str(([0; tmat(:,3)]\ranges.')),'m/s'),'128Hz',strcat('v_{est} = ',num2str(([0; tmat(:,4)]\ranges.')),'m/s')); 
%% SP3-V

filename = 'sidex_vs_sp3_vertical.txt';
fout = read_moveout_asc(filename);
TT = array2timetable(fout,'SampleRate',2048);
for i = 1:length(TT.Properties.VariableNames)
    TT.Properties.VariableNames(i) = {[num2str((i-1)*10) 'm']};
end

smallTT = TT(1:512,fliplr([2 4 6 8 10 12 14 16 18 20 22]));
figure
stackedplot(smallTT)
set(gca,'fontsize',15)
title('Vertical Source, SP3-V')
grid on

%% SP5-V

filename = 'sidex_vs_sp5_vertical.txt';
fout = read_moveout_asc(filename);
TT = array2timetable(fout,'SampleRate',2048);
for i = 1:length(TT.Properties.VariableNames)
    TT.Properties.VariableNames(i) = {[num2str((i-1)*10) 'm']};
end

smallTT = TT(1:512,fliplr([2 4 6 8 10 12 14 16 18 20 22]));
figure
stackedplot(smallTT)
set(gca,'fontsize',15)
title('Vertical Source, SP5-V')
grid on