
%% set directory
prefix = '/Volumes/RUIC_Backup/SIDEX20_data/sidex_2020-01-27/';
directory = dir([prefix 'Sidex*.txt']);
date = directory(5).name(7:14); % records the date that the data was collected
data = [];
FS = 1000;
good_chn = [1:16]; %choose chns to be plotted
first_file = 1025; %choose the first file to be read
cycle = 0;

%% plot a portion of the directory

num_files = 1; %choose the number of datafiles to be plotted
data = [];
count = 0;
datenumlist = zeros(num_files,1);
t = [];

for i = first_file:length(directory)
    i
    time = directory(i).name([7:14 16:21]);
    timestring = [time(1:8) time(9:10),':',time(11:12),':',time(13:14)];
    datenumlist(i) = datenum(timestring, 'YYYYmmddHH:MM:SS');
    filename = [prefix directory(i).name];
    try 
        M = dlmread(filename, ',', 2, 0);
    catch 
        warning(['Could not read file ',num2str(i),'. Skipped to next file.']);
    end
    data = [data; M(:,1:16)];
    t = [t linspace(datenumlist(i),datenumlist(i)+1.1574e-08*size(M,1),size(M,1))];

    count = count + 1;

    if count >= num_files || i == length(directory) %stop if reached the number of minute of data wanted to plot or if reach the end of the directory
        break;
    end
end

%data(1,:) = [];

first_file = i+1;
cycle = cycle+1;

t = datetime(t,'ConvertFrom','datenum');
time = 0:1/FS:size(data,1)/FS - 1/FS; %plot data in seconds

%% Filter data if wanted

f = [0 0.01 0.02 0.2 0.21 1]; %5-50 Hz 
%f = [0 0.19 0.2 0.6 0.61 1]; %50-150 Hz
%f = [0 0.01 0.02 0.98 0.99 1]; %150-250 Hz
a = [0 0 1 1 0 0];
ord = 500;

b = firpm(ord,f,a);
hd = dfilt.dffir(b);
zplane(hd)
data = filter(hd,data);

%% Plotting

figure
for j = 1:12
plot(t,data(:,j)+(j-1)*0.75*max(max(abs(data(:,1:12)))))
hold on
xlabel('Time')
ylabel('Amplitude')
%datetick('x','HH:MM:SS:FFF','keepticks');
xlim([t(1) t(end)]);
title(date);
set(gca,'fontsize',20);
grid on
end


%% Spectrogram

Fc = 0;
FS = 1000;
L = 128; % Window length (1ms)
R = 64; % Overlap percentage
NFFT = L; 
time = 0:1/FS:size(data,1)/FS - 1/FS;

for chn = 1
    chn
%output spectrogram estimate of data
data1 = data;
data1(isnan(data1(:,good_chn(chn))),good_chn(chn)) = 0;
figure
spectrogram(data1(:,good_chn(chn)),hanning(L),R,NFFT,FS);
title(['Chn = ', num2str(good_chn(chn))])
set(gca,'fontsize',20);
end