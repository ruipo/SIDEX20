

% INPUTS
%
% data		time series (time x depth)
% p		element positions matrix,  [x_e1,y_e1,z_e1;x_e2...]
% FS		sample frequency (Hz)
% elev		desired look angles in elevation in degrees
% az    desired look angles in azimuth in degrees
% c		sound speed
% f_range     frequency range of interest
% window    frequency window weighting
% overlap    percent overlap [0 1]
% weighting     spatial array weighting
% NFFT
%
% OUTPUTS
%
% beamform_output	beamformed result (time x elev x az x freq)

function[beamform_output,t,t_end] = beamform_3D(data,p,FS,elev,az,c,f_range,NFFT,window,overlap,weighting)

% Define variables
%window = window./(FS*norm(window,2)^2);
%window = window./sum(window);
N = size(p,1);
beam_elev = (90-elev).*(pi/180);
beam_az = az.*(pi/180);
win_len = length(window);
t_end = size(data,1)/FS;

for m = 2:size(data,2)
    window = [window,window(:,1)];
end

% Format data
% data_len = 2^nextpow2(size(data,1));
% zero_len = data_len - size(data,1);
% zero_mat = zeros(zero_len,size(data,2));
% data = [data;zero_mat];
window_start = round(win_len-win_len*overlap);
num_window = round(size(data,1)/window_start);
beamform_output = zeros(num_window,length(beam_elev),length(beam_az),NFFT);
t = zeros(num_window,1);

% FFT Data

f1 = f_range(1);
f2 = f_range(end);
w = exp(-1i*2*pi*(f2-f1)/(NFFT*FS));
a = exp(1i*2*pi*f1/FS);

ts_f_mat = zeros(NFFT,N,num_window);
for l = 1:num_window
    ts_f_mat(:,:,l) = (1/(sqrt(FS)*norm(window(:,1),2)))*sqrt(2)*czt(window.*data(l*window_start-window_start+1:l*window_start-window_start+win_len,:),NFFT,w,a);
    t(l) = ((l+1)*window_start-window_start+1)/FS;
end

% Start beamforming
flist = linspace(f_range(1),f_range(end),NFFT)';

% calculate k
k = 2*pi*flist./c;

% linear window
if strcmp(weighting,'uniform')
    win = ones(1,N);
    win = win./norm(win,2);
end

% icex window
if strcmp(weighting,'icex_hanning')
    win = hanning(42)';
    win(2) = []; win(3) = []; win(4) = []; win(5) = []; win(6) = [];
    win(end-1) = []; win(end-2) = []; win(end-3) = []; win(end-4) = []; win(end-5) = [];
    win = win./norm(win,2);
end

% simi xarray window
if strcmp(weighting,'simi_xarray_hanning')
    win = [0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9820 0.9293 0.7371 0.3694 0.9293  0.7371 0.2248  0.8992 0.7371  0.3162  0.9293 0.7371 0.2594 0.4984 0.6375 0.6341 0.4559 0.6014 0.5083 2.8988e-05 0.4980];  
    win = win./norm(win,2);
end

% Hamming window
 if strcmp(weighting,'hanning')
    win = hanning(N)';
    win = win./norm(win,2);
end

% build steering vectors
for j = 1:length(beam_az)
    for mm = 1:length(beam_elev)
        %mm
        % form steering vector
        steer = exp(1i * k *(sin(beam_elev(mm))*cos(beam_az(j))*p(:,1)'+sin(beam_elev(mm))*sin(beam_az(j))*p(:,2)'+cos(beam_elev(mm))*p(:,3)'));

        % apply weighting
        steer = steer.*(ones(size(k,1),1)*win);

        % beamform
        for l = 1:num_window
            b_elem = sum((conj(steer).*ts_f_mat(:,:,l)).');
            beamform_output(l,mm,j,:) = abs(b_elem).^2;
        end
    end
end

%t = t - t(1);
%t_end = t_end - t(1);
end