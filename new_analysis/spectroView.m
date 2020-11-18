
function [S_mat,f,t] = spectroView(data, FS, win_len)

L = win_len; % Window length (0.1s)
R = L/2; % Overlap percentage
NFFT = L; 
S_mat = [];

for chn = 1:size(data,2)
    
    [S,f,t] = spectrogram(data(:,chn),hanning(L),R,NFFT,FS);
    S_mat(:,:,chn) = S;
    
end

end
