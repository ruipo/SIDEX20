
freq3 = linspace(YMAX,YMIN,358);
kx3 = linspace(XMIN/1000, XMAX/1000, 151);

figure
h = pcolor(kx3,freq3, flipud(DATA3.'));
set(h,'Edgecolor','None')
colormap jet
set(gca,'fontsize',20)
xlabel('Wavenumber (1/m)')
ylabel('Frequency (Hz)')
title('1m Ice, SD = 0.1m; \omega-k Diagram')
caxis([-100 40]);


%%
freqmax = zeros(length(kx3),1);
flist = fliplr(freq3);
for i = 1:size(DATA3,1)
    [~,ind] = findpeaks(DATA3(i,:),'MinPeakProminence',1);
    if ~isempty(ind)
        freqmax(i) = flist(max(ind));
    end
end

%%

freqmax_1_10 = [1.953 2.929 3.906 4.882 5.859 6.835];
kx_1_10 = [0.1669 0.2131 0.2453 0.2654 0.2855 0.3116];

kx_200_500 = kx3(63:138);
freqmax_200_500 = freqmax(63:138);
%%

kx = [kx_1_10 kx_10_200 kx_200_500];
freqmax = [freqmax_1_10.'; freqmax_10_200; freqmax_200_500];

kxs = linspace(kx(1),kx(end),100);
freqmaxs = interp1(kx,freqmax,kxs,'spline');

slope = zeros(length(freqmaxs),1);
for j = 2:length(freqmaxs)-1
    slope(j-1) = (freqmaxs(j+1) - freqmaxs(j-1))/(kxs(j+1) - kxs(j-1));
end

figure

subplot(1,2,1)
title('Mode 1 \omega-k Diagram')
plot(kxs,freqmaxs,'linewidth',2);
grid on
xlabel('Wavenumber (1/m)')
ylabel('Frequency (Hz)')
set(gca,'fontsize',20)

subplot(1,2,2)
title('Estimated Group Velocity')
hold on
semilogx(freqmaxs(1:end-2),slope(1:end-2),'ro');
grid on
ylabel('Group Velocity (m/s)')
%yticks([0 25 50 75 100 125 150 175 200 225 250])
%set(gca,'Yticklabel',[0 25 50 75 100 125 150 175 200 225 250])
%xticks([0 50 100 150 200 250 300 350 400 450 500])
%set(gca,'Xticklabel',[0 50 100 150 200 250 300 350 400 450 500])
xlabel('Frequency (Hz)')
set(gca,'fontsize',20)
h1 = lsline;
h1.Color = 'k';
h1.LineWidth = 2;

