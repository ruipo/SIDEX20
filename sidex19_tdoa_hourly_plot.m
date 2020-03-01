% Made for SIDEX data collected in 2018-2019
% Requires the data directory structure of ~/sidex/sidex_DATE/sidex_TIME.txt
% with required code in ~/sidex/code/ folder

% RUN code in ~/sidex/ folder
% can be modified to analyzed similar datasets
% RUN after sidex19_tdoa_hourly_run.m

%% Compile all results
set(0,'defaultFigureVisible','off')
xpos = [-8.2; -19.6; 23]; % position of geophones
ypos = [-23.9; 10.8; -11.9];

n = 50; % set boundary box for plot (n = 50 means a 100mby100m square)
x = -n:1:n;
y = x;

c_mos = zeros(2*n+1,2*n+1);

folder = dir('~/sidex/*stats'); %folder where event localization results are stored
for fold = 1:length(folder) %2:17% 18:38; %39:48;
    fold
    
    event_loc_all = [];
    timestamp_all = [];
    
    cd(folder(fold).name) % go into a folder
    hrs = dir(['~/sidex/' folder(fold).name '/hour_*']); % for each hour of data results within this folder
    
    for h = 1:length(hrs) % sort hours to be chronological
        if length(hrs(h).name)<15
            prev = hrs(h).name(6:end);

            hrs(h).name(6) = '0';
            hrs(h).name(7:15) = prev;
        end

        [~,ind] = sort({hrs.name});
        hrs = hrs(ind);
    end
    
    for h = 1:length(hrs) % For each hour
        
        try load(['~/sidex/' folder(fold).name '/' hrs(h).name '/tdoa_allaxis.mat']) % load in localization data (events found all axes in this case)
    
            event_loc_all = [event_loc_all loc_est_all];
            timestamp_all = [timestamp_all timestamp];
            
        catch
            disp('no event in this hour')
                        
        end
    end
    
    if ~isempty(timestamp_all) % if there are events in this folder
    c_tot = zeros(2*n+1,2*n+1);
    cplot = zeros(2*n+1,2*n+1); 
    date_start = datetime(timestamp_all(1),'convertFrom','datenum');
    date_end = datetime(timestamp_all(end),'convertFrom','datenum');
    
    timelist = datenum(linspace(date_start,date_end,minutes(date_end-date_start))); %create a timeseries from start of first event to end of last event by 1min increments

    xlim([-n n]);
    ylim([-n n]);
    
    event_num = 1;

    d = timelist(2)-timelist(1);
    for t = 1:length(timelist)-1 % For each timestep in timeseries

        event_loc_p = [];

        while timelist(t) <= timestamp_all(event_num) && timestamp_all(event_num) < timelist(t)+d % if events occurred in this time step
            
            event_loc_p = [event_loc_p event_loc_all(:,event_num)]; %add these events to a list

            event_num = event_num+1;
        end
        
        if isempty(event_loc_p) % if this list is empty, then c for this timestep is 0; cplot from previous timestep is faded 50% by multiplying it by 0.5
            c = zeros(2*n+1,2*n+1);
            cplot =cplot*0.5;
        else % else, color in the location of these events in c, colormap depends on number of events that occurred at any location. 
            disp([num2str(t) '/' num2str(length(timelist))])
            indx = event_loc_p(1,:)<=n;
            event_loc_p = event_loc_p(:,indx);
            indx = event_loc_p(1,:)>=-n;
            event_loc_p = event_loc_p(:,indx);
            indy = event_loc_p(2,:)<=n;
            event_loc_p = event_loc_p(:,indy);
            indy = event_loc_p(2,:)>=-n;
            event_loc_p = event_loc_p(:,indy);
            event_loc_p = [event_loc_p [n-0.5;n-0.5] [-n+0.5;-n+0.5]];
            c = hist3(event_loc_p.','Nbins',[2*n 2*n],'CdataMode','auto');
            c(1,1) = c(1,1)-1;
            c(end,end) = c(end,end)-1;
            c(:,end+1) = 0;
            c(end+1,:) = 0;

            cplot = cplot*0.75+c; % cplot is 75% of cplot from last timestep + c of this timestep.
            c_tot = c_tot+c; % ctot is the culmulatvie sum of all c's from each timestep (records the number of events at each location for the entire folder)
        end
        
        if max(max(cplot>0.1)) == 1 % save cplot of this timestep (1 min total)
        figure('units','normalized','outerposition',[0 0 1 1])
        title(datestr(datetime(timelist(t),'convertFrom','datenum')))
        hold on
        xlabel('X position (m)')
        ylabel('Y position (m)')
        set(gca,'fontsize',20)
        caxis([0 3])
        colormap cool
        h = pcolor(x,y,cplot);
        set(h, 'edgecolor','none')
        colorbar
        plot(xpos,ypos,'ko','MarkerFaceColor', 'k','MarkerSize',10);
        set(gca,'layer','top')
        xL = xlim;
        yL = ylim;
        line([0 0], yL,'color','black');
        line(xL, [0 0],'color','black');
        pause(0.1)
        saveas(gcf,['/home/ruic/Desktop/sidex_tdoa_figs_dec/' datestr(datetime(timelist(t),'convertFrom','datenum')) '.png']);
        %pause
        close all
        end
    end
        
    figure('units','normalized','outerposition',[0 0 1 1]) %save ctot of entire folder (daily total)
    h1 = pcolor(x,y,c_tot);
    colormap cool
    set(h1,'edgecolor','none')
    colorbar
    caxis([0 50])
    xlabel('X position (m)')
    ylabel('Y position (m)')
    title(folder(fold).name(1:10))
    set(gca,'fontsize',20)
    grid on
    xlim([-n n]);
    ylim([-n n])
    hold on
    plot(xpos,ypos,'ko','MarkerFaceColor', 'k','MarkerSize',10);
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'color','black');
    line(xL, [0 0],'color','black');
    saveas(gcf,['/home/ruic/Desktop/sidex_tdoa_figs_daily/' folder(fold).name(1:10) '.png']);
        
    close all
    
    c_mos = c_mos+c_tot; %add c_tot to c_mos to get monthly total
    end
    
    cd('~/sidex/')
end

save('/home/ruic/Desktop/sidex_tdoa_figs_daily/dec_c_mos','c_mos');

set(0,'defaultFigureVisible','on')
        
        
%% make videos

testf = dir('~/Desktop/sidex_tdoa_figs_dec/*.png'); % path to folder with daily/minute figures

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([testf(1).folder 'dec_minute.avi']); % name output avi file here!
  writerObj.FrameRate = 5;
  open(writerObj);
  
  for frame = 1 : length(testf)
    % Construct an output image file name.
    outputBaseFileName = sprintf(testf(frame).name);
    outputFullFileName = fullfile(testf(frame).folder, outputBaseFileName);
    % Read the image in from disk.
    thisFrame = imread(outputFullFileName);
    % Convert the image into a "movie frame" structure.
    recalledMovie(frame) = im2frame(thisFrame);
    % Write this frame out to a new video file.
    writeVideo(writerObj, thisFrame);
  end
  close(writerObj);
  % Get rid of old image and plot.
