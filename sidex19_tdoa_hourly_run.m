% Made for SIDEX data collected in 2018-2019
% Requires the data directory structure of ~/sidex/sidex_DATE/sidex_TIME.txt
% with required code in ~/sidex/code/ folder

% RUN code in ~/sidex/ folder
% can be modified to analyzed similar datasets
% RUN after sidex19_eSelect_hourly_run.m

disp('WARNING WILL OVERWRITE EXISTING FOLDERS IF ANALYSIS IS ALREADY DONE!!! PROCEED?') % run warning, press enter to continue

pause

set(0,'DefaultFigureVisible','off')
folder = dir('~/sidex/sidex*');
for fold = 10:length(folder)
    
    % set directory
    prefix = ['~/sidex/' folder(fold).name '/'];
    directory = dir([prefix 'Sidex_201*.txt']);
    date = directory(end).name(7:14); % records the date that the data was collected
    FS = 500;
    num_chn = 16;
    good_chn = [1:6 10:12];%choose chns to be plotted
    first_file = 1; %choose the first file to be read

    foldername = [prefix(end-10:end-1) '_stats'];

    cycle = 1;
    hr = 1;
    while first_file < length(directory) 
        data = zeros(1,num_chn);
        datenumlist = zeros(1000,1);
        t = [];

        for i = first_file:length(directory) % for each data file in directory
            time = directory(i).name(16:21); % record time of file
            timestring = [date(1:4),'-',date(5:6),'-',date(7:8),' ',time(1:2),':',time(3:4),':',time(5:6)];
            datenumlist(i) = datenum(timestring, 'yyyy-mm-dd HH:MM:SS');
            filename = [prefix directory(i).name];
            try 
                M = dlmread(filename, ',', 2, 0); % read in data in file
            catch 
                warning(['Could not read file ',num2str(i),'. Skipped to next file.']);
            end
            data = [data(:,1:16); M(:,1:16)];
            t = [t linspace(datenumlist(i),datenumlist(i)+6.9444e-04-(6.9444e-04/(FS*60)),size(M,1))];

            if t(end) >= t(1) + 6.9444e-04*65 || i == length(directory) %stop if imported an hour's worth of data
                disp(datetime(t(end),'ConvertFrom','datenum'))
                break;
            end

        end
        
        foldername = [prefix(end-10:end-1) '_stats']; % go into folder containing selected events for this hour
        cd('~/sidex/')
        cd(foldername)
        name = ['hour_',num2str(hr),'_eSelect'];
        cd(name)
        save(['~/sidex/' foldername '/' name '/files'],'first_file','i');
        
        hr = hr+1;
        first_file = i+1;
        
        load('e_Select.mat'); % load in event selection data
        
        events = sum(~isnan(tally_start),1);
        events_s = find(events==3); % event==3 takes only events found on all 3 axes. 
        
        if isempty(events_s) % if event_s not empty, localize the events
            continue
        else
            cd ..
            cd ..
            cd('code')
            
            xpos = [-8.2; -19.6; 23];
            ypos = [-23.9; 10.8; -11.9];
            
            c_list = [1200];
            N = 1000;
            plotting = 0;
            
            num_event = length(events_s);
            loc_est_all = zeros(2,num_event);
            timestamp = zeros(1,num_event);
            
            data(1,:) = [];
            for chn = 1:size(data,2)
                data(:,chn) = data(:,chn)-mean(data(:,chn));
                data(:,chn) = data(:,chn)./std(data(:,chn));
            end
        
        
            xdata = data(:,[1 4 10]).';
            ydata = data(:,[2 5 11]).';
            zdata = data(:,[3 6 12]).';
            
            
            for n = 1:num_event
                disp([num2str(n) ' / ' num2str(num_event)])
                
                loc_est = zeros(2,1);
                ss = min(tally_start(:,events_s(n)));
                es = max(tally_end(:,events_s(n)));
                
                [xloc_est,~,xerr] = loc_est_hyp(xdata,xpos,ypos,ss,es,FS,c_list,N,plotting);
                %pause(1)
                [yloc_est,~,yerr] = loc_est_hyp(ydata,xpos,ypos,ss,es,FS,c_list,N,plotting);
                %pause(1)
                [zloc_est,~,zerr] = loc_est_hyp(zdata,xpos,ypos,ss,es,FS,c_list,N,plotting);
                
                if xerr<1E-15 || isnan(xerr)
                    xerr = 100000;
                end
                
                if yerr<1E-15 || isnan(yerr)
                    yerr = 100000;
                end
                
                if zerr<1E-15 || isnan(zerr)
                    zerr = 100000;
                end
            
                w = [1/xerr 1/yerr 1/zerr]./norm([1/xerr 1/yerr 1/zerr]); % weight estimate of event location with err of x,y,z axes
                
                loc_est(1) = w*[xloc_est(1);yloc_est(1);zloc_est(1)];
                loc_est(2) = w*[xloc_est(2);yloc_est(2);zloc_est(2)]; 
                
                loc_est_all(:,n) = loc_est;                
                timestamp(n) = t(round((es+ss)/2));

            end
            
            temp = sum(isnan(loc_est_all));
            inds = find(temp>=1);
            loc_est_all(:,inds) = [];
       
            timestamp(inds) = [];
            
            cd ..
        
            save(['~/sidex/' foldername '/' name '/tdoa_allaxis'],'loc_est_all','timestamp'); % save localization data to folders
            
        end
        
        
    end
end
    
set(0,'DefaultFigureVisible','on')

% CAN NOW RUN toda_hourly_plot.m to visualize results.

%% Compile all results
xpos = [-8.2; -19.6; 23];
ypos = [-23.9; 10.8; -11.9];
event_loc_all = [];
timestamp_all = [];
figure(1)
plot(xpos,ypos,'ko');
xlabel('X position (m)')
xlabel('Y position (m)')
grid on
xlim([-50 50]);
ylim([-50 50])


folder = dir('~/sidex/*stats');
for fold = 2:length(folder)
    fold
    
    cd(folder(fold).name)
    hrs = dir(['~/sidex/' folder(fold).name '/hour_*']);
    
    for h = 1:length(hrs)
        if length(hrs(h).name)<15
            prev = hrs(h).name(6:end);

            hrs(h).name(6) = '0';
            hrs(h).name(7:15) = prev;
        end

        [~,ind] = sort({hrs.name});
        hrs = hrs(ind);
    end
    
    for h = 1:length(hrs)
        h
        try load(['~/sidex/' folder(fold).name '/' hrs(h).name '/tdoa_allaxis.mat'])
    
            event_loc_all = [event_loc_all loc_est_all];
            timestamp_all = [timestamp_all timestamp];

            
            for e = 1:size(loc_est_all,2)
                figure(1)
                plot(xpos,ypos,'ko');
                xlabel('X position (m)')
                xlabel('Y position (m)')
                plot(loc_est_all(1,e),loc_est_all(2,e),'*b')
                title([folder(fold).name(1:10) ' ' hrs(h).name(1:4) ' ' hrs(h).name(6:7)])
                pause(0.1)
                hold on
                set(gca,'fontsize',20)
                grid on
                xlim([-50 50]);
                ylim([-50 50]);
         
            end
            
            
            hold off
            pause
            
        catch
            disp('no event in this hour')
            figure(1)
            title([folder(fold).name(1:10) ' ' hrs(h).name(1:4) ' ' hrs(h).name(6:7)])
            pause(2)
            
            
        end
    end
    
    cd('~/sidex/')
end
     
%% Plotting
xpos = [-8.2; -19.6; 23];
ypos = [-23.9; 10.8; -11.9];

figure
hold on
plot(xpos,ypos,'ko');
xlabel('X position (m)')
xlabel('Y position (m)')
xlim([-250 250]);
ylim([-250 250]);
xL = xlim;
yL = ylim;
line([0 0], yL,'color','black');
line(xL, [0 0],'color','black');
grid on

for e = 1:size(event_loc_all,2)
    plot(event_loc_all(1,e),event_loc_all(2,e),'*b')
end

set(gca,'fontsize',20);
        