% Made for SIDEX data collected in 2018-2019
% Requires the data directory structure of ~/sidex/sidex_DATE/sidex_TIME.txt
% with required code in ~/sidex/code/ folder

% RUN code in ~/sidex/ folder
% can be modified to analyzed similar datasets

disp('WARNING WILL OVERWRITE EXISTING FOLDERS IF ANALYSIS IS ALREADY DONE!!! PROCEED?') % run warning, press enter to continue

pause

folder = dir('~/sidex/sidex*');
for fold = 1:length(folder)
    
    timestamp = [];
    
    % set directory
    prefix = ['~/sidex/' folder(fold).name '/'];
    directory = dir([prefix 'Sidex_201*.txt']);
    date = directory(end).name(7:14); % records the date that the data was collected
    FS = 500;
    num_chn = 16;
    good_chn = [1:6 10:12];%choose chns to be plotted
    first_file = 1; %choose the first file to be read

    foldername = [prefix(end-10:end-1) '_stats'];
    num_event = []; % record number of event detections in each axis
    num_event_zaxis = [];
    num_event_xyaxis = [];
    num_event_allaxis = [];

    cycle = 1;
    while first_file < length(directory)
        data = zeros(1,16);
        datenumlist = zeros(1000,1);
        t = [];

        for i = first_file:length(directory)% for each data file in directory
            time = directory(i).name(16:21);% record time of file
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

        first_file = i+1;
        
        data(1,:) = [];
        for chn = 1:size(data,2) % demean and normalize data
            data(:,chn) = data(:,chn)-mean(data(:,chn));
            data(:,chn) = data(:,chn)./std(data(:,chn));
        end
        
        %separate data into x, y, z axes
        xdata = data(:,[1 4 10]).'; 
        ydata = data(:,[2 5 11]).';
        zdata = data(:,[3 6 12]).';
        
        cd('code')
       
        % Event detection

        lta_size_s = 1; %LTA window in seconds
        sta_size_s = 0.05; %STA window in seconds
        pem_s = 0.05;
        pet_s = 0.05;
        ratio_thres = 1.5;
        weight_thres = 1;

        % x-dir detection
        xnum_chn = size(xdata,1);
        [xstart_sample,xend_sample,xchn_triggered_mat,~,~,~] = sl_eSelect(xdata,FS,lta_size_s,sta_size_s,pem_s,pet_s,ratio_thres,weight_thres);

        % y-dir detection
        ynum_chn = size(ydata,1);
        [ystart_sample,yend_sample,ychn_triggered_mat,~,~,~] = sl_eSelect(ydata,FS,lta_size_s,sta_size_s,pem_s,pet_s,ratio_thres,weight_thres);

        % z-dir detection
        znum_chn = size(zdata,1);
        [zstart_sample,zend_sample,zchn_triggered_mat,~,~,~] = sl_eSelect(zdata,FS,lta_size_s,sta_size_s,pem_s,pet_s,ratio_thres,weight_thres+1);

        cd ..
        
        % Record each event, count number of axes that they are seen on
        if length(xend_sample)~=length(xstart_sample)
        xend_sample = xend_sample(1:min([length(xstart_sample),length(xend_sample)]));
        xstart_sample = xstart_sample(1:min([length(xstart_sample),length(xend_sample)]));
        end
        
        if length(yend_sample)~=length(ystart_sample)
        yend_sample = yend_sample(1:min([length(ystart_sample),length(yend_sample)]));
        ystart_sample = ystart_sample(1:min([length(ystart_sample),length(yend_sample)]));
        end
        
        if length(zend_sample)~=length(zstart_sample)
        zend_sample = zend_sample(1:min([length(zstart_sample),length(zend_sample)]));
        zstart_sample = zstart_sample(1:min([length(zstart_sample),length(zend_sample)]));
        end
        
        u = union(xstart_sample,ystart_sample);
        u2 = union(u,zstart_sample);
        
        for ent = fliplr(2:length(u2))
            if u2(ent)-u2(ent-1) <= 2*sta_size_s*FS
                u2(ent) = [];
            end
        end
            
        
        m_num = length(u2);
        tally_start = zeros(3,m_num);
        tally_end = zeros(3,m_num);
        
        for entry = 1:m_num
            
            x_s = abs(xstart_sample - u2(entry));
            xolap = x_s <= 2*sta_size_s*FS;
            indx = find(xolap==1);

            y_s = abs(ystart_sample - u2(entry));
            yolap = y_s <= 2*sta_size_s*FS;
            indy = find(yolap==1);

            z_s = abs(zstart_sample - u2(entry));
            zolap = z_s <= 2*sta_size_s*FS;
            indz = find(zolap==1);

            if isempty(indx)
                tally_start(1,entry) = NaN;
                tally_end(1,entry) = NaN;

            else 
                tally_start(1,entry) = xstart_sample(min(indx));
                tally_end(1,entry) = xend_sample(min(indx));
            end 

            if isempty(indy)
                tally_start(2,entry) = NaN;
                tally_end(2,entry) = NaN;

            else 
                tally_start(2,entry) = ystart_sample(min(indy));
                tally_end(2,entry) = yend_sample(min(indy));
            end 


            if isempty(indz)
                tally_start(3,entry) = NaN;
                tally_end(3,entry) = NaN;

            else 
                tally_start(3,entry) = zstart_sample(min(indz));
                tally_end(3,entry) = zend_sample(min(indz));
            end 
        end 
        
        t_start = nan(size(tally_start));
        t_end = nan(size(tally_end));
        
        for i = 1:3
            t_start(i,~isnan(tally_start(i,:))) = t(tally_start(i,~isnan(tally_start(i,:))));
            t_end(i,~isnan(tally_end(i,:))) = t(tally_end(i,~isnan(tally_end(i,:))));
        end
       
        tallyxy = tally_start(1:2,:);
        xyaxis = sum(~isnan(tallyxy),1);
        xyaxis = find(xyaxis==2);
        num_xyaxis = length(xyaxis);
        
        allaxis = sum(~isnan(tally_start),1);
        allaxis = find(allaxis==3);
        num_allaxis = length(allaxis);
        
        zaxis = ~isnan(tally_start(3,:));
        zaxis = find(zaxis==1);
        num_zaxis = length(zaxis);
        
        cd(foldername)
        name = ['hour_',num2str(cycle),'_eSelect'];
        mkdir(name)
        cd ..
        
        % save results to folder
        save(['~/sidex/' foldername '/' name '/e_Select'],'xstart_sample','xend_sample','ystart_sample','yend_sample','zstart_sample','zend_sample','tally_start','tally_end','t_start','t_end')
        
        num_event(cycle) = size(tally_start,2);
        num_event_allaxis(cycle) = num_allaxis;
        num_event_zaxis(cycle) = num_zaxis;
        num_event_xyaxis(cycle) = num_xyaxis;
        timestamp(cycle) = t(round(length(t)/2));
        
        cycle = cycle + 1;
        
    end
    
    eventtimelist = datetime(timestamp,'ConvertFrom','datenum');
    
    % count event totals in this folder (daily)
    save(['~/sidex/' foldername '/num_events'],'num_event','num_event_allaxis','num_event_zaxis','num_event_xyaxis','timestamp','eventtimelist') 
        
end
        
disp('DONE!')

% CAN NOW RUN sidex_tdoa_hourly_run.m to localized detected events!

%% Compile all results

num_event_all = [];
num_event_allaxis_all = [];
num_event_xyaxis_all = [];
num_event_zaxis_all = [];
eventtimelist_all = [];
timestamp_all = [];

folder = dir('~/sidex/*stats');
for fold = 2:length(folder)
    fold
    
    cd(folder(fold).name)
    load('num_events.mat')
    
    num_event_all = [num_event_all num_event];
    num_event_allaxis_all = [num_event_allaxis_all num_event_allaxis];
    num_event_xyaxis_all = [num_event_xyaxis_all num_event_xyaxis];
    num_event_zaxis_all = [num_event_zaxis_all num_event_zaxis];
    eventtimelist_all = [eventtimelist_all eventtimelist];
    timestamp_all = [timestamp_all timestamp];
    
    cd('~/sidex/')
end

%% Plotting

figure
subplot(3,1,1)
stem(eventtimelist_all,num_event_allaxis_all)
%xlabel('Date')
%ylabel('Number of Detections')
title('Detection on all axes')
grid on
set(gca,'fontsize',20)
set(gca, 'YScale', 'log')
subplot(3,1,2)
stem(eventtimelist_all,num_event_xyaxis_all)
%xlabel('Date')
ylabel('Number of Detections')
title('Detection on x-y axes')
grid on
set(gca,'fontsize',20)
set(gca, 'YScale', 'log')
subplot(3,1,3)
stem(eventtimelist_all,num_event_zaxis_all)
xlabel('Date')
%ylabel('Number of Detections')
title('Detection on z axis')
grid on
set(gca,'fontsize',20)
set(gca, 'YScale', 'log')



    