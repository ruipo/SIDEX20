close all
path(path,'~/mural_gitlab/sidex_src/src/sidex_eselect_localization/')
 standalone_GPS =[
     71.3357 -156.3982
     71.3299 -156.4016
   71.3345 -156.4165
   71.3333 -156.4081]

SPHEROID = wgs84Ellipsoid('m');

lat_start = standalone_GPS(1,1);
long_start =standalone_GPS(1,2);

map_grid_x = []

% for each point in the grid, calculate x,y grid location:


% assume flatness for now, it makes life easier

[y,x]= geodetic2ned(standalone_GPS(:,1),standalone_GPS(:,2),0,lat_start,long_start,0,SPHEROID);

% import csv as individual coluns
FS=4060;
data=Untitled;
data(isnan(data))=0;
tvec = (1:length(Untitled(:,1)))/FS;
tmat=tvec*0;
offset_vec=[0,0,0,0]%[0,2,2,6]%[1,1,3,6];%[4,0,0,4];% offsets in s per channel
fft_length=1000;
figure(3)
hold on
figure(2)
hold on
scaling=0.03;
%%
elnum=[1,1,1,2,2,2,3,3,3,4,4,4];
t = 0:1/FS:size(data,1)/FS - 1/FS;
[b,a]=butter(6,100/(1000/2),'low');
[d,c]=butter(6,10/(1000/2),'high');
data_filt=nan*data;
for ii=1:length(data(1,:))
    offset_ind=offset_vec(elnum(ii))*FS;
    data_filt(1:end-offset_ind,ii)=data(offset_ind+1:end,ii);%filter(b,a,data(offset_ind+1:end,ii));
    
    data_filt(:,ii)=filter(b,a,data_filt(:,ii));
    data_filt(:,ii)=filter(d,c,data_filt(:,ii));
    figure(3)
    hold on
    plot((1:length(data_filt(:,ii)))/FS,data_filt(:,ii)/max(data_filt(:,ii))+ii)
    figure(2)
    
%     %pause
%     subplot(5,3,ii)
%     specgram(data_filt(:,ii),fft_length,FS,fft_length,fft_length-1)
%     title(num2str(ii))
%     caxis([-70,-30])
%     axis([0,6,0,200])
%     subplot(5,3,13:15)
%     hold on
%     plot((1:length(data(:,ii)))/FS,data_filt(:,ii)/max(data_filt(:,ii))+ii)
%     axis([0,6,0,14])
%     % filter
%     
end
pause
% TDOA:
% find each event segment:
%%
%figure
figure
 [PKS,LOCS] = findpeaks((data_filt(:,4)),'MinPeakHeight',1e-4,'MinPeakDistance',1000);
 findpeaks((data_filt(:,4)),'MinPeakHeight',1e-4,'MinPeakDistance',1000)


 standalone_GPS =[
     71.3357 -156.3982
     71.3299 -156.4016
   71.3345 -156.4165
   71.3333 -156.4081]
%  standalone_GPS =[
%    71.3345 -156.4165
%    71.3357 -156.3982
%    71.3299 -156.4016
%    71.3333 -156.4081]

geophone_GPS = standalone_GPS;
SPHEROID =wgs84Ellipsoid('m');
[geophone_y,geophone_x]=geodetic2ned(geophone_GPS(:,1),geophone_GPS(:,2),0,geophone_GPS(4,1),geophone_GPS(4,2),0,SPHEROID)
 geophone_x=geophone_x -30;% adjust to coordinates for other array
 geophone_y=geophone_y+30;
 figure(10)
hold on
plot(geophone_x,geophone_y,'b-')
clear loc_est*

clear diff*

for pkind=1:length(PKS)
    % event segment:
    st_ind = max(LOCS(pkind)-1000,1);
    e_ind = min(LOCS(pkind)+10000,size(data_filt(:,1)));
    event_data=data_filt(st_ind:e_ind,:);
    
%event_data_z =  data_filt(LOCS(pkind)-100:LOCS(pkind)+1000,:);
% z_data = [event_data_z(:,1),event_data_z(:,4),event_data_z(:,7),event_data_z(:,10)];
% x_data = [event_data(:,2),event_data(:,5),event_data(:,8),event_data(:,11)];
% y_data = [event_data(:,3),event_data(:,6),event_data(:,9),event_data(:,12)];
 y_data = [event_data(:,3),event_data(:,6),event_data(:,9),event_data(:,12)];
 z_data = [event_data(:,1),event_data(:,4),event_data(:,7),event_data(:,10)];
 x_data = [event_data(:,2),event_data(:,5),event_data(:,8),event_data(:,11)];

% adjust 
% next, run 
clists=850
ii=1;
jj=1;
c_z=150
jj=pkind
[loc_est_z(jj,1:2),c_est,err] = loc_est_hyp_constc(z_data',geophone_x,geophone_y,1,length(event_data),1000,c_z,500,0)
    
ii=pkind;
c_xy=100
[loc_est_x(ii,1:2),c_est,err] = loc_est_hyp_constc(x_data',geophone_x,geophone_y,1,length(event_data),1000,c_xy,500,0)

[loc_est_y(ii,1:2),c_est,err] = loc_est_hyp_constc(y_data',geophone_x,geophone_y,1,length(event_data),1000,c_xy,500,0)
    
    diff_xz(pkind) = pdist([real(loc_est_x(ii,:));real(loc_est_z(jj,:))],'euclidean');
    diff_xy(pkind)=pdist([real(loc_est_x(ii,:));loc_est_y(ii,:)],'euclidean');
    diff_yz(pkind)=pdist([real(loc_est_z(jj,:));loc_est_y(ii,:)],'euclidean');
    ii=ii+1;
    
    jj=jj+1;


%figure
%pcolor(200:10:410,200:10:400,diff_xz+diff_xy+diff_yz)




figure(10)
 clf
 hold on
 plot(geophone_x,geophone_y,'b-')
 plot(loc_est_z(pkind,1),loc_est_z(pkind,2),'r*')
 plot(loc_est_x(pkind,1),loc_est_x(pkind,2),'rx')
 plot(loc_est_y(pkind,1),loc_est_y(pkind,2),'r.')


%axis([-100,50,-50,100])
end

figure(10)
%geophone_x=geophone_x -30;% adjust to coordinates for other array
% geophone_y=geophone_y+30;
 figure(10)
 clf
hold on
plot(geophone_x,geophone_y,'b-')
plot(loc_est_z(:,1),loc_est_z(:,2),'r*')
plot(loc_est_x(:,1),loc_est_x(:,2),'rx')
plot(loc_est_y(:,1),loc_est_y(:,2),'r.')
%%
