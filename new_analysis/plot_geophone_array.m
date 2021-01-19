%% map geophone and hydrophon locations
load('/Users/Rui/Documents/Graduate/Research/SIDEX/SIDEX20/new_analysis/2020_locations.mat')

[crack_long,ind] = sort(crack_long);
crack_lat = crack_lat(ind);

%geophone_GPS = [71.332947,-156.407381; 71.33338,-156.409125; 71.333504,-156.407474; 71.333275,-156.408282];
%standalone_GPS =[71.3357 -156.3982; 71.3299 -156.4016; 71.3345 -156.4165; 71.3333 -156.4081];

figure
%worldmap([71.325 71.34],[-156.42 -156.39])
plot(geophone_GPS(:,2),geophone_GPS(:,1),'b.', 'MarkerSize',30);
hold on
%h = scatterm(itp_lat(1756:2017),itp_lon(1756:2017),86,'c','d','filled');
%h = scatterm(itp_lat(2018:2064),itp_lon(2018:2064),86,'b','^','filled');
%h = scatterm(itp_lat(2065:2102),itp_lon(2065:2102),66,'g','o','filled');
%h = scatterm(itp_lat(2103:2122),itp_lon(2103:2122),88,color1,'*');
plot(standalone_GPS(:,2),standalone_GPS(:,1),'b^', 'MarkerSize',8,'MarkerFaceColor','b');
plot(geophone_GPS(3,2),geophone_GPS(3,1),'rx', 'MarkerSize',8,'linewidth',1.5);
plot(crack_long,crack_lat,'k','linewidth',1)
%scatterm(73.25,-149.4,555,'r','p','filled');
grid on
xlim([-156.42 -156.395])
ylim([71.329 71.337])
legend('Cabled Geophones','Stand-alone Geohpones','Hydrophones')
set(gca,'fontsize',20)
xlabel('Longitude')
ylabel('Latitude')
%scaleruler
%set(gca,'Visible','off')
lat_lon_proportions()

%% Alaska

load coastlines

figure
worldmap([55 75.336],[-175.42 -135.395])
geoshow(coastlat,coastlon,'Color','k')
hold on
scatterm(geophone_GPS(3,1),geophone_GPS(3,2),88,'r','p','filled');
set(gca,'fontsize',20)

%% Calibration events locations

% run sidex20_calib_run to get lat long of calibration events
hold on
plot(long_event,lat_event,'r*','MarkerSize',8)

