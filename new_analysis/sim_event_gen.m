
function[tdoa_sim,data] = sim_event_gen(xloc,yloc,c,xpos,ypos,FS)
% Simulate data from trasient event that occurs at xloc, yloc, and received
% at geophone locations xpos, ypos at sampling frequency FS. Assume
% propagation speed of c.

d = [];
t = [];
tdoa_sim = [];
% find distance and arrival time to each receiver
for i = 1:length(xpos)
    d(i) = sqrt((xpos(i)-xloc)^2 + (ypos(i)-yloc)^2);
    t(i) = d(i)/c;
end

% find simulated tdoa between each receiver
for j = 1:length(t)
    for k = 1:length(t)
        tdoa_sim(k,j) = t(k)-t(j);
    end
end


% find maximum travel time between receivers
for j = 1:length(xpos)
    for k = 1:length(xpos)
        dlist(j,k) = sqrt((xpos(j)-xpos(k))^2 + (ypos(j)-ypos(k))^2);
        max_t = max(max(dlist))/c;
    end
end

% simulate time series; pretend event is 1 second long and total time
% series is 5 seconds long or 5 times as long as max_t if max_t > 1 sec
event_sim = randn(FS,1);

if length(event_sim) >= max_t*FS
    data = zeros(length(xpos),500*length(event_sim));
elseif length(event_sim) < max_t*FS
    data = zeros(length(xpos),10*round(max_t*FS));
end
dl = size(data,2);
time = 0:1/FS:dl/FS - 1/FS;

% add event into simulated data time series
for n = 1:size(tdoa_sim,1)
    data(n,round(dl/2)+round(tdoa_sim(1,n)*FS)-round(FS/2):round(dl/2)+round(tdoa_sim(1,n)*FS)+round(FS/2)-1) = (event_sim + 0.1*randn(FS,1));
end

%%%% plot simulated data

figure
for ch = 1:size(data,1)
plot(time,data(ch,:)+(ch-1)*max(max(abs(data))))
hold on
xlabel('time')
ylabel('amplitude')
xlim([time(1) time(end)]);
title(['SImulated data for event at (', num2str(xloc), ',', num2str(yloc), ')']);
set(gca,'fontsize',25);
grid on
end
end