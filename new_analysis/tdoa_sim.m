function[tdoa_sim] = tdoa_sim(xloc,yloc,c,xpos,ypos)
% Returns the simulate tdoa_at given an event location, propagation speed,
% and the xy pos of the geophone receivers.

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

end