function[tdoa_sim] = tdoa_sim(xloc,yloc,c,xpos,ypos)
% Returns the simulate tdoa_at given an event location, propagation speed,
% and the xy pos of the geophone receivers.

%d = [];
%t = [];
tdoa_sim = zeros(length(xpos),length(ypos));
% find distance and arrival time to each receiver
%for i = 1:length(xpos)
    d = sqrt((xpos-xloc).^2 + (ypos-yloc).^2);
    t = d./c;
%end

% find simulated tdoa between each receiver
for j = 1:length(t)
    for k = j:length(t)
        tdoa_sim(k,j) = t(k)-t(j);
        tdoa_sim(j,k) = -tdoa_sim(k,j);
    end
end

end