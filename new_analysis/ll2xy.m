function [x,y] = ll2xy(lat,long,x0,y0)

% converts lat/long coordinates to x/y given the lat/long that corresponds to [0,0]

% INPUTS
% lat - lat coordinate to be converted
% long - long coordinate to be converted
% x0 - lat coordinate that corresponds to x=0 in xy plane
% y0 - long coordinate that corresponds to y=0 in xy plane

% OUTPUTS
% x - x coordinate of input lat/long
% y - y coordinate of input lat/long


x = zeros(length(lat),1);
y = x;


for i = 1:length(lat)
    
    if lat(i) < x0
        y(i) = -lldistkm([lat(i) y0],[x0 y0])*1000;
    else
        y(i) = lldistkm([lat(i) y0],[x0 y0])*1000;
    end
    
    if long(i) < y0
        x(i) = -lldistkm([x0 long(i)],[x0 y0])*1000;
    else
        x(i) = lldistkm([x0 long(i)],[x0 y0])*1000;
    end
end


end