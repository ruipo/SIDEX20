

function [xp,yp] = draw_circle(x,y,r,N)
% Function that returns x,y coords of a circle give the center and a radius

% x and y are the coordinates of the center of the circle
% r is the radius of the circle
% N is the number of points to use
% 0.01 is the angle step, bigger values will draw the circle faster but you might notice imperfections (not very smooth)

ang=linspace(0,2*pi,N); 
xp=x+r*cos(ang);
yp=y+r*sin(ang);

end