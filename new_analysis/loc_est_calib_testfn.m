%-----------------------------%
% Written by: Rui Chen
% Email: ruic@mit.edu
% Last updated: 02/28/2020
%-----------------------------%

function [loc_est,c_est,err,tdoa_mat] = loc_est_calib_testfn(zdata,xdata,ydata,xpos,ypos,start_sample,end_sample,FS,clist,N,plotting,calib_act)
% Estimates the location of events recorded on the z-axis channels. 
% Input data from x/y axis as zdata if want to estimate location of event on those axes. 

%-----------------------------------------------------------------------------------------------------------------------------------------------%
% INPUTS
% zdata - zaxis data [timeseries, channel];
% xdata - xaxis data [timeseries, channel];
% ydata - yaxis data [timeseries, channel];
% xpos - xpos of the geophones ex:[x1; x2; x3];
% ypos - ypos of the geophones ex:[y1; y2; y3];
% start_sample - start sample number of the event (int);
% end_sample - end sample number of the event (int);
% FS - sampling frequency in Hz (int);
% clist - list of propagation speeds to test;
% N - number of points generated to plot hyperbola (int) (Higher N -> more accurate prediction, runs more slowly, typically set to 500 or 1000)
% plotting - 1 (show intermediate plots); 0 (do not show intermediate plots);
% calib_act - actual location of calibration events [x y];

% OUTPUTS
% loc_est - estimated event location [x y];
% c_est - estimated propagation speed;
% err - error compared to simulated source at estimated event location;
% tdoa_mat - time difference of arrival matrix of input data;
%------------------------------------------------------------------------------------------------------------------------------------------------%

warning('off','all')  

% MPD using z,x-y axis data
mpd_mat = zeros(500,2,size(zdata,1));
for chn = 1:size(zdata,1)
    datain = [zdata(chn,start_sample:end_sample); xdata(chn,start_sample:end_sample); ydata(chn,start_sample:end_sample)];
    [~, ~, X_est, Y_est] = MPD(datain.',500);

    mpd_mat(:,1,chn) = X_est+xpos(chn);
    mpd_mat(:,2,chn) = Y_est+ypos(chn);
end

%get actual locations of the calibration events from the input calibration data
x_act = calib_act(1);
y_act = calib_act(2);

error = zeros(length(clist),1); %initialize errorlist
c_estlist = zeros(length(clist),1); %initialize c_estlist
loc_estlist = zeros(2,length(clist)); %initialize loc_estlist

for ccount = 1:length(clist) %loop throught list of propagation speeds
    
    if mod(ccount,100) == 0 
        ccount % COMMENT OUT IF YOU DO NOT WANT TO SEE CCOUNT NUMBER DURING RUN;
    end
    
    c0 = clist(ccount); %set propagation speed
    
    [tdoa_mat] = tdoa_sidex(zdata,start_sample,end_sample,FS,'hilbert'); %calculate tdoa matrix based on the propagation speed. 
    
    %------------------------------------------------------------------------------------------------------------------------------------------------%
    %find midpoints, focal radii, and rotation angles of all receiver pairs
    x_mat = [];
    y_mat = [];
    
    for i = 1:length(xpos)-1
        for j = i+1:length(xpos) %for reach receiver pair

            p1 = [xpos(i),ypos(i)]; %position of 1st receiver
            p2 = [xpos(j),ypos(j)]; %position of 2nd receiver

            m = (p1(:) + p2(:))/2; %midpoint
            c = sqrt((m(1)-p1(1))^2 + (m(2)-p1(2))^2); %distance from midpoint to a focus
            P = c0*tdoa_mat(i,j); %focal radii

            if abs(p2(1)-p1(1)) >= abs(p2(2)-p1(2)) %check orientation of hyperbola
                orientation = 1; %horizontal
                disp(orientation)
                loc_pre_rotation = [c 0]; %location of right receiver before rotation

                if p1(1)-m(1) >= 0 %check to see which receiver is the right receiver
                    loc_pre_translation = [(p1(1)-m(1)) (p1(2)-m(2))];
                elseif p2(1)-m(1) >= 0
                    loc_pre_translation = [(p2(1)-m(1)) (p2(2)-m(2))];
                end

                if loc_pre_translation(2) <= 0 %check to see what direction was the rotation and calculate the angle. (clockwise = positve ang)
                    ang = -acos(dot(loc_pre_rotation,loc_pre_translation)/(norm(loc_pre_rotation)*norm(loc_pre_translation)));
                elseif loc_pre_translation(2) > 0
                    ang = acos(dot(loc_pre_rotation,loc_pre_translation)/(norm(loc_pre_rotation)*norm(loc_pre_translation)));
                end

                a = abs(P)/2; %calculate a and b parameters for hyperbola.
                b = sqrt(c^2-a^2);

                points = create_hyperbola(a,b,ang,m(1),m(2),orientation,N); %create horizontal hyperbola.
                
                %determine which of the two hyperbolas to plot based on what sensor the event first arrived. 
                if P > 0 && p2(1)-m(1) >= 0
                    x = points(1,:);
                    y = points(2,:);
                elseif P > 0 && p2(1)-m(1) < 0
                    x = points(3,:);
                    y = points(4,:);
                elseif P < 0 && p2(1)-m(1) >= 0
                    x = points(3,:);
                    y = points(4,:);
                elseif P < 0 && p2(1)-m(1) < 0
                    x = points(1,:);
                    y = points(2,:);
                else
                    x = points(1,:);
                    y = points(2,:);
                end

                if plotting ~= 0 %plot hyperbola
                    figure(gcf)
                    plot(p1(1),p1(2),'ko');
                    hold on
                    plot(p2(1),p2(2),'ko');
                    plot(x,y,'r');
                    xL = xlim;
                    yL = ylim;
                    line([0 0], yL,'color','black');
                    line(xL, [0 0],'color','black');
                    pause
                    grid on
                end

            %------------------------------------------------------------------------------------------------------------------------------------------------%
            
            %repeat previous part for vertical hyperbola
            elseif abs(p2(1)-p1(1)) < abs(p2(2)-p1(2)) 
                orientation = 2; 
                disp(orientation)
                loc_pre_rotation = [0 c];

                if p1(2)-m(2) >= 0
                    loc_pre_translation = [(p1(1)-m(1)) (p1(2)-m(2))];
                elseif p2(2)-m(2) >= 0
                    loc_pre_translation = [(p2(1)-m(1)) (p2(2)-m(2))];
                end

                if loc_pre_translation(1) <= 0
                    ang = acos(dot(loc_pre_rotation,loc_pre_translation)/(norm(loc_pre_rotation)*norm(loc_pre_translation)));
                elseif loc_pre_translation(2) > 0
                    ang = -acos(dot(loc_pre_rotation,loc_pre_translation)/(norm(loc_pre_rotation)*norm(loc_pre_translation)));
                end

                a = abs(P)/2;
                b = sqrt(c^2-a^2);

                points = create_hyperbola(a,b,ang,m(1),m(2),orientation,N);
                
                if P > 0 && p2(2)-m(2) >= 0
                    x = points(1,:);
                    y = points(2,:);
                elseif P > 0 && p2(2)-m(2) < 0
                    x = points(3,:);
                    y = points(4,:);
                elseif P < 0 && p2(2)-m(2) >= 0
                    x = points(3,:);
                    y = points(4,:);
                elseif P < 0 && p2(2)-m(2) < 0
                    x = points(1,:);
                    y = points(2,:);
                else
                    x = points(1,:);
                    y = points(2,:);
                end

                if plotting ~= 0
                    figure(gcf);
                    plot(p1(1),p1(2),'ko');
                    hold on
                    plot(p2(1),p2(2),'ko');
                    plot(x,y,'r');
                    xL = xlim;
                    yL = ylim;
                    line([0 0], yL,'color','black');
                    line(xL, [0 0],'color','black');
                    pause
                    grid on
                end

            end
            
           %add all plotted hyperbola coordinates into a matrix 
           x_mat = [x_mat;x]; 
           y_mat = [y_mat;y];

        end
    end
    
    %-------------------------------------------------------------------------------------------------------------------------------------------------%
    % add in MPD lines
    
    for chn = 1:size(zdata,1)
        if plotting ~= 0
            figure(gcf);
            plot(mpd_mat(:,1,chn),mpd_mat(:,2,chn),'b');
            pause
        end
        x_mat = [x_mat; mpd_mat(:,1,chn).']; 
        y_mat = [y_mat; mpd_mat(:,2,chn).'];
    end
    

   %-------------------------------------------------------------------------------------------------------------------------------------------------%
   gs = -350;
   ge = 350;
   dg = 0.1;
   gridlist = gs:dg:ge;
   %gridlist2 = -350:1:350;
   sigma = 50;
   
   prob_grid = ones(length(gridlist));
   figure
   xlim([-350 350])
   ylim([-350 350])
   
    
   for n = 1:size(x_mat,1)
   
       interp_grid = zeros(length(gridlist));
       
       if issorted(x_mat(n,:)) || issorted(fliplr(x_mat(n,:)))
           y_interp = interp1(x_mat(n,:), y_mat(n,:), gridlist, 'spline');
           check1 = gs<y_interp;
           check2 = y_interp<ge;
           check = check1.*check2;
           xgridlist = gridlist(logical(check));
           ygridlist = y_interp(logical(check));
           for ent = 1:length(xgridlist)
               [~,ent1] = min(abs(gridlist-xgridlist(ent)));
               [~,ent2] = min(abs(gridlist-ygridlist(ent)));
               interp_grid(:,ent1) = 1/(sigma*sqrt(2*pi))*exp(-0.5*((([1:1:length(gridlist)]-ent2)./sigma).^2));
           end
          
       elseif issorted(y_mat(n,:)) || issorted(fliplr(y_mat(n,:)))
           x_interp = interp1(y_mat(n,:), x_mat(n,:), gridlist, 'spline');
           check1 = gs<x_interp;
           check2 = x_interp<ge;
           check = check1.*check2;
           ygridlist = gridlist(logical(check));
           xgridlist = x_interp(logical(check));
           for ent = 1:length(ygridlist)
               [~,ent1] = min(abs(gridlist-ygridlist(ent)));
               [~,ent2] = min(abs(gridlist-xgridlist(ent)));
               interp_grid(ent1,:) = 1/(sigma*sqrt(2*pi))*exp(-0.5*((([1:1:length(gridlist)]-ent2)./sigma).^2));
           end
           
       else
           check1 = gs<x_mat(n,:);
           check2 = x_mat(n,:)<ge;
           check3 = gs<y_mat(n,:);
           check4 = y_mat(n,:)<ge;
           check = check1.*check2.*check3.*check4;
           
           xregion = x_mat(n,:);
           xregion = xregion(logical(check));
           xregion2 = movmean(xregion,2);
           yregion = y_mat(n,:);
           yregion = yregion(logical(check));
           yregion2 = movmean(yregion,2);
           
           xregions = [xregion xregion2];
           yregions = [yregion yregion2];
           
           for ent = 1:length(xregions)
               [~,ent1] = min(abs(gridlist-xregions(ent)));
               [~,ent2] = min(abs(gridlist-yregions(ent)));
               interp_grid(:,ent1) = (interp_grid(:,ent1).' + 1/(sigma*sqrt(2*pi))*exp(-0.5*((([1:1:length(gridlist)]-ent2)./sigma).^2)))./2;             
           end
           
       end
       
       interp_grid = interp_grid./sum(interp_grid(:));
       prob_grid = prob_grid.*interp_grid;
       %prob_grid = prob_grid./sum(prob_grid(:));
       imagesc(gridlist,gridlist,prob_grid*100)
       set(gca,'YDir','normal')
       colormap jet
       pause;
       %clf;
       
   end
  


   %------------------------------------------------------------------------------------------------------------------------------------------------%
   %calculate the intersections of all hyperbolas and circles
    x_intlist = [];
    y_intlist = [];
    
    %for each receiver pair
    for n = 1:size(x_mat,1)-1
        for k = n+1:size(x_mat,1) 

            [x_int,y_int] = intersections(x_mat(n,:),y_mat(n,:),x_mat(k,:),y_mat(k,:),true); %find the intersections of each hyperbola with each other hyperbola.
            
            if length(x_int) > 10
                x_int = [];
                y_int = [];
            end
            
            %record the coordinates of all intersection points
            x_intlist = [x_intlist;x_int]; 
            y_intlist = [y_intlist;y_int];
            
        end
    end
    
    %------------------------------------------------------------------------------------------------------------------------------------------------%
    %estimate event location based on the intersection locations
    if ~isempty(x_intlist) && ~isempty(y_intlist)
        
        figure
        x_hist = histfit(real(x_intlist((x_intlist<350 & x_intlist>-350))),round(length(x_intlist)/2),'kernel'); %fit a kernel distribution to histogram of x-intersection values
        [xpeaks,xloclist] = findpeaks(x_hist(2).YData); %find location of the peaks of the kernel distribution
        [~,lx] = max(xpeaks); %x-estimate is the location of tallest peak
        locx = xloclist(lx);
        
        figure
        y_hist = histfit(real(y_intlist((y_intlist<350 & y_intlist>-350))),round(length(y_intlist)/2),'kernel'); %fit a kernel distribution to histogram of y-intersection values
        [ypeaks,yloclist] = findpeaks(y_hist(2).YData); %find location of the peaks of the kernel distribution
        [~,ly] = max(ypeaks);%y-estimate is the location of the tallest peak
        locy = yloclist(ly);
            
        loc_est = [x_hist(2).XData(locx); y_hist(2).XData(locy)];
            
        close
        close
    
    else
        
        loc_est = [nan;nan];
    end
    
    
    if isempty(xloclist) %if no peaks in kernel function is found, use mode method for x,y estimates
            loc_est = [mode(round(x_intlist,3)) loc_est];
    end
        
    if isempty(yloclist)
            loc_est = [loc_est mode(round(y_intlist,1))];
    end
    
    %------------------------------------------------------------------------------------------------------------------------------------------------%
    %calculates the error between estimation and actual event locations
    %selects the propagation speed/estimation location with the lowest error
    
    error(ccount) = sqrt((x_act - loc_est(1))^2 + (y_act - loc_est(2))^2); 
    c_estlist(ccount) = c0;
    loc_estlist(1,ccount) = loc_est(1);
    loc_estlist(2,ccount) = loc_est(2);

    if plotting ~= 0
        figure(gcf);
        plot(loc_est(1),loc_est(2),'r*');
        hold off

        xlabel('X position (m)')
        xlabel('Y position (m)')
        title(['propagation speed = ', num2str(c0), ' m/s']);
        xlim([-350 350]);
        ylim([-150 350]);
        
        disp(['x estimate = ', num2str(loc_est(1)), '; y estimate = ',num2str(loc_est(2)) '.']);
        disp(['propagation speend estimate = ', num2str(c0), 'm/s.']);
        disp(['estimate error = ', num2str(error(ccount)), '.']);
        
        pause
    end
% STOP EARLY TO SAVE COMPUTATION TIME    
    if error(ccount) < 0.1
        break
    end
    
end

%find the minimum error
%find corresponding propagation speed and location estimate
[err,eloc] = min(error(1:ccount));
c_est = c_estlist(eloc); 
loc_est = loc_estlist(:,eloc);

warning('on','all')

end

