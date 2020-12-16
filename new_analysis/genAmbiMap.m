function [x_est,y_est,ambi_map] = genAmbiMap(gridlist, x_mat, y_mat, SNR_list, ang_list)
     
   ambi_map = zeros(length(gridlist));
   figure
   xlim([gridlist(1) gridlist(end)])
   ylim([gridlist(1) gridlist(end)])
    
   for n = 1:size(x_mat,1)
       n
       % determine sigma of Gaussian using SNR
%        if SNR_list(n) < 20
%            sigma = 15-SNR_list(n);
%        else
%            sigma = 5;
%        end
       sigma = -50/(1+20*exp(-0.16*SNR_list(n)))+53;
       sigma = round(sigma);
       disp(sigma)
       if sigma >= 45
           continue
       end
         
       % initialize interp_grid matrix
       interp_grid = zeros(length(gridlist));
       
       % Hyperbola has no repeated x values
       if issorted(x_mat(n,:)) || issorted(fliplr(x_mat(n,:)))
           y_interp = interp1(x_mat(n,:), y_mat(n,:), gridlist, 'spline'); % interpolate x,y to gridlist and y_interp
           % segment interpolated values to within bounds of gridlist
           check1 = gridlist(1)<y_interp;
           check2 = y_interp<gridlist(end);
           check = check1.*check2;
           xgridlist = gridlist(logical(check));
           ygridlist = y_interp(logical(check));
           % fill in interp_grid matrix
           for ent = 1:length(xgridlist)
               [~,ent1] = min(abs(gridlist-xgridlist(ent)));
               [~,ent2] = min(abs(gridlist-ygridlist(ent)));
               %interp_grid(:,ent1) = 1/(sigma*sqrt(2*pi))*exp(-0.5*((([1:1:length(gridlist)]-ent2)./sigma).^2));
               interp_grid = max(interp_grid, customGauss([length(gridlist) length(gridlist)], sigma, sigma, 0, 0, -1/(1+10*exp(-0.077*sigma))+1, [ent2-floor(length(gridlist)/2) ent1-floor(length(gridlist)/2)]));

           end
           
       % Hyperbola has no repeated y values  
       elseif issorted(y_mat(n,:)) || issorted(fliplr(y_mat(n,:)))
           x_interp = interp1(y_mat(n,:), x_mat(n,:), gridlist, 'spline');
           check1 = gridlist(1)<x_interp;
           check2 = x_interp<gridlist(end);
           check = check1.*check2;
           ygridlist = gridlist(logical(check));
           xgridlist = x_interp(logical(check));
           for ent = 1:length(ygridlist)
               [~,ent1] = min(abs(gridlist-ygridlist(ent)));
               [~,ent2] = min(abs(gridlist-xgridlist(ent)));
               %interp_grid(ent1,:) = 1/(sigma*sqrt(2*pi))*exp(-0.5*((([1:1:length(gridlist)]-ent2)./sigma).^2));
               interp_grid = max(interp_grid, customGauss([length(gridlist) length(gridlist)], sigma, sigma, 0, 0, -1/(1+10*exp(-0.077*sigma))+1, [ent1-floor(length(gridlist)/2) ent2-floor(length(gridlist)/2)]));


           end
       
       % Hyperbola has repeated x or y values    
       else
           % rotate hyperbola to new coordinate system so there are no repreated x or y values
           xrot = x_mat(n,:).*cos(rad2deg(ang_list(n)))-y_mat(n,:).*sin(rad2deg(ang_list(n)));
           yrot = x_mat(n,:).*sin(rad2deg(ang_list(n)))+y_mat(n,:).*cos(rad2deg(ang_list(n)));
           gridlistrot = gridlist(1):(gridlist(2)-gridlist(1))*0.8:gridlist(end); % use a slightly finer discritization step so when we undo the rotation, we get a fine enough grid
           
           % if rotated hyperbola has no repeated x values
           if issorted(xrot) || issorted(fliplr(xrot))
               y_interprot = interp1(xrot, yrot, gridlistrot, 'spline'); % interpolate
               % rotate back to original coordinate system
               gridlistbrot = gridlistrot.*cos(rad2deg(-ang_list(n)))-y_interprot.*sin(rad2deg(-ang_list(n)));
               y_interp = gridlistrot.*sin(rad2deg(-ang_list(n)))+y_interprot.*cos(rad2deg(-ang_list(n)));
               % clip to within bounds of gridlist
               check1 = gridlist(1)<y_interp;
               check2 = y_interp<gridlist(end);
               check3 = gridlist(1)<gridlistbrot;
               check4 = gridlistbrot<gridlist(end);
               check = check1.*check2.*check3.*check4;
               xgridlist = gridlistbrot(logical(check));
               ygridlist = y_interp(logical(check));
               % fill in interp_grid matrix
               for ent = 1:length(xgridlist)
                   [~,ent1] = min(abs(gridlist-xgridlist(ent)));
                   [~,ent2] = min(abs(gridlist-ygridlist(ent)));
                   % use max of interp_grid and new values becuase sometimes the new values close to zero will over write older values otherwise
                   %interp_grid(:,ent1) = max(interp_grid(:,ent1).', 1/(sigma*sqrt(2*pi))*exp(-0.5*((([1:1:length(gridlist)]-ent2)./sigma).^2)));
                   interp_grid = max(interp_grid, customGauss([length(gridlist) length(gridlist)], sigma, sigma, 0, 0, -1/(1+10*exp(-0.077*sigma))+1, [ent2-floor(length(gridlist)/2) ent1-floor(length(gridlist)/2)]));


               end
               
           % if rotated hyperbola has no repeated y values
           elseif issorted(yrot) || issorted(fliplr(yrot))
               x_interprot = interp1(yrot, xrot, gridlistrot, 'spline');
               x_interp = x_interprot.*cos(rad2deg(-ang_list(n)))-gridlistrot.*sin(rad2deg(-ang_list(n)));
               gridlistbrot = x_interprot.*sin(rad2deg(-ang_list(n)))+gridlistrot.*cos(rad2deg(-ang_list(n)));
               check1 = gridlist(1)<x_interp;
               check2 = x_interp<gridlist(end);
               check3 = gridlist(1)<gridlistbrot;
               check4 = gridlistbrot<gridlist(end);
               check = check1.*check2.*check3.*check4;
               ygridlist = gridlistbrot(logical(check));
               xgridlist = x_interp(logical(check));
               for ent = 1:length(ygridlist)
                   [~,ent1] = min(abs(gridlist-ygridlist(ent)));
                   [~,ent2] = min(abs(gridlist-xgridlist(ent)));
                   %interp_grid(ent1,:) = max(interp_grid(ent1,:), 1/(sigma*sqrt(2*pi))*exp(-0.5*((([1:1:length(gridlist)]-ent2)./sigma).^2)));
                   interp_grid = max(interp_grid, customGauss([length(gridlist) length(gridlist)], sigma, sigma, 0, 0, -1/(1+10*exp(-0.077*sigma))+1, [ent1-floor(length(gridlist)/2) ent2-floor(length(gridlist)/2)]));


               end
           end
       end
       
       % Normalize interp_grid before adding to ambi_map
       %interp_grid = interp_grid./sum(interp_grid(:));
       ambi_map = ambi_map + interp_grid;
       %prob_grid = prob_grid./sum(prob_grid(:));

       %pause;
       %clf;
       
   end

   roi = ambi_map>max(max(ambi_map))*0.9;
   imagesc(gridlist,gridlist,10*log10(ambi_map))
   set(gca,'YDir','normal')
   colormap jet
   colorbar
   caxis([0 10]);
   hold on
   
   if max(any(roi))>0
       stats = regionprops(roi);
       centroid = stats.Centroid;
       bbox = stats.BoundingBox;
       bxmin = gridlist(ceil(bbox(1)));
       bymin = gridlist(ceil(bbox(2)));
       x_est = gridlist(round(centroid(1)));
       y_est = gridlist(round(centroid(2)));

       plot(x_est,y_est,'k*')
       rectangle('position',[bxmin bymin bbox(3)-1 bbox(4)-1])
   else
       x_est = NaN;
       y_est = NaN;
   end
       
end