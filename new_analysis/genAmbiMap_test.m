function [ambi_map] = genAmbiMap_test(gridlist, x_mat, y_mat, SNR_list, ang_list)

   ambi_map = zeros(length(gridlist));
   figure
   xlim([gridlist(1) gridlist(end)])
   ylim([gridlist(1) gridlist(end)])
   theta = linspace(0,2*pi,100);
   r = 0.1:0.5:500;

   for n = 1:size(x_mat,1)
       n
       % initialize interp_grid matrix
       interp_grid = zeros(length(gridlist));
       
       for ind = 1:length(r)
           [y_int,x_int] = intersections(x_mat(n,:),y_mat(n,:),r(ind)*cos(theta),r(ind)*sin(theta),true);
           if ~isempty(x_int) && length(x_int)<5
               for ent = 1:length(x_int)
                   [~,xind] = min(abs(gridlist-x_int));
                   [~,yind] = min(abs(gridlist-y_int));
                   interp_grid = max(interp_grid, customGauss([length(gridlist) length(gridlist)], 10, 10, 0, 0, 1, [xind-floor(length(gridlist)/2) yind-floor(length(gridlist)/2)]));
               end
               
           else
               continue
           end
           
       end
       
       ambi_map = ambi_map + interp_grid;
       imagesc(gridlist,gridlist,ambi_map)
       set(gca,'YDir','normal')
       colormap jet
       pause;
       clf;
   end
   
     

end