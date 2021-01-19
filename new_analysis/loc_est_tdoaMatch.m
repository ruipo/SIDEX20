function [loc_est,c_est,err,rms_err_mat] = loc_est_tdoaMatch(zdata,xpos,ypos,start_sample,end_sample,FS,xrange,yrange,ds,clist,plotting,calib_act,ploti)

error = zeros(length(clist),1); %initialize errorlist
loc_estlist = zeros(2,length(clist)); %initialize loc_estlist

% get actual locations of the calibration events from the input calibration data
x_act = calib_act(1);
y_act = calib_act(2);

% initialize test space and rms error matrix
xlist = xrange(1):ds:xrange(end);
ylist = yrange(1):ds:yrange(end);
rms_err_mat = zeros(length(xlist),length(ylist));

for ccount = 1:length(clist) %loop throught list of propagation speeds
    
    if mod(ccount,10) == 0 
        ccount % COMMENT OUT IF YOU DO NOT WANT TO SEE CCOUNT NUMBER DURING RUN;
    end
    
    c0 = clist(ccount); %set propagation speed
    
    [tdoa_mat] = tdoa_sidex(zdata,start_sample,end_sample,FS,'hilbert'); %calculate tdoa matrix based on the propagation speed. 
    
    for xind = 1:length(xlist)
        for yind = 1:length(ylist)
            
            xtest = xlist(xind);
            ytest = ylist(yind);
            
            tdoa_test = tdoa_sim(xtest,ytest,c0,xpos,ypos);
            
            rms_err_mat(xind,yind) = sqrt(mean(mean((tdoa_mat-tdoa_test).^2)));
            
        end
    end
    
    Nk = sum(sum(1./rms_err_mat));
    prob_mat = (1./rms_err_mat)./(Nk);
    mk = mean(mean(prob_mat));
    
    if plotting ~= 0
        figure('units','normalized','outerposition',[0 0 1 1])
        imagesc(xlist,ylist,10*log10(prob_mat.'./mk));
        xL = xlim;
        yL = ylim;
        colormap jet
        colorbar
        hold on

        xlabel('X Position (m)')
        ylabel('Y Position (m)')
        xlim([-350 350]);
        ylim([-150 350]);
        set(gca,'YDir','normal');
        %pause

    end
    
    [xestind,yestind] = find(rms_err_mat == min(rms_err_mat(:)));
    loc_estlist(1,ccount) = xlist(xestind);
    loc_estlist(2,ccount) = ylist(yestind);
    error(ccount) = sqrt((x_act - xlist(xestind))^2 + (y_act - ylist(yestind))^2); 
    
    
end
[err,eloc] = min(error(1:ccount));
c_est = clist(eloc); 
loc_est = loc_estlist(:,eloc);

if plotting ~= 0
    figure(gcf)
    path = '/Users/Rui/Documents/Graduate/Research/SIDEX/SIDEX20/new_analysis/';        
    plot(x_act, y_act,'k*','MarkerSize',8,'linewidth',1.5);
    plot(loc_est(1), loc_est(2),'m*','MarkerSize',8,'linewidth',1.5);
    plot(xpos,ypos,'k.', 'MarkerSize',30);
    caxis([-2 10]);
    legend('True Location','Estimated Location','Cabled Geophones')
    set(gca,'fontsize',20);
    title(['Propagation Speed = ' num2str(c0) ' m/s']);
    saveas(gcf,[path,'loc_est_tdoaMatch/' num2str(ploti) '.png']);
    close all
end



end






