
function[err_mat] = tdoa_cmp(xlist, ylist, clist,xpos, ypos, tdoa_mat)

err_mat = zeros(length(xlist),length(ylist),length(clist));
l = length(xpos);

for xx = 1:length(xlist)
    xx
    x = xlist(xx);
    for yy = 1:length(ylist)
        y = ylist(yy);
        for cc = 1:length(clist)
            c = clist(cc);
            tdoa_sim = zeros(l);
            
            % find distance and arrival time to each receiver
            d = sqrt((x-xpos).^2 + (y-ypos).^2);
            t = d./c;

            % find simulated tdoa between each receiver
            for j = 1:l
                for k = j:l
                    tdoa_sim(k,j) = (t(k)-t(j));
                    tdoa_sim(j,k) = -(t(k)-t(j));
                end
            end
            
            err_mat(xx,yy,cc) = norm(tdoa_mat-tdoa_sim);
        end
    end
end