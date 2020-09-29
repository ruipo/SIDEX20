
function [fout] = read_moveout_asc(filename)
    
    f = readmatrix(filename,'NumHeaderLines',6,'CommentStyle',{'/#'});
    f(:,6:end) = [];
    f(any(isnan(f),2),:) = []; 
    fout = reshape(f.',2045,101);
end

