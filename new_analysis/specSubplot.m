function [] = specSubplot(input,tvec,fvec,row,col,ent)

    subplot(row,col,ent)
    imagesc(tvec,fvec,20*log10(abs(input./1E-6)))
    set(gca,'YDir','normal')
    colormap jet
    caxis([30 80])
    colorbar
    ylim([0 200])
    set(gca,'fontsize',20)
    

end
