% Ser range of imaging
zRangeMin = 0.012;
zRangeMax = 0.035;
zIndMin = round(zRangeMin/dz)+1;
zIndMax = round(zRangeMax/dz)+1;
zIndRangeShow = zIndMin:zIndMax;

filtSizeL = [round(data.lambda2L/dz/2), round(data.lambda2L/dx/2)];
filtSizeT = [round(data.lambda2T/dz/2), round(data.lambda2T/dx/2)];

% Show images
crange = max(abs(BFdata_LdL(:)));
if (ID == 'DV15143')
    crange = crange/2;
end

% LdL mode - zero angle
BFdata_LdL_f = wiener2(abs(data.LdL), filtSizeL);
env = pipeMask.*abs(BFdata_LdL_f);
figure;
imagesc(x_axis(xIndRange)*1e3,z_axis(zIndRangeShow)*1e3,env(zIndRangeShow,:));
caxis([0 crange]/2); colormap(gca,"hot"); axis image;
if (flag)
    hold on;
    plot(x_axis(xIndRange)*1e3,frontDist(xIndRange)*1e3,'w:','LineWidth',1);
    plot(x_axis(xIndRange)*1e3,backDist(xIndRange)*1e3,'w:','LineWidth',1);
    hold off;
end
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
title(sprintf('L-L Mode: %.1f deg', rad2deg(TXangles(2))),'FontName','Courier','FontSize',14,'FontWeight','bold');
set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-32 -24 -16 -8 0 8 16 24 32]);
saveas(gcf,sprintf('images/%s/das256_[%s]_frame%d_LdL_Env.fig',ID,ID,frameNum));
exportgraphics(gcf, sprintf('images/%s/das256_[%s]_frame%d_LdL_Env.png',ID,ID,frameNum), 'ContentType', 'image', 'BackgroundColor', 'white');

% LLdLL mode - zero angle
figure;
BFdata_LLdLL_f = wiener2(BFdata_LLdLL, filtSizeL);
env = pipeMask.*abs(BFdata_LLdLL_f);%
imagesc(x_axis(xIndRange)*1e3,z_axis(zIndRangeShow)*1e3,env(zIndRangeShow,:));
caxis([0 crange]/2); colormap(gca,"hot"); axis image;
if (flag)
    hold on;
    plot(x_axis(xIndRange)*1e3,frontDist(xIndRange)*1e3,'w:','LineWidth',1);
    plot(x_axis(xIndRange)*1e3,backDist(xIndRange)*1e3,'w:','LineWidth',1);
    hold off;
end
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
title(sprintf('LL-LL Mode: %.1f deg', rad2deg(TXangles(2))),'FontName','Courier','FontSize',14,'FontWeight','bold');%
set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-32 -24 -16 -8 0 8 16 24 32]);
saveas(gcf,sprintf('images/%s/das256_[%s]_frame%d_LLdLL_Env.fig',ID,ID,frameNum));
exportgraphics(gcf, sprintf('images/%s/das256_[%s]_frame%d_LLdLL_Env.png',ID,ID,frameNum), 'ContentType', 'image', 'BackgroundColor', 'white');

% TdT mode - negative angle
BFdata_TdT_neg_f = wiener2(abs(gradient(data.TdT_neg)), filtSizeT);
env = pipeMask.*BFdata_TdT_neg_f;
figure;
imagesc(x_axis(xIndRange)*1e3,z_axis(zIndRangeShow)*1e3,env(zIndRangeShow,:));
caxis([0 crange]/2); colormap(gca,"hot"); axis image;
if (flag)
    hold on;
    plot(x_axis(xIndRange)*1e3,frontDist(xIndRange)*1e3,'w:','LineWidth',1);
    plot(x_axis(xIndRange)*1e3,backDist(xIndRange)*1e3,'w:','LineWidth',1);
    hold off;
end
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
title(sprintf('T-T Mode: %.1f deg', rad2deg(TXangles(1))),'FontName','Courier','FontSize',14,'FontWeight','bold');%
set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-32 -24 -16 -8 0 8 16 24 32]);
saveas(gcf,sprintf('images/%s/das256_[%s]_frame%d_TdT_neg_Env.fig',ID,ID,frameNum));
exportgraphics(gcf, sprintf('images/%s/das256_[%s]_frame%d_TdT_neg_Env.png',ID,ID,frameNum), 'ContentType', 'image', 'BackgroundColor', 'white');

% TdT mode - positive angle
BFdata_TdT_pos_f = wiener2(abs(gradient(data.TdT_pos)), filtSizeT);
env = pipeMask.*BFdata_TdT_pos_f;
figure;
imagesc(x_axis(xIndRange)*1e3,z_axis(zIndRangeShow)*1e3,env(zIndRangeShow,:));
caxis([0 crange]/2); colormap(gca,"hot"); axis image;
if (flag)
    hold on;
    plot(x_axis(xIndRange)*1e3,frontDist(xIndRange)*1e3,'w:','LineWidth',1);
    plot(x_axis(xIndRange)*1e3,backDist(xIndRange)*1e3,'w:','LineWidth',1);
    hold off;
end
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
title(sprintf('T-T Mode: %.1f deg', rad2deg(TXangles(3))),'FontName','Courier','FontSize',14,'FontWeight','bold');%
set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-32 -24 -16 -8 0 8 16 24 32]);
saveas(gcf,sprintf('images/%s/das256_[%s]_frame%d_TdT_pos_Env.fig',ID,ID,frameNum));
exportgraphics(gcf, sprintf('images/%s/das256_[%s]_frame%d_TdT_pos_Env.png',ID,ID,frameNum), 'ContentType', 'image', 'BackgroundColor', 'white');

% TTdT mode - coherent compounding 
BFdata_TTdT_f = wiener2(abs(gradient(data.TTdT_neg+data.TTdT_pos)), filtSizeT);
env = pipeMask.*BFdata_TTdT_f/2;
figure;
imagesc(x_axis(xIndRange)*1e3,z_axis(zIndRangeShow)*1e3,env(zIndRangeShow,:));
caxis([0 crange]/2); colormap(gca,"hot"); axis image;
if (flag)
    hold on;
    plot(x_axis(xIndRange)*1e3,frontDist(xIndRange)*1e3,'w:','LineWidth',1);
    plot(x_axis(xIndRange)*1e3,backDist(xIndRange)*1e3,'w:','LineWidth',1);
    hold off;
end
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
title(sprintf('TT-T Mode (compounded)'),'FontName','Courier','FontSize',14,'FontWeight','bold');%
set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-32 -24 -16 -8 0 8 16 24 32]);
saveas(gcf,sprintf('images/%s/das256_[%s]_frame%d_TTdT_comp_Env.fig',ID,ID,frameNum));
exportgraphics(gcf, sprintf('images/%s/das256_[%s]_frame%d_TTdT_comp_Env.png',ID,ID,frameNum), 'ContentType', 'image', 'BackgroundColor', 'white');

% TdTT mode - compounded 
BFdata_TdTT_f = wiener2(abs(gradient(data.TdTT_neg+data.TdTT_pos)), filtSizeT);
env = pipeMask.*BFdata_TdTT_f/2;
figure;
imagesc(x_axis(xIndRange)*1e3,z_axis(zIndRangeShow)*1e3,env(zIndRangeShow,:));
caxis([0 crange]/2); colormap(gca,"hot"); axis image;
if (flag)
    hold on;
    plot(x_axis(xIndRange)*1e3,frontDist(xIndRange)*1e3,'w:','LineWidth',1);
    plot(x_axis(xIndRange)*1e3,backDist(xIndRange)*1e3,'w:','LineWidth',1);
    hold off;
end
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
title(sprintf('T-TT Mode (compounded)'),'FontName','Courier','FontSize',14,'FontWeight','bold');%
set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-32 -24 -16 -8 0 8 16 24 32]);
saveas(gcf,sprintf('images/%s/das256_[%s]_frame%d_TdTT_comp_Env.fig',ID,ID,frameNum));
exportgraphics(gcf, sprintf('images/%s/das256_[%s]_frame%d_TdTT_comp_Env.png',ID,ID,frameNum), 'ContentType', 'image', 'BackgroundColor', 'white');

% TTdTT mode negative angle
BFdata_TTdTT_neg_f = wiener2(abs(gradient(data.TTdTT_neg)), filtSizeT);
env = pipeMask.*BFdata_TTdTT_neg_f;
figure;
imagesc(x_axis(xIndRange)*1e3,z_axis(zIndRangeShow)*1e3,env(zIndRangeShow,:));
caxis([0 crange]/2); colormap(gca,"hot"); axis image;
if (flag)
    hold on;
    plot(x_axis(xIndRange)*1e3,frontDist(xIndRange)*1e3,'w:','LineWidth',1);
    plot(x_axis(xIndRange)*1e3,backDist(xIndRange)*1e3,'w:','LineWidth',1);
    hold off;
end
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
title(sprintf('TT-TT Mode: %.1f deg',rad2deg(TXangles(1))),'FontName','Courier','FontSize',14,'FontWeight','bold');%
set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-32 -24 -16 -8 0 8 16 24 32]);
saveas(gcf,sprintf('images/%s/das256_[%s]_frame%d_TTdTT_neg_Env.fig',ID,ID,frameNum));
exportgraphics(gcf, sprintf('images/%s/das256_[%s]_frame%d_TTdTT_neg_Env.png',ID,ID,frameNum), 'ContentType', 'image', 'BackgroundColor', 'white');

% TTdTT mode positive angle
BFdata_TTdTT_pos_f = wiener2(abs(gradient(data.TTdTT_pos)), filtSizeT);
env = pipeMask.*BFdata_TTdTT_pos_f;
figure;
imagesc(x_axis(xIndRange)*1e3,z_axis(zIndRangeShow)*1e3,env(zIndRangeShow,:));
caxis([0 crange]/2); colormap(gca,"hot"); axis image;
if (flag)
    hold on;
    plot(x_axis(xIndRange)*1e3,frontDist(xIndRange)*1e3,'w:','LineWidth',1);
    plot(x_axis(xIndRange)*1e3,backDist(xIndRange)*1e3,'w:','LineWidth',1);
    hold off;
end
set(gca,'FontName','Courier','FontSize',12,'FontWeight','demi');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('x (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
ylabel('z (mm)','FontName','Courier','FontSize',14,'FontWeight','bold');
title(sprintf('TT-TT Mode: %.1f deg',rad2deg(TXangles(3))),'FontName','Courier','FontSize',14,'FontWeight','bold');%
set(gca,'YTick',[5 15 25 35 45]); set(gca,'XTick',[-32 -24 -16 -8 0 8 16 24 32]);
saveas(gcf,sprintf('images/%s/das256_[%s]_frame%d_TTdTT_pos_Env.fig',ID,ID,frameNum));
exportgraphics(gcf, sprintf('images/%s/das256_[%s]_frame%d_TTdTT_pos_Env.png',ID,ID,frameNum), 'ContentType', 'image', 'BackgroundColor', 'white');







