%% load images

files = dir('Stimuli/*/*/*.png');
IM = {};
for f=1:length(files)
    IM{f} = imread(fullfile(files(f).folder,files(f).name));
    %IM{f} = IM{f};
end

%%
f=figure(1);clf
f.Position=[f.Position(1:2) 1000 1200];
f.Resize='off';
f.PaperPositionMode='auto';
f.PaperOrientation='portrait';

a1=subplot(2,2,1);
montage(IM(121:240), 'Size', [12,10], 'BorderSize', 1)
title('A    Texform objects','Units','Normalized','Position',[0 1.01 0],'FontSize',16,'HorizontalAlignment','left');

a2=subplot(2,2,2);
montage(IM(1:120), 'Size', [12,10], 'BorderSize', 1)
title('B    Intact objects','Units','Normalized','Position',[0 1.01 0],'FontSize',16,'HorizontalAlignment','left');

left = a1.Position(1)+a1.Position(3);
width = a2.Position(1)-left;
bottom = a1.Position(2)+.0035;
height = (a1.Position(4)-.007)*.25;

annotation('textbox',[left,bottom+3*height,width,height],'String',{'Big','Animals'},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',16);
annotation('textbox',[left,bottom+2*height,width,height],'String',{'Big','Objects'},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',16);
annotation('textbox',[left,bottom+height,width,height],'String',{'Small','Animals'},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',16);
annotation('textbox',[left,bottom,width,height],'String',{'Small','Objects'},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',16);

rng(187356);
idx = randsample(120,10);

subplot(7,1,4)
montage(IM(120+idx), 'Size', [1,10], 'BorderSize', 1);hold on
title('C    Example 5Hz sequence (texforms)','Units','Normalized','Position',[0 1.01 0],'FontSize',16,'HorizontalAlignment','left');

a=subplot(7,1,5);
m=montage(IM(idx), 'Size', [1,10], 'BorderSize', 1);hold on
title('D    Example 5Hz sequence (intact objects)','Units','Normalized','Position',[0 1.01 0],'FontSize',16,'HorizontalAlignment','left');

for i = 1:2
    a = subplot(7,1,i+3);
    x = movmean(linspace(a.XLim(1),a.XLim(2),11),2);
    text(x(2:end),repmat(a.YLim(2),1,10),'200ms','HorizontalAlignment','center','VerticalAlignment','top');
    plot(x(2:end),repmat(mean(a.YLim),1,10),'.','MarkerSize',10,'MarkerFaceColor',.99*[1 1 1],'MarkerEdgeColor',.99*[1 1 1])
    plot(x(8),mean(a.YLim),'.','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r');drawnow
    text(x(8),min(a.YLim),'button press','HorizontalAlignment','center','VerticalAlignment','bottom')
end

%% save
fn = sprintf('figures/figure_design');
fprintf('Saving.');print(gcf,'-dpng','-r500',fn);fprintf('.')
im=imread([fn '.png']);fprintf('.')
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');fprintf('.Done\n')