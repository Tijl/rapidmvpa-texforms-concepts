function plot_decoding_all_combined()

    %% load stats
    x1=load('results/stats_decoding_pairwise.mat','stats');
    x2=load('results/stats_cross_decoding_pairwise.mat','stats');
    x3=load('results/stats_decoding_difference.mat','stats');
    
    %%
    res_cell={};cc=clock();mm='';
    for s=1:20
        fn = sprintf('results/sub-%02i_channel_searchlight_multiclass.mat',s);
        try
            load(fn,'res')
            res = cosmo_stack(res);
            res.sa.s = repmat(s,size(res.sa.c1));
            res_cell{end+1} = res;
        catch
        end
        mm=cosmo_show_progress(cc,s/20,[],mm);
    end
    res_searchlight = cosmo_stack(res_cell);
    layout=cosmo_meeg_find_layout(res_searchlight);

    %% one plot per contrast
    co=viridis(5);
    co2 = plasma(4);
    bfthresh = 10;
    tfcethresh = 1.6449;
	timevect = x1.stats{1}.timevect;
    for c1=1:3 %image / animacy / size
        
        %%
        f=figure(1);clf
        f.Position=[f.Position(1:2) 900 1100];
        f.Resize='off';f.PaperPositionMode='auto';
        plotnr=0;
        for plot2=1:6 %texform / object / difference / maps / cross texform / cross object / 
            if plot2~=4 %format: 1 big axes with 4 small axes for stats
                if plot2<3
                    c2 = plot2;
                    stats = x1.stats(c1,c2,:);
                elseif plot2>4
                    c2 = plot2-4;
                    stats = x2.stats(c1,c2,:);
                else
                    c2 = plot2-2;
                    stats = x3.stats(c1,:);
                end
                plotnr = plotnr+1;%ceil(c2/2)*4-mod(c2,2)-2;
                a = subplot(3,2,plotnr);hold on
                a.Position = a.Position+[0 .1 0 -.09];
                a.XLim = [min(timevect) max(timevect)];
                if plot2==3 % difference plot
                    a.YLim = [-.01 .06];
                    a.YTick=[-.01 :.01:.055];
                    %line at chance
                    plot(timevect,0*timevect,'k--');
                    ylabel('accuracy difference')
                else
                    a.YLim = [.49 .565];
                    a.YTick=[.5 :.02:.6];
                    %line at chance
                    plot(timevect,.5+0*timevect,'k--');
                    ylabel('decoding accuracy')
                end
                
                % error area
                for c3=1:4
                    s = stats{c3};
                    fill([timevect fliplr(timevect)],[s.mu-s.se fliplr(s.mu+s.se)],...
                        co(c3,:),'FaceAlpha',.2,'EdgeAlpha',0);
                end
                % means
                h=[];
                for c3=1:4
                    s = stats{c3};
                    c3t=s.c3label;
                    h(5-c3)=plot(timevect,s.mu,'Color',co(c3,:),'DisplayName',c3t,'LineWidth',1.5);
                end
                
                leg=legend(h);leg.Box='off';
                if plot2<3
                    c1t=s.c1label;
                    c2t=s.c2label; 
                    t=title(sprintf('%s    %ss: %s decoding',char('A'+plotnr-1),strrep(c2t,'object','intact object'),c1t));
                elseif plot2>4
                    c2t=s.c2label;
                    t=title(sprintf('%s    cross-decoding (test on %ss)',char('A'+plotnr-1),strrep(c2t,'object','intact object')));
                else
                    t=title(sprintf('%s    decoding difference (intact-texform)',char('A'+plotnr-1)));
                end
                t.HorizontalAlignment='left';
                t.Units = 'normalized';
                t.Position = [-.12,1.01,0];
                a.FontSize=8;
                t.FontSize=15;
                a.XTick=-100:100:1000;
                a.TickDir='out';
                %stats
                for c3=1:4
                    ab=.02;
                    ah=.025;
                    a2 = axes('Position',[a.Position(1:2)-[0 (5-c3)*ah+ab] a.Position(3) ah]);
                    hold on;
                    a2.XLim=a.XLim;
                    a2.XTick=a.XTick;
                    if c3==1
                        a2.XTickLabel = a.XTick;
                        xlabel('time (ms)')
                    else
                        a2.XTickLabel = [];
                    end
                    a2.YLim=[0.5,3.5];
                    a2.YTick=[1:3];
                    a2.YTickLabel = {'BF < 1/3','BF > 10','p < 0.05'};
                    a2.FontSize = 8;
                    a2.YAxis.FontSize=8;
                    a2.XAxis.FontSize=a.FontSize;
                    a2.TickDir='out';a2.TickLength=[.005 0];
                    for i=.5:1:3.5
                        plot(timevect,0*timevect+i,':','Color','k')
                    end
                    ms = 10;
                    s = stats{c3};
                    x = s.bf<1/3;
                    plot(timevect(x),0*timevect(x)+1,'.','Color',co(c3,:),'MarkerSize',ms)
                    x = s.bf>10;
                    plot(timevect(x),0*timevect(x)+2,'.','Color',co(c3,:),'MarkerSize',ms)
                    x = s.tfce_zval>1.6449;
                    plot(timevect(x),0*timevect(x)+3,'.','Color',co(c3,:),'MarkerSize',ms)
                end
            else %topomaps
                plotnr = plotnr+1;
                a = subplot(3,2,plotnr);axis off;
                a.Position = a.Position+[0 .1 0 -.09];
                t=title(sprintf('%s    channel decoding differences at 5Hz',char('A'+plotnr-1)));
                t.HorizontalAlignment='left';
                t.Units = 'normalized';
                t.Position = [-.12,1.01,0];
                t.FontSize=15;
                pos = a.Position;
                pos(2) = pos(2)+.06;
                aw = pos(3)/3;
                ah = pos(4)*.5;
                ahb = .02;
                
                timepoints = [0:50:400];
                xs1 = cosmo_average_samples(cosmo_slice(res_searchlight,res_searchlight.sa.c1==c1 & res_searchlight.sa.c2==1 & res_searchlight.sa.c3==4),'split_by',{});
                xs2 = cosmo_average_samples(cosmo_slice(res_searchlight,res_searchlight.sa.c1==c1 & res_searchlight.sa.c2==2 & res_searchlight.sa.c3==4),'split_by',{});
                ds = xs1;
                ds.samples = xs2.samples-xs1.samples;
                co2 = viridis();
                % map to FT struct for visualization
                ft = cosmo_map2meeg(ds);
                
                for t=1:9
                    a = axes('Position',[pos(1)+aw*mod(t-1,3) pos(2)-(ahb+ah)*floor((t-1)/3) aw ah]);
                    cfg = [];
                    if c1==1
                        cfg.zlim = [0 prctile(ds.samples,99)];
                    else
                        cfg.zlim = [0 prctile(ds.samples,99)];
                    end
                    cfg.layout = layout;
                    cfg.style = 'straight';
                    cfg.comment = 'no';
                    cfg.xlim = timepoints(t) + [-10 10];
                    cfg.colormap = co2;
                    ft_topoplotER(cfg, ft);
                    text(mean(a.XLim),min(a.YLim),sprintf('%i ms',timepoints(t)),'VerticalAlignment','top','HorizontalAlignment','center')                    
                    if t==7
                        c = colorbar('Position',[a.Position(1)-.025 a.Position(2)+.5*ahb .01 3*ah]);
                        c.Label.String={'accuracy difference','intact - texform'};
                    end
                end
            end
            drawnow;
        end
        %% save
        fn = sprintf('figures/figure_decoding_combined_%s',c1t);
        print(gcf,'-dpng','-r500',fn)
        im=imread([fn '.png']);
        [i,j]=find(mean(im,3)<255);margin=2;
        imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');
    end