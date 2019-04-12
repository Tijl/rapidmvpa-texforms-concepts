function plot_decoding_pairwise()

    %% load stats
    load('results/stats_decoding_pairwise.mat','stats')
    
    %% plot per category
    
    co=viridis(5);
    co2 = plasma(4);
    bfthresh = 3;
	timevect = stats{1}.timevect;
    for c1=1:3
        f=figure(1);clf
        f.Position=[f.Position(1:2) 900 280];
        f.Resize='off';f.PaperPositionMode='auto';
        plotnr=0;
        for c2=1:2
            plotnr = plotnr+1;
            a=subplot(1,2,plotnr); hold on
            
            a.XLim = minmax(timevect);
            a.YLim = [.45 .565];
            st = .485; %upper limit
            lh = diff(a.YLim)/60;%.003; %spacing between dot lines
            li = -diff(a.YLim)/13; %spacing between triplets
            ms = 2; %dot size
            zc = .9*[1 1 1]; %empty space bf colour         
            
            %line at chance
            plot(timevect,.5+0*timevect,'k--');
            
            % error area
            for c3=1:4
                s = stats{c1,c2,c3};
                fill([timevect fliplr(timevect)],[s.mu-s.se fliplr(s.mu+s.se)],...
                    co(c3,:),'FaceAlpha',.2,'EdgeAlpha',0);
            end
            
            % means
            h=[];
            for c3=1:4
                s = stats{c1,c2,c3};
                c3t=s.c3label;
                h(5-c3)=plot(timevect,s.mu,'Color',co(c3,:),'DisplayName',c3t,'LineWidth',1.5);
            end
            
            if ismember(plotnr,[1 3 5])
                %connecting lines
                for z=-2:1
                    i=3;
                    plot([-190 min(timevect)],[.47+z*.007 (st+lh*z+li*(i-1))],'-o','Color',zc.*.9,'LineWidth',ms,'Clipping','off','MarkerSize',ms,'MarkerFaceColor',zc*.9);
                end
            end
            
            % BF
            for c3=1:4
                s = stats{c1,c2,c3};
                x = zeros(size(s.bf));
                x(s.bf>bfthresh)=1;
                x(s.bf<1)=-1;
                x(s.bf<(1/bfthresh))=-2;
                c = 5-c3;
                for z=-2:1
                    plot(minmax(timevect),(st+lh*z+li*(c-1))*[1 1],'-','Color',zc,'LineWidth',ms);
                end
                plot(timevect,st+lh*x+li*(c-1),'o','Color',co(c3,:),'MarkerSize',ms);
                plot(timevect(x>0 | x<-1), st+lh*x(x>0 | x<-1)+li*(c-1),'o','Color',co(c3,:),'MarkerSize',ms,'MarkerFaceColor',co(c3,:));
                tt=text(max(timevect),st+li*(c-1),['  BF ' s.c3label],'VerticalAlignment','middle','FontSize',8,'Color',co(c3,:));
            end
            
            leg=legend(h);leg.Box='off';
            c1t=s.c1label;
            c2t=s.c2label;
            t=title(sprintf('%s    %ss: %s decoding',char('A'+plotnr-1),strrep(c2t,'object','intact object'),c1t));
            t.HorizontalAlignment='left';
            t.Units = 'normalized';
            t.Position = [-.12,1.01,0];
            a.FontSize=11;
            t.FontSize=15;
            xlabel('time (ms)')
            ylabel('decoding accuracy')
            a.YTick=[.5 :.02:.6];
            a.XTick=-100:100:1000;
            a.YLabel.Position=[a.YLabel.Position(1) .53 -1];
            
            if ismember(plotnr,[1 3 5])
                % bf box
                a2 = axes('Position',[a.Position(1)-.085 a.Position(2) .06 .27],'box','on');
                %leg.Position = a2.Position+[0 a2.Position(4) 0 -.12];
                a2.XTick=[];a2.YTick=[];
                a2.YLim=[.05 1.1];a2.XLim=[0 1];hold on
                locs = linspace(0.1,.9,5);
                mlocs = movmean(locs,2);mlocs = mlocs(2:end);
                tt = fliplr({sprintf('BF>%i',bfthresh),'BF>1','BF<1',sprintf('BF<1/%i',bfthresh)});
                bfc = co(2,:);
                for y=1:length(locs)
                    if y==3
                        plot([0.05 .95],locs(y)*[1 1],'-.','Color',bfc,'LineWidth',.8)
                    else
                        plot([0.05 .95],locs(y)*[1 1],'Color',bfc,'LineWidth',.8)
                    end
                    if y<5
                        if ismember(y,[2 3])
                            plot(.2,mlocs(y),'o','Color',bfc,'MarkerSize',6)
                        else
                            plot(.2,mlocs(y),'o','Color',bfc,'MarkerSize',6,'MarkerFaceColor',bfc)
                        end
                        t=text(1/3,mlocs(y),tt{y},'Color',bfc,'FontSize',8);

                    end
                end
                t=text(.05,1,'Bayes Factor:','FontSize',8,'Color',bfc);
            end
            
        end
        
        % save
        fn = sprintf('figures/figure_decoding_pairwise_%s',c1t);
        print(gcf,'-dpng','-r500',fn)
        im=imread([fn '.png']);
        [i,j]=find(mean(im,3)<255);margin=2;
        imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');
    
    end