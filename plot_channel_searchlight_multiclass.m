
%%
addpath('~/CoSMoMVPA/mvpa')
addpath('~/fieldtrip');ft_defaults

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
res = cosmo_stack(res_cell);

%%
layout=cosmo_meeg_find_layout(res);
timepoints = 0:50:450;
for c1=1:3
    for c2=1:2
        for c3=1:4
            %%
            x = cosmo_average_samples(cosmo_slice(res,res.sa.c1==c1 & res.sa.c2==c2 & res.sa.c3==c3),'split_by',{});
            % map to FT struct for visualization
            ft = cosmo_map2meeg(x);
                
            plotnr=0;
            f=figure(1);clf
            drawnow
            f.Position=[f.Position(1:2) 1000 500];
            f.Resize='off';
            f.PaperPositionMode='auto';f.PaperOrientation='portrait';
            for t=1:length(timepoints)
                plotnr=plotnr+1;
                subplot(2,5,plotnr)
                co = viridis();
                % show figure with plots for each sensor
                cfg = [];
                if c1==1
                    cfg.zlim = [1/120 .01];
                else
                    cfg.zlim = [1/2 .53];
                end
                cfg.layout = layout;
                cfg.style = 'straight';
                cfg.comment = 'no';
                cfg.xlim = timepoints(t) + [-10 10];
                cfg.colormap = co;
                ft_topoplotER(cfg, ft);
                title(sprintf('%i ms',timepoints(t)))
            end
            %% save
            fn = sprintf('figures/topomaps/figure_channel_searchlight_%s_%s_%s',x.sa.c1label{1},x.sa.c2label{1},x.sa.c3label{1});
            print(gcf,'-dpng','-r500',fn)
            im=imread([fn '.png']);
            [i,j]=find(mean(im,3)<255);margin=2;
            imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');
        end
    end
end
