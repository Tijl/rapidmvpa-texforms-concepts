function stats_cross_decoding_pairwise()

    if isempty(which('cosmo_wtf'))
        addpath('~/CoSMoMVPA/mvpa')
    end
    
    r=.707;
    
    %% load data
    fprintf('Loading data\n')
    nsubjects=20;
    res_all={};cc = clock();mm='';
    for f=1:nsubjects
        fn = sprintf('results/sub-%02i_cross_decoding_pairwise.mat',f);
        load(fn,'res')
        res_avg={};
        for i = 1:length(res)
            if ~isempty(res{i})
                res_avg{i} = cosmo_average_samples(res{i},'split_by',{});
            end
        end
        res_avg = cosmo_stack(res_avg);
        res_avg.sa.snum = f*ones(size(res_avg.samples,1),1);
        res_all{f} = res_avg;
        mm=cosmo_show_progress(cc,f/nsubjects,sprintf('%i/%i',f,nsubjects),mm);
    end
    
    res_all=cosmo_stack(res_all);
    
    %% compute stats
    
    fprintf('Computing stats\n')
    stats = {};
    timevect = res_all.a.fdim.values{1};
    cc = clock();mm='';
    for c1=1:3
        for c2=1:2
            for c3=1:4
                idx = res_all.sa.c1==c1 & res_all.sa.c2==c2 & res_all.sa.c3==c3;
                res_idx = cosmo_slice(res_all,idx);
                x = res_idx.samples;

                opt={};
                opt.niter = 1000;
                opt.h0_mean = .5;
                res_idx.sa.targets = ones(size(res_idx.sa.snum));
                res_idx.sa.chunks = cumsum(ones(size(res_idx.sa.snum)));
                ds_tfce = cosmo_montecarlo_cluster_stat(res_idx,cosmo_cluster_neighborhood(res_idx),opt);
                                
                s = struct();
                s.n = size(x,1);
                s.mu = mean(x);
                s.se = std(x)./sqrt(s.n);
                h0mean = .5;
                s.tstat = (s.mu-h0mean)./s.se;
                s.bf = t1smpbf(s.tstat,s.n,r);
                s.tfce_zval = ds_tfce.samples;
                s.c1 = c1;
                s.c2 = c2;
                s.c3 = c3;
                s.c1label = res_idx.sa.c1label{1};
                s.c2label = res_idx.sa.c2label{1};
                s.c3label = res_idx.sa.c3label{1};
                s.timevect = timevect;
                stats{c1,c2,c3} = s;
            end
        end
        mm = cosmo_show_progress(cc,c1/3,'',mm);
    end
    
    fprintf('Saving\n')
    save('results/stats_cross_decoding_pairwise.mat','stats')
    fprintf('Done\n')
    