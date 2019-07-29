function run_channel_searchlight_multiclass(varargin)

    %%
    if ismac
        if isempty(which('cosmo_wtf'))
            addpath('~/CoSMoMVPA/mvpa')
        end
        nproc = 2;
    else %on HPC
        addpath('../CoSMoMVPA/mvpa');
        addpath('../fieldtrip')
        % start cluster, give it a unique directory
        % starting a pool can fail when 2 procs are requesting simultaneous
        % thus try again after a second until success
        pool=[];
        while isempty(pool) 
            try
                pc = parcluster('local');
                pc.JobStorageLocation=tempdir;
                pool=parpool(pc);
            catch err
                disp(err)
                delete(gcp('nocreate'));
                pause(1)
            end
        end
        nproc=cosmo_parallel_get_nproc_available();
    end
    ft_defaults;

    opt = struct();
    opt = cosmo_structjoin(opt,varargin);
    subjectnr = opt.subject;
    
    %%
    fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',subjectnr);
    outfn = sprintf('results/sub-%02i_channel_searchlight_multiclass.mat',subjectnr);
    fprintf('loading %s\n',fn);tic
    load(fn,'ds')
    fprintf('loading data finished in %i seconds\n',ceil(toc))

%% slice to make this easier
% idx = ~mod(ds.sa.stimnumber,10);
% ds = cosmo_slice(ds,idx);
% idx = ~mod(ds.fa.time,20);
% ds = cosmo_dim_prune(cosmo_slice(ds,idx,2));

%% setup
targets = {ds.sa.stimnumber,...
    double(ismember(ds.sa.stimnumber,[0:29 60:89])),... %animacy
    double(ds.sa.stimnumber<60),... %big vs small
    };
targetlabels = {'image','animacy','size'};
stimconditionlabels = {'texform','object'};
durationconditionlabels = {'60Hz','30Hz','20Hz','5Hz'};

[c1,c2,c3]=meshgrid(1:length(targetlabels),1:length(stimconditionlabels),1:length(durationconditionlabels));
targetxcondi=[c1(:) c2(:) c3(:)];

%per stimcondition x durationcondition
res=cell(length(targetxcondi),1);
cctc = clock();
for tc = 1:length(targetxcondi)
    c1=targetxcondi(tc,1);
    c2=targetxcondi(tc,2);
    c3=targetxcondi(tc,3);
    mmtc = sprintf('%7s - %7s - %4s',targetlabels{c1},stimconditionlabels{c2},durationconditionlabels{c3});
    
    fprintf('\n\n%s decoding...\n',mmtc)
    ds.sa.targets = targets{c1};
    ds.sa.chunks = ds.sa.streamnumber;
    idxc2 = ds.sa.stimcondition==c2-1;
    idxc3 = ds.sa.durationcondition==c3-1;
    dsb = cosmo_slice(ds,idxc2 & idxc3);
    
    nh1 = cosmo_meeg_chan_neighborhood(dsb, 'count', 4);
    nh2 = cosmo_interval_neighborhood(dsb,'time','radius',0);
    nh = cosmo_cross_neighborhood(ds,{nh1,nh2});

    ma = struct();
    ma.classifier = @cosmo_classify_lda;
    ma.output = 'accuracy';
    ma.partitions = cosmo_nfold_partitioner(dsb);
    ma.nproc = nproc;

    r = cosmo_searchlight(dsb,nh,@cosmo_crossvalidation_measure,ma);
    
    sa = {};
    sa.c1 = repmat(c1,size(r.samples,1),1);
    sa.c2 = repmat(c2,size(r.samples,1),1);
    sa.c3 = repmat(c3,size(r.samples,1),1);
    sa.c1label = repmat(targetlabels(c1),size(r.samples,1),1);
    sa.c2label = repmat(stimconditionlabels(c2),size(r.samples,1),1);
    sa.c3label = repmat(durationconditionlabels(c3),size(r.samples,1),1);
    
    r.sa = cosmo_structjoin(r.sa,sa);    
    res{tc} = r;
    
    save(outfn,'res','-v7.3');
    
    fprintf('%s finished. Overall progress:\n',mmtc)
    cosmo_show_progress(cctc,tc/length(targetxcondi),sprintf('%i/%i',tc,length(targetxcondi)),'');
end
