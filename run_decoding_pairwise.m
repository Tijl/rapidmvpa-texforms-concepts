function run_decoding_pairwise(varargin)

    %%
    if ismac
        if isempty(which('cosmo_wtf'))
            addpath('~/CoSMoMVPA/mvpa')
        end
        nproc = 2;
    else %on HPC
        addpath('../CoSMoMVPA/mvpa');
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

    opt = struct();
    opt = cosmo_structjoin(opt,varargin);
    subjectnr = opt.subject;
    
    %%
    fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',subjectnr);
    outfn = sprintf('results/sub-%02i_decoding_pairwise.mat',subjectnr);
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
    
    nh = cosmo_interval_neighborhood(dsb,'time','radius',0);

    ma = struct();
    ma.classifier = @cosmo_classify_lda;
    
    %create partitioning scheme
    ut = unique(dsb.sa.targets);
    % all chunks to leave out
    uc = unique(dsb.sa.chunks);
    sa=struct('target1',[],'target2',[],'leftoutchunk',[],'leftoutexemplar1',[],'leftoutexemplar2',[]);
    ma.output = 'fold_accuracy';
    ma.partitions = struct();
    ma.partitions.train_indices = {};
    ma.partitions.test_indices = {};
    if length(ut)==2
        % for categorical contrasts: exemplar-by-sequence
        dsb.sa.cvtargets = dsb.sa.stimnumber;
        % ue1 and ue2 are the unique exemplars (to leave out in the test set)
        ue1 = unique(dsb.sa.cvtargets(dsb.sa.targets==ut(1)));
        ue2 = unique(dsb.sa.cvtargets(dsb.sa.targets==ut(2)));
        % leave combinations of exemplar pairs out once
        ue = [ue1 ue2];
        for j=1:length(uc) % for each chunk to leave out
            idx_chunk = dsb.sa.chunks==uc(j); % find chunk to leave out
            for k=1:size(ue,1) % for each exemplar pair to leave out
                % store left out chunk and exemplar in result
                sa.target1(end+1,1) = ut(1);
                sa.target2(end+1,1) = ut(2);
                sa.leftoutchunk(end+1,1) = uc(j);
                sa.leftoutexemplar1(end+1,1) = ue(k,1);
                sa.leftoutexemplar2(end+1,1) = ue(k,2);
                % set partitions
                idx_ex = ismember(dsb.sa.cvtargets,ue(k,:));
                ma.partitions.train_indices{1,end+1} = find(~idx_chunk & ~idx_ex);
                ma.partitions.test_indices{1,end+1} = find(idx_chunk & idx_ex);
            end
        end
    else
        % all pairwise combinations
        combs = combnk(unique(dsb.sa.targets,'rows'),2);
        ma.check_partitions = false;
        for i=1:length(combs) % for each pair
            idx_ex = ismember(dsb.sa.targets,combs(i,1)) | ismember(dsb.sa.targets,combs(i,2));
            for j=1:length(uc) % for each chunk to leave out
                idx_chunk = dsb.sa.chunks==uc(j); % find chunk to leave out
                % store left out chunk and exemplar in result
                sa.target1(end+1,1) = combs(i,1);
                sa.target2(end+1,1) = combs(i,2);
                sa.leftoutchunk(end+1,1) = uc(j);
                sa.leftoutexemplar1(end+1,1) = combs(i,1);
                sa.leftoutexemplar2(end+1,1) = combs(i,2);
                % set partitions
                ma.partitions.train_indices{1,end+1} = find(~idx_chunk & idx_ex);
                ma.partitions.test_indices{1,end+1} = find(idx_chunk & idx_ex);
            end
        end
    end
    ma.nproc = nproc;

    r = cosmo_searchlight(dsb,nh,@cosmo_crossvalidation_measure,ma);
    
    sa.c1 = repmat(c1,size(r.samples,1),1);
    sa.c2 = repmat(c2,size(r.samples,1),1);
    sa.c3 = repmat(c3,size(r.samples,1),1);
    sa.c1label = repmat(targetlabels(c1),size(r.samples,1),1);
    sa.c2label = repmat(stimconditionlabels(c2),size(r.samples,1),1);
    sa.c3label = repmat(durationconditionlabels(c3),size(r.samples,1),1);
    
    r.sa = cosmo_structjoin(r.sa,sa);
    r2 = cosmo_average_samples(r,'split_by',{'target1','target2'});
    
    res{tc} = r2;
    
    save(outfn,'res','-v7.3');
    
    fprintf('%s finished. Overall progress:\n',mmtc)
    cosmo_show_progress(cctc,tc/length(targetxcondi),sprintf('%i/%i',tc,length(targetxcondi)),'');
end
