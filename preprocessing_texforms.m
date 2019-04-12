function preprocessing_texforms(partid)

    %% eeglab
    if ~ismac
        addpath('../CoSMoMVPA/mvpa')
        addpath('../eeglab')
    end
    eeglab
    
    %% get files

    datapath = 'data';
    
    contfn = sprintf('%s/derivatives/eeglab/sub-%02i_task-rsvp_continuous.set',datapath,partid);
    if isfile(contfn)
        fprintf('Using %s\n',contfn)
    	EEG_cont = pop_loadset(contfn);
    else
        % load EEG file
        EEG_raw = pop_loadbv(sprintf('%s/sub-%02i/eeg/',datapath,partid), sprintf('sub-%02i_task-rsvp_eeg.vhdr',partid));
        EEG_raw = eeg_checkset(EEG_raw);
        EEG_raw.setname = partid;
        EEG_raw = eeg_checkset(EEG_raw);

        % high pass filter
        EEG_raw = pop_eegfiltnew(EEG_raw, 0.1,[]);

        % low pass filter
        EEG_raw = pop_eegfiltnew(EEG_raw, [],100);

        % downsample
        EEG_raw = pop_resample( EEG_raw, 250);
        EEG_raw = eeg_checkset(EEG_raw);

        % create eventlist
        EEG_cont = pop_creabasiceventlist( EEG_raw , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' });
        EEG_cont = eeg_checkset(EEG_cont);
        
        pop_saveset(EEG_cont,contfn);
    end
    
    %% add eventinfo to events
    eventsfncsv = sprintf('data/sub-%02i/eeg/sub-%02i_task-rsvp_events.csv',partid,partid);
    eventsfntsv = sprintf('data/sub-%02i/eeg/sub-%02i_task-rsvp_events.tsv',partid,partid);
    eventlist = readtable(eventsfncsv);
    %duration is a bids thing
    eventlist.stimduration = eventlist.duration;
    eventlist.duration=[];
    
    idx = strcmp({EEG_cont.event.codelabel},'E1');
    onset = vertcat(EEG_cont.event(idx).latency);
    duration = ones(size(onset));

    neweventlist = [table(onset,duration,'VariableNames',{'onset','duration'}) eventlist];
    
    tmp = [tempname '.csv'];
    writetable(neweventlist,tmp,'Delimiter','\t')
    movefile(tmp,eventsfntsv)
    
    %% run binlister for bins and extract bin-based epochs
    EEG_cont = pop_binlister( EEG_cont , 'BDF', sprintf('makebins.txt'), 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG', 'ExportEL', tempname);
    EEG_cont = pop_epochbin( EEG_cont , [-100  1000], 'none'); % no baseline correct
    EEG_cont = eeg_checkset(EEG_cont);
    
    %% convert to cosmo
    ds = cosmo_flatten(permute(EEG_cont.data,[3 1 2]),{'chan','time'},{{EEG_cont.chanlocs.labels},EEG_cont.times},2);
    ds.a.meeg=struct(); %or cosmo thinks it's not a meeg ds 
    ds.sa = table2struct(eventlist,'ToScalar',true);
    cosmo_check_dataset(ds,'meeg');
    
    %% save epochs
    save(sprintf('%s/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',datapath,partid),'ds','-v7.3')
    
end
