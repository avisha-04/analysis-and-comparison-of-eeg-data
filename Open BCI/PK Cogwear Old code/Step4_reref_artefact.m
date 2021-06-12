% Clear Memory and MatLab Processes
close all; clear; clc;

% Loop through files
folder = 'C:\\Users\\yangg\\Documents\\PK Cogwear\\filtered\\';
files = dir(fullfile(folder, '*.csv'));

% initialize artefact tracker matrix
artefact_tracker =zeros(length(files), 3);

for k = 1:length(files)
    baseFileName = files(k).name;
    FileName = fullfile(folder, baseFileName);
    
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    sfreq = 200
    EEG = pop_importdata('dataformat','ascii','nbchan',0,'data',FileName,'srate',sfreq,'pnts',0,'xmin',0,'chanlocs','C:\\Users\\yangg\\Documents\\PK Cogwear\\cogwear_chanloc_04.ced');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','filt','gui','off'); 
    EEG = eeg_checkset( EEG );
    EEG = pop_reref( EEG, 2);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','reref','gui','off'); 

    % Normalizing w.r.t. the Baselines

    %eeg_baseline = EEG.data(:, 2375:2500);

    %EEG.data = EEG.data - mean(eeg_baseline, 2) ./ std(eeg_baseline, [], 2);   

    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    pop_processMARA ( ALLEEG,EEG,CURRENTSET );

    % Artefact Rejection
    prob_threshold = 0.7;

    % Run ICA
    EEG = pop_runica(EEG, 'extended', 1, 'icatype', 'runica');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );

     % Process MARA
    [~, MARA_info] = MARA(EEG);
    [min_probability, temp] = min(MARA_info.posterior_artefactprob);
    MARA_info.posterior_artefactprob(temp(1)) = MARA_info.posterior_artefactprob(temp) - 0.0001;

    artefacts = MARA_info.posterior_artefactprob >= prob_threshold...
        & MARA_info.posterior_artefactprob >= min_probability;
    
    artefact_tracker(k, :) =  artefacts;
    
    % Reject IC's
    EEG = pop_subcomp( EEG, find(artefacts, length(artefacts)-1) );

    % Output
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'data_pre_pruned', 'gui', 'off');

    % Filter
    EEG = pop_eegfiltnew(EEG, 1,45,826,0,[],0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','rect_filt','gui','off');
    
    % Close EEGLAB and all windows
    close all;   
    
    % Save output file
    output = EEG.data;
    outfile = strcat('C:\\Users\\yangg\\Documents\\PK Cogwear\\processed\\',baseFileName);
    writematrix(output, outfile);
    
end

% Save artefact tracker
writematrix(artefact_tracker, 'C:\\Users\\yangg\\Documents\\PK Cogwear\\artefacts.csv');
save('filesartefact.mat','files')
