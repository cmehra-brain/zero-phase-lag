%% functional networks PLV + wPLI methods paper

% This script takes pre-processed EEG data and produces functional
% connectivity matrices

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript: 
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------


% all these methods use the hilbert transform
% PLV + wPLI in 5 frequency bins
% 2 seconds of data padding either sides of 20 seconds resting state eyes
% closed data. This padding is removed once filtering and hilbert transform
% are done.

%% clean and set up drives
clear
clc

%set up folders
path2data = '/home/.../EEG_files/EEG_cropped_T1_DK_1200';
addpath([path2data,'/'])
addpath('/home/.../EEG_files/Script/FunNets/')

path2output_PLV_1_4 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/PLV/1_4Hz';
path2output_PLV_4_8 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/PLV/4_8Hz';
path2output_PLV_8_13 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/PLV/8_13Hz';
path2output_PLV_13_20 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/PLV/13_20Hz';
path2output_PLV_20_32 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/PLV/20_32Hz';

path2output_wPLI_1_4 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/wPLI/1_4Hz';
path2output_wPLI_4_8 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/wPLI/4_8Hz';
path2output_wPLI_8_13 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/wPLI/8_13Hz';
path2output_wPLI_13_20 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/wPLI/13_20Hz';
path2output_wPLI_20_32 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/wPLI/20_32Hz';

ROI_files = fullfile(path2data, '/', 'ROIdata*');
EEG_list = dir(ROI_files); %directory lists files and folders in the current folder

%% loop through each cropped EEG file
for i = 1:length(EEG_list)
    %%
    disp(i)
    cd(path2data);
    load(EEG_list(i).name);
    
    % --- 1. Ensure that each row is a EEG recording ---
    % if the rows are longer than the columns --> flip them around
    if size(ROIdata,1)>size(ROIdata,2)
        ROIdata=ROIdata';
    end
    
    %%
    % --- 2. Cut data in 20 segments --------
    T_min=4;    % minimum time required for a segment
    T_seg = 24; % length in seconds for each segment
    [ROIdata,n_seg]=cut_segments(ROIdata,sf,T_seg,T_min);
    %n_seg tells you how many segments of 20 seconds
    %this cuts your segments into the full recording into 20 seconds.
    % returns ROI data as a 3D matrix. each segment is in its own dimension.
    
    %visualise raw EEG data after - realies on eeg_visual.m
    %script from Petroula
    %eeg_visual(ROIdata(:,:,1),sf)
    
      
    %% FILTERING AND CALCULATING PLV + wPLI IN ONE STEP 1-4Hz
    
    networks_1_4_PLV = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_1_4_wPLI = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    
    %networks_1_4_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    
    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,4,sf,250,'low'); %creates filtered epochs of the right length
        PLV=plv_undirected_with_zero_lag_no_surrogates(mysegmentsignal); % calculate PLV for each 20 second segment.
        %no surrogates used | padding in data removed after doing the hilbert transform
        networks_1_4_PLV(:,:,count_segments) = PLV;
        
        wPLI=weighted_phase_lag_index_padding(mysegmentsignal); % calculate PLV for each 20 second segment.
        networks_1_4_wPLI(:,:,count_segments) = wPLI; %calculate the PLV for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform
        
        
    end
    
    %save PLV
    cd(path2output_PLV_1_4)
    save(['nets_',EEG_list(i).name],'networks_1_4_PLV');

    %save wPLI
    cd(path2output_wPLI_1_4)
    save(['nets_',EEG_list(i).name],'networks_1_4_wPLI');
    
    %% 4-8 Hz
    
    networks_4_8_PLV = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_4_8_wPLI = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    
    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,4,sf,250,'high');
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,8,sf,250,'low'); %creates filtered epochs of the right length
        PLV=plv_undirected_with_zero_lag_no_surrogates(mysegmentsignal); % calculate PLV for each 20 second segment.
        %no surrogates used | padding in data removed after doing the hilbert transform
        networks_4_8_PLV(:,:,count_segments) = PLV;
        
        wPLI=weighted_phase_lag_index_padding(mysegmentsignal); % calculate PLV for each 20 second segment.
        networks_4_8_wPLI(:,:,count_segments) = wPLI; %calculate the PLV for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform
    end
    
    %save PLV
    cd(path2output_PLV_4_8)
    save(['nets_',EEG_list(i).name],'networks_4_8_PLV');
    
    %save wPLI
    cd(path2output_wPLI_4_8)
    save(['nets_',EEG_list(i).name],'networks_4_8_wPLI');
    
    %% 8-13 Hz
    
    networks_8_13_PLV = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_8_13_wPLI = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    
    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,8,sf,250,'high');
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,13,sf,250,'low'); %creates filtered epochs of the right length
        PLV=plv_undirected_with_zero_lag_no_surrogates(mysegmentsignal); % calculate PLV for each 20 second segment.
        %no surrogates used | padding in data removed after doing the hilbert transform
        networks_8_13_PLV(:,:,count_segments) = PLV;
        
        wPLI=weighted_phase_lag_index_padding(mysegmentsignal); % calculate PLV for each 20 second segment.
        networks_8_13_wPLI(:,:,count_segments) = wPLI; %calculate the PLV for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform
    end
    
    %save PLV
    cd(path2output_PLV_8_13)
    save(['nets_',EEG_list(i).name],'networks_8_13_PLV');
    
    %save wPLI
    cd(path2output_wPLI_8_13)
    save(['nets_',EEG_list(i).name],'networks_8_13_wPLI');
    
    %% 13-20 Hz
    
    networks_13_20_PLV = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_13_20_wPLI = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    
    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,13,sf,250,'high');
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,20,sf,250,'low'); %creates filtered epochs of the right length
        PLV=plv_undirected_with_zero_lag_no_surrogates(mysegmentsignal); % calculate PLV for each 20 second segment.
        %no surrogates used | padding in data removed after doing the hilbert transform
        networks_13_20_PLV(:,:,count_segments) = PLV;
        
        wPLI=weighted_phase_lag_index_padding(mysegmentsignal); % calculate PLV for each 20 second segment.
        networks_13_20_wPLI(:,:,count_segments) = wPLI; %calculate the PLV for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform
    end
    
    %save PLV
    cd(path2output_PLV_13_20)
    save(['nets_',EEG_list(i).name],'networks_13_20_PLV');
    
    %save wPLI
    cd(path2output_wPLI_13_20)
    save(['nets_',EEG_list(i).name],'networks_13_20_wPLI');
    
    %% 20-32 Hz
    
    networks_20_32_PLV = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_20_32_wPLI = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    
    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,20,sf,250,'high');
        PLV=plv_undirected_with_zero_lag_no_surrogates(mysegmentsignal); % calculate PLV for each 20 second segment.
        %no surrogates used | padding in data removed after doing the hilbert transform
        networks_20_32_PLV(:,:,count_segments) = PLV;
        
        wPLI=weighted_phase_lag_index_padding(mysegmentsignal); % calculate PLV for each 20 second segment.
        networks_20_32_wPLI(:,:,count_segments) = wPLI; %calculate the PLV for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform
    end
    
    %save PLV
    cd(path2output_PLV_20_32)
    save(['nets_',EEG_list(i).name],'networks_20_32_PLV');
    
    %save wPLI
    cd(path2output_wPLI_20_32)
    save(['nets_',EEG_list(i).name],'networks_20_32_wPLI');
end
