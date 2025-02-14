%% functional networks AEC + AEC orthogonalised methods paper
% This script takes pre-processed EEG data and produces functional
% connectivity matrices

% all these methods use the hilbert transform
% in 5 frequency bins
% 2 seconds of data padding either sides of 20 seconds resting state eyes
% closed data. This padding is removed once filtering and hilbert transform
% are done.

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript:
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------


%% clean and set up drives
clear
clc

%set up folders
addpath('/home/.../EEG_files/Script/FunNets/')

%input data
path2data = '/home/.../EEG_files/EEG_cropped_T1_DK_1200';
addpath([path2data,'/'])
ROI_files = fullfile(path2data, '/', 'ROIdata*');
EEG_list = dir(ROI_files); %directory lists files and folders in the current folder

%output data
path2output_AEC_1_4 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/AEC/1_4Hz';
path2output_AEC_4_8 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/AEC/4_8Hz';
path2output_AEC_8_13 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/AEC/8_13Hz';
path2output_AEC_13_20 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/AEC/13_20Hz';
path2output_AEC_20_32 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/AEC/20_32Hz';

path2output_orth_AEC_1_4 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/orth_AEC/1_4Hz';
path2output_orth_AEC_4_8 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/orth_AEC/4_8Hz';
path2output_orth_AEC_8_13 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/orth_AEC/8_13Hz';
path2output_orth_AEC_13_20 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/orth_AEC/13_20Hz';
path2output_orth_AEC_20_32 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/orth_AEC/20_32Hz';


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

    %% FILTERING AND CALCULATING AEC + orth_AEC IN ONE STEP 1-4Hz

    networks_1_4_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_1_4_orth_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal = FIR1_filtfilt(mysegmentsignal,4,sf,250,'low'); %creates filtered epochs of the right length
        AEC = AEC_connectivity(mysegmentsignal); % calculate AEC for each 20 second segment.
        %padding in data removed after doing the hilbert transform
        networks_1_4_AEC(:,:,count_segments) = AEC;

        orth_AEC = amplitude_envelope_correlation_orth(mysegmentsignal); % calculate AEC for each 20 second segment.
        networks_1_4_orth_AEC(:,:,count_segments) = orth_AEC; %calculate the AEC for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform

        %         AEC = AEC_connectivity(mysegmentsignal);
        %         networks_1_4_AEC(:,:,count_segments) = AEC;

    end

    %save AEC
    cd(path2output_AEC_1_4)
    save(['nets_',EEG_list(i).name],'networks_1_4_AEC');

    %save orth_AEC
    cd(path2output_orth_AEC_1_4)
    save(['nets_',EEG_list(i).name],'networks_1_4_orth_AEC');

    %% 4-8 Hz

    networks_4_8_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_4_8_orth_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,4,sf,250,'high');
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,8,sf,250,'low'); %creates filtered epochs of the right length
        AEC=AEC_connectivity(mysegmentsignal); % calculate AEC for each 20 second segment.
        %padding in data removed after doing the hilbert transform
        networks_4_8_AEC(:,:,count_segments) = AEC;

        orth_AEC=amplitude_envelope_correlation_orth(mysegmentsignal); % calculate AEC for each 20 second segment.
        networks_4_8_orth_AEC(:,:,count_segments) = orth_AEC; %calculate the AEC for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform
    end

    %save AEC
    cd(path2output_AEC_4_8)
    save(['nets_',EEG_list(i).name],'networks_4_8_AEC');

    %save orth_AEC
    cd(path2output_orth_AEC_4_8)
    save(['nets_',EEG_list(i).name],'networks_4_8_orth_AEC');

    %% 8-13 Hz

    networks_8_13_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_8_13_orth_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,8,sf,250,'high');
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,13,sf,250,'low'); %creates filtered epochs of the right length
        AEC=AEC_connectivity(mysegmentsignal); % calculate AEC for each 20 second segment.
        %padding in data removed after doing the hilbert transform
        networks_8_13_AEC(:,:,count_segments) = AEC;

        orth_AEC = amplitude_envelope_correlation_orth(mysegmentsignal); % calculate AEC for each 20 second segment.
        networks_8_13_orth_AEC(:,:,count_segments) = orth_AEC; %calculate the AEC for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform
    end

    %save AEC
    cd(path2output_AEC_8_13)
    save(['nets_',EEG_list(i).name],'networks_8_13_AEC');

    %save orth_AEC
    cd(path2output_orth_AEC_8_13)
    save(['nets_',EEG_list(i).name],'networks_8_13_orth_AEC');

    %% 13-20 Hz

    networks_13_20_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_13_20_orth_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,13,sf,250,'high');
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,20,sf,250,'low'); %creates filtered epochs of the right length
        AEC=AEC_connectivity(mysegmentsignal); % calculate AEC for each 20 second segment.
        %padding in data removed after doing the hilbert transform
        networks_13_20_AEC(:,:,count_segments) = AEC;

        orth_AEC = amplitude_envelope_correlation_orth(mysegmentsignal); % calculate AEC for each 20 second segment.
        networks_13_20_orth_AEC(:,:,count_segments) = orth_AEC; %calculate the AEC for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform
    end

    %save AEC
    cd(path2output_AEC_13_20)
    save(['nets_',EEG_list(i).name],'networks_13_20_AEC');

    %save orth_AEC
    cd(path2output_orth_AEC_13_20)
    save(['nets_',EEG_list(i).name],'networks_13_20_orth_AEC');

    %% 20-32 Hz

    networks_20_32_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_20_32_orth_AEC = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,20,sf,250,'high');
        AEC=AEC_connectivity(mysegmentsignal); % calculate AEC for each 20 second segment.
        %padding in data removed after doing the hilbert transform
        networks_20_32_AEC(:,:,count_segments) = AEC;

        orth_AEC=amplitude_envelope_correlation_orth(mysegmentsignal); % calculate AEC for each 20 second segment.
        networks_20_32_orth_AEC(:,:,count_segments) = orth_AEC; %calculate the AEC for each 20 second segment
        %this code would have already cropped out the 2 secs at the start
        %and at the end to adjust for hilbert transform
    end

    %save AEC
    cd(path2output_AEC_20_32)
    save(['nets_',EEG_list(i).name],'networks_20_32_AEC');

    %save orth_AEC
    cd(path2output_orth_AEC_20_32)
    save(['nets_',EEG_list(i).name],'networks_20_32_orth_AEC');
end
