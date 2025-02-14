%% code to calculate the functional network as dervied from coherence and imaginary coherency

% This script takes pre-processed EEG data and produces functional
% connectivity matrices

%decision to use bandpass filtering with this too, so that you have the FC
%methods acting on all data points

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript:
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------

%% clean and set up drives
clear
clc

%set up folders
path2data = '/home/.../EEG_files/EEG_cropped_T1_DK_1200';
addpath([path2data,'/'])
addpath('/home/.../EEG_files/Script/FunNets/')

path2output_COH_1_4 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/COH/1_4Hz';
path2output_COH_4_8 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/COH/4_8Hz';
path2output_COH_8_13 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/COH/8_13Hz';
path2output_COH_13_20 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/COH/13_20Hz';
path2output_COH_20_32 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/COH/20_32Hz';

path2output_img_COH_1_4 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/img_COH/1_4Hz';
path2output_img_COH_4_8 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/img_COH/4_8Hz';
path2output_img_COH_8_13 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/img_COH/8_13Hz';
path2output_img_COH_13_20 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/img_COH/13_20Hz';
path2output_img_COH_20_32 = '/home/.../EEG_files/EEG_networks_DK_1200_methodspaper/img_COH/20_32Hz';

ROI_files = fullfile(path2data, '/', 'ROIdata*');
EEG_list = dir(ROI_files); %directory lists files and folders in the current folder
%%
%loop through each cropped EEG file
for i = 1:length(EEG_list)
    disp(i)
    load(EEG_list(i).name);

    % --- Ensure that each row is a EEG recording ---
    % if the rows are longer than the columns --> flip them around
    if size(ROIdata,1)>size(ROIdata,2)
        ROIdata=ROIdata';
    end
    %%
    % --- Cut data in 20 segments ---
    T_min=4;    % minimum time required for a segment
    T_seg = 24; % length in seconds for each segment
    [ROIdata,n_seg]=cut_segments(ROIdata,sf,T_seg,T_min);
    %n_seg tells you how many segments of 20 seconds
    %this cuts your segments into the full recording into 20 seconds.
    % returns ROI data as a 3D matrix. each segment is in its own dimension.

    %% 1-4 Hz
    % --- Calculate functional network for each 20 second segment one at a time, using COHERENCE and imaginary coherence ---

    networks_1_4_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_1_4_img_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,4,sf,250,'low'); %creates filtered epochs with padding
        mysegmentsignal = mysegmentsignal(:, 500:5500-1); %remove padding after bandpass filtering
        [network_coh,network_img_coh]=coherence_img_coherence(mysegmentsignal);

        networks_1_4_coh(:,:,count_segments) = network_coh;
        networks_1_4_img_coh(:,:,count_segments) = network_img_coh;

    end

    %save COH
    cd(path2output_COH_1_4)
    save(['nets_',EEG_list(i).name],'networks_1_4_coh');

    %save img_coh
    cd(path2output_img_COH_1_4)
    save(['nets_',EEG_list(i).name],'networks_1_4_img_coh');


    %% 4-8 Hz
    % --- Calculate functional network for each 20 second segment one at a time, using COHERENCE and imaginary coherence ---

    networks_4_8_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_4_8_img_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,4,sf,250,'high'); %creates filtered epochs with padding
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,8,sf,250,'low'); %creates filtered epochs with padding
        mysegmentsignal = mysegmentsignal(:, 500:5500-1); %remove padding after bandpass filtering
        [network_coh,network_img_coh]=coherence_img_coherence(mysegmentsignal);

        networks_4_8_coh(:,:,count_segments) = network_coh;
        networks_4_8_img_coh(:,:,count_segments) = network_img_coh;

    end

    %save COH
    cd(path2output_COH_4_8)
    save(['nets_',EEG_list(i).name],'networks_4_8_coh');

    %save img_coh
    cd(path2output_img_COH_4_8)
    save(['nets_',EEG_list(i).name],'networks_4_8_img_coh');

    %% 8-13 Hz
    % --- Calculate functional network for each 20 second segment one at a time, using COHERENCE and imaginary coherence ---

    networks_8_13_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_8_13_img_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,8,sf,250,'high'); %creates filtered epochs with padding
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,13,sf,250,'low'); %creates filtered epochs with padding
        mysegmentsignal = mysegmentsignal(:, 500:5500-1); %remove padding after bandpass filtering
        [network_coh,network_img_coh]=coherence_img_coherence(mysegmentsignal);

        networks_8_13_coh(:,:,count_segments) = network_coh;
        networks_8_13_img_coh(:,:,count_segments) = network_img_coh;

    end

    %save COH
    cd(path2output_COH_8_13)
    save(['nets_',EEG_list(i).name],'networks_8_13_coh');

    %save img_coh
    cd(path2output_img_COH_8_13)
    save(['nets_',EEG_list(i).name],'networks_8_13_img_coh');

    %% 13-20 Hz
    % --- Calculate functional network for each 20 second segment one at a time, using COHERENCE and imaginary coherence ---

    networks_13_20_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_13_20_img_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,13,sf,250,'high'); %creates filtered epochs with padding
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,20,sf,250,'low'); %creates filtered epochs with padding
        mysegmentsignal = mysegmentsignal(:, 500:5500-1); %remove padding after bandpass filtering
        [network_coh,network_img_coh]=coherence_img_coherence(mysegmentsignal);

        networks_13_20_coh(:,:,count_segments) = network_coh;
        networks_13_20_img_coh(:,:,count_segments) = network_img_coh;

    end

    %save COH
    cd(path2output_COH_13_20)
    save(['nets_',EEG_list(i).name],'networks_13_20_coh');

    %save img_coh
    cd(path2output_img_COH_13_20)
    save(['nets_',EEG_list(i).name],'networks_13_20_img_coh');

    %% 20-32 Hz
    % --- Calculate functional network for each 20 second segment one at a time, using COHERENCE and imaginary coherence ---

    networks_20_32_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);
    networks_20_32_img_coh = nan*ones(size(ROIdata,1),size(ROIdata,1),n_seg);

    for count_segments = 1:n_seg
        mysegmentsignal = squeeze(ROIdata(:,:,count_segments)); %squeeze makes 3D into 2D by removing dimensions of length 1
        mysegmentsignal=FIR1_filtfilt(mysegmentsignal,20,sf,250,'high'); %creates filtered epochs with padding
        mysegmentsignal = mysegmentsignal(:, 500:5500-1); %remove padding after bandpass filtering
        [network_coh,network_img_coh]=coherence_img_coherence(mysegmentsignal);

        networks_20_32_coh(:,:,count_segments) = network_coh;
        networks_20_32_img_coh(:,:,count_segments) = network_img_coh;

    end

    %save COH
    cd(path2output_COH_20_32)
    save(['nets_',EEG_list(i).name],'networks_20_32_coh');

    %save img_coh
    cd(path2output_img_COH_20_32)
    save(['nets_',EEG_list(i).name],'networks_20_32_img_coh');
end