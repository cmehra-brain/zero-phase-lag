%% This script produces table of proportion of near-zero-phase-lag FC vs Euclidean distance per participant in output

% data output will be in the following format
% Row 1: participant one, frequency band 1, proportions at each distance
% Row 2: participant one, frequency 2, proportions at each distance
% etc

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript:
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------

%%
clear
clc

%add paths to scripts
addpath('...')

%% load data

% lag data per participant
MainEEG = '.../pos_neg_lag/'; %load matrices of average phase lag between all ROIs per participant. Matrices in each frequency band stored in separate folders
EEG_list = dir(fullfile(MainEEG,'**/nets_ROIdata*'));

%distance paths
euclidean = '.../euclidean_distance/'; %Euc distance derived from DWI images, n = 50
euclidean_list = dir(fullfile(euclidean,'*_distMat.mat'));

%% Calculate average Euclidean distance between region pairs
%I have decided to do just one across all age groups, as the mean euclidean
%distance between ROIs only changes from 68.34mm from the youngest
%participant, to 72.4850 as the maxiumum value.
%because you are using average euclidean distance across participants, this
%allows you to extrapolate to participants without diffusion imaging data. therefore, all 153
%participants are included in this analysis.

cd(euclidean)

euclidean_all = nan*ones(68,68,length(euclidean_list));

for w = 1:length(euclidean_list) %bring up each set of networks one at a time
    load(euclidean_list(w).name);
    euclidean_all(:,:,w) = D;
end

mean_euclidean = mean(euclidean_all,3);
euc_vector = mean_euclidean(triu(true(size(mean_euclidean)) ,1));

%% specify the minimum Euclidean distances 
euc = [0 10 20 30 40 50 60 70 80 90];

%% toggle definition of what near zero is
nz = 0.3;
%% specify participant IDs which you want to include
selected_participants = {...};

%% loop through each participant

%make output variable
output = [];

for i = 1:length(selected_participants)
    disp(i)
    participant_ID = convertCharsToStrings(selected_participants{i});

    %make an empty table to store data per participant
    proportion_distance_per_freq = nan*ones(5,length(euc)+2);
    %add participant ID
    proportion_distance_per_freq(:,1) = participant_ID;
    %add label the frequency bands in order
    proportion_distance_per_freq(1,2) = "1";
    proportion_distance_per_freq(2,2) = "4";
    proportion_distance_per_freq(3,2) = "8";
    proportion_distance_per_freq(4,2) = "13";
    proportion_distance_per_freq(5,2) = "20";

    for j = 1:length(EEG_list)
        pre_ID = erase(string(EEG_list(j).name), ".mat");
        EEG_ID = erase(string(pre_ID), "nets_ROIdata ");

        if eq(participant_ID,EEG_ID)

            %if it matches, then act
            cd(EEG_list(j).folder);
            data = load(EEG_list(j).name);
            fnames = fieldnames(data);
            EEG_network = data.(fnames{1});

            %get rid of the temporarily needed variables
            clear data
            clear fnames

            %Ensure only the first network is chosen.
            EEG_network = EEG_network(:,:,1); %only keep the first, get rid of the rest

            % vectorise the matrix
            EEG_vector = EEG_network(triu(true(size(EEG_network)) ,1));

            %% only keep region pairs a certain distance appart

            proportion_distance = nan*ones(1,length(euc));

            for k = 1:length(euc)
                EEG_vector_per_distance = EEG_vector;
                EEG_vector_per_distance(euc_vector < euc(k)) = [];

                % then calculate the proportion of near-zero-phase delay at
                % each minimum euclidean distance

                near_zero = sum(EEG_vector_per_distance >= -nz & EEG_vector_per_distance <= nz)/length(EEG_vector_per_distance)*100;
                near_pi = sum(EEG_vector_per_distance>=(pi-nz))/length(EEG_vector_per_distance)*100;
                near_neg_pi = sum(EEG_vector_per_distance<=(-pi+nz))/length(EEG_vector_per_distance)*100;
                near_penalizable = near_zero + near_pi + near_neg_pi;

                %produces a vector of proportion of near-zero-phase delay
                %at each minimum euclidean distance

                proportion_distance(1, k) = near_penalizable;

            end

            %store it based on each frequency band

            FC_freq = erase(string(EEG_list(j).folder), ".../pos_neg_lag"); %matrices of phase lag included in subfolders within this

            if FC_freq == "/1_4Hz"
                proportion_distance_per_freq(1,3:end) = proportion_distance;

            elseif FC_freq == "/4_8Hz"
                proportion_distance_per_freq(2,3:end) = proportion_distance;

            elseif FC_freq == "/8_13Hz"
                proportion_distance_per_freq(3,3:end) = proportion_distance;

            elseif FC_freq == "/13_20Hz"
                proportion_distance_per_freq(4,3:end) = proportion_distance;

            elseif FC_freq == "/20_32Hz"
                proportion_distance_per_freq(5,3:end) = proportion_distance;
            end

        end
    end

    output = [output; proportion_distance_per_freq];
end
