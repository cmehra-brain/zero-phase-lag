%This script produces histograms of phase distribution between -pi and +pi
%radians. The histogram combines FC across all subjects.

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript: 
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------

%% clear variables and add paths as required
clear
clc

addpath('...') %add path to scripts

%% Define participants

all_ages = {...}; %cell containing participant IDs included in analysis

%% STATE WHICH PARTICIPANTS YOU WANT TO SELECT FOR THIS ANALYSIS
selected_participants = all_ages;
%%
%EEGs paths - define and load
MainEEG = '.../pos_neg_lag/'; %load matrices of average phase lag between all ROIs per participant. Matrices in each frequency band stored in separate folders
EEG_list = dir(fullfile(MainEEG,'**/nets_ROIdata*'));

%distance paths
euclidean = '.../euclidean_distance/'; %Euc distance derived from DWI images, n = 50

%load distance paths
euclidean_list = dir(fullfile(euclidean,'*_distMat.mat'));

%edit names
euclidean_list = struct2table(euclidean_list);
euclidean_list.euc_ID = erase(euclidean_list.name,'_distMat.mat');
euclidean_list.euc_ID = convertCharsToStrings(euclidean_list.euc_ID);
euclidean_list = table2struct(euclidean_list);

%% Create variables for outputs for each frequency band
lag_1_4 = [];
lag_4_8 = [];
lag_8_13 = [];
lag_13_20 = [];
lag_20_32 = [];

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

%% Make various histograms included in the paper

for i = 1:length(selected_participants)
    %disp(i)
    ID_of_interest = convertCharsToStrings(selected_participants{i});

    for j = 1:length(EEG_list)
        pre_ID = erase(string(EEG_list(j).name), ".mat");
        EEG_ID = erase(string(pre_ID), "nets_ROIdata ");

        if eq(ID_of_interest,EEG_ID)

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

            %% Toggle below based on the region-pairs you want included in analysis

            %% ONLY THE INTERHEMISPHEREIC HOMOTOPIC CONNECTIONS
            % q = [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34];
            % p = [35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68];
            % idx=sub2ind(size(EEG_network),q,p);
            % out=EEG_network(idx);
            % vector = out';

            %% INTERHEMISPHEREIC HETEROTOPIC CONNECTIONS
            % q = [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34];
            % p = [40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	35	36	37	38	39];
            % idx=sub2ind(size(EEG_network),q,p);
            % out=EEG_network(idx);
            % vector = out';

            %% ALL THE CONNECTIONS
            % vector = EEG_network(triu(true(size(EEG_network)) ,1));

            %% Only include FC at minimum Euclidean distance above (and including) 80mm
            % vector = EEG_network(triu(true(size(EEG_network)) ,1));
            % euc = 80;
            % vector(euc_vector < euc) = []; %needs to be less than, as you are INCLUDING region pairs with the distance stated in euc

            %%
            FC_freq = erase(string(EEG_list(j).folder), ".../pos_neg_lag"); %matrices of phase lag included in subfolders within this

            if FC_freq == "/1_4Hz"
                lag_1_4 = [lag_1_4; vector];

            elseif FC_freq == "/4_8Hz"
                lag_4_8 = [lag_4_8; vector];

            elseif FC_freq == "/8_13Hz"
                lag_8_13 = [lag_8_13; vector];

            elseif FC_freq == "/13_20Hz"
                lag_13_20 = [lag_13_20; vector];

            elseif FC_freq == "/20_32Hz"
                lag_20_32 = [lag_20_32; vector];

            end
        end
    end
end
%%
%put outputs from each frequency in a structure

mats = {lag_1_4,lag_4_8,lag_8_13,lag_13_20,	lag_20_32};
names_of_mats = {'lag_1_4',	'lag_4_8',	'lag_8_13',	'lag_13_20','lag_20_32'};
freq = {'A', 'B','C','D','E'};
dummies = {nan, nan, nan, nan, nan};
structure_mean_mat = struct('matrices', mats, 'method_freq',names_of_mats, 'mean_network' , dummies, 'freq_category', freq, 'near_zero', dummies, 'near_pi', dummies, 'near_neg_pi', dummies, 'near_penalizable', dummies);

for k = 1:length(structure_mean_mat)
    mat_method_freq = structure_mean_mat(k).matrices;
end

%% what percentage of functional connectivity is at zero or near zero phase delay:

nz = 0.3; %CAN TOGGLE DEFINITION

for a = 1:length(structure_mean_mat)
    X = structure_mean_mat(a).matrices;
    structure_mean_mat(a).near_zero = sum(X>=-nz & X<=nz)/length(X)*100;
    structure_mean_mat(a).near_pi = sum(X>=(pi-nz))/length(X)*100;
    structure_mean_mat(a).near_neg_pi = sum(X<=(-pi+nz))/length(X)*100;
    structure_mean_mat(a).near_penalizable = structure_mean_mat(a).near_zero + structure_mean_mat(a).near_pi + structure_mean_mat(a).near_neg_pi;
end

structure_mean_mat = struct2table(structure_mean_mat);

%% plot probability density function

structure_mean_mat = table2struct(structure_mean_mat);

for z = 1:length(structure_mean_mat)
    freq_interest = structure_mean_mat(z).matrices;
    subplot(1,5,z);
    histogram(freq_interest, 'BinWidth', 0.1)
    %ylim([0 20000]);
    title(structure_mean_mat(z).method_freq)
    %xlabel('Phase lag between ROIs (radians)');
end

