%% script to correlate DWI matrix with EEG matrices with each method in each frequency band
% also contains code to obtrain WM tract distance and proxy of WM transmission time for post hoc
% analysis

% outputs of this script can be toggled to give the following:
% 1) per participant structure - function concordance for all FC methods in
% all frequency bands
% 2) average WM tract distance between region pairs per person (i.e. one average per pair of regions, as many tracts join one pair of regions) 
% 3) proxy of WM transmission time for each region-pair per person

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

%% DWI paths
path2DWI = '...DWI_methods_paper/included_study/'; %participants from 2 sites
DWI_files = fullfile(path2DWI, '**', '*.mat');
DWI_list = dir(DWI_files);

%%
%EEGs paths 
MainEEG = '...EEG_networks_DK_1200_methodspaper/full_nets'; % define a main directory
EEG_list = dir(fullfile(MainEEG,'**/nets_ROIdata*'));

%sort rows of EEG
EEG_list = struct2table(EEG_list); % convert the struct array to a table
EEG_list = sortrows(EEG_list, {'name', 'folder'}); % sort by name followed by folder type
EEG_list = table2struct(EEG_list); % change it back to struct array if necessary

path2output = '...';
output = nan*ones(50,31);
%last row of output is for the median

%%
for i = 1:length(DWI_list) %bring up each set of networks one at a time
    %%
    disp(i)
    cd(DWI_list(i).folder);
    load(DWI_list(i).name);

    %% DIFFUSION IMAGING METRIC OF INTEREST

    % Xpln = Log-transformed, then normalised streamline count
    % Lp = Each edge represents the median length (in mm) of the streamlines connecting parcels i and j.
    % Hp = Each edge represents the median HMOA of the streamlines connecting parcels i and j.
    % (TO USE TRANSMISSION TIME, MUST SELECT HP)

    diffusion = Xpln; %TOGGLE to choose diffusion imaging metric

    %%
    %make an ID number
    DWI_ID = str2num(DWI_list(i).ID);
    output(i,1) = DWI_ID;

    % then go through each EEG file and see if it matches
    for j = 1:length(EEG_list)
        pre_ID = erase(string(EEG_list(j).name), ".mat");
        EEG_ID = str2double(erase(string(pre_ID), "nets_ROIdata "));

        if eq(DWI_ID,EEG_ID)

            %if it matches, load the EEG
            cd(EEG_list(j).folder);
            data = load(EEG_list(j).name);
            fnames = fieldnames(data);
            EEG_network = data.(fnames{1});

            %get rid of the temporarily needed variables
            clear data
            clear fnames

            %Ensure only the first network is chosen.
            EEG_network = EEG_network(:,:,1); %only keep the first, get rid of the rest

            %select only the TOP triange for each
            EEG_cor = EEG_network(triu(true(size(EEG_network)) ,1));
            DWI_cor = diffusion(triu(true(size(diffusion)) ,1));
            % length_cor = Lp(triu(true(size(Lp)) ,1)); %length needed for
            % proxy of WM transmission time analysis

            %IF YOU WANT TO KEEP ONLY THE 'MONOSYNAPTIC' CONNECTIONS
            % EEG_cor(DWI_cor == 0) = [];
            % DWI_cor(DWI_cor == 0) = [];

            EEG_cor(isnan(DWI_cor)) = [];
            % length_cor(isnan(DWI_cor)) = [];%length needed for
            % proxy of WM transmission time analysis
            DWI_cor(isnan(DWI_cor)) = [];

            % %% IF YOU ARE INTERESTED IN A PROXY OF TRANSMISSION TIME
            % %TIME = DISTANCE / SPEED
            %
            % %Set up distance vector
            % distance = Lp(triu(true(size(Lp)) ,1));
            % distance(isnan(distance)) = [];
            %
            % %time equation:
            % time_cor = distance./DWI_cor; %MUST HAVE SELECTED Hp as the
            % metric of interest on line 56
            %%
            FC_freq = erase(string(EEG_list(j).folder), ".../EEG_networks_DK_1200_methodspaper/full_nets/");

            %% Structure - function concordance calculation
            %correlate the two matrices in a spearmnan correlation
            SF_cor = corr(EEG_cor,DWI_cor, 'Type','Spearman');

            %% USE THIS SECTION IF YOU JUST WANT TO LOOK AT THE TOP X PERCENT OF STRUCTURAL CONNECTIONS

            %sort rows by DWI weights, largest at the top
            % both_mats = [DWI_cor EEG_cor];
            % both_mats_tab = array2table(both_mats);
            % both_mats_tab.Properties.VariableNames = {'DWI', 'EEG'};
            % both_mats_tab = sortrows(both_mats_tab, {'DWI'}, 'descend');
            %
            % top_x_percent = 10; %CHANGE - percentage of top DWI strenghts you want to keep
            % top_percent_number = floor(height(both_mats_tab(:,1))/(100/top_x_percent));
            % top_edges_DWI = both_mats_tab(1:top_percent_number,:);
            % top_edges_DWI = top_edges_DWI{:,:}; %array to matrix
            %
            % [SF_cor, pval]  = corr(top_edges_DWI, 'Type','Spearman');
            % SF_cor = SF_cor(1,2);
            % pval = pval(1,2);

            %%

            if FC_freq == "COH/1_4Hz"
                output(i,2) = SF_cor;

            elseif FC_freq == "COH/4_8Hz"
                output(i,3) = SF_cor;

            elseif FC_freq == "COH/8_13Hz"
                output(i,4) = SF_cor;

            elseif FC_freq == "COH/13_20Hz"
                output(i,5) = SF_cor;

            elseif FC_freq == "COH/20_32Hz"
                output(i,6) = SF_cor;

            elseif FC_freq == "img_COH/1_4Hz"
                output(i,7) = SF_cor;

            elseif FC_freq == "img_COH/4_8Hz"
                output(i,8) = SF_cor;

            elseif FC_freq == "img_COH/8_13Hz"
                output(i,9) = SF_cor;

            elseif FC_freq == "img_COH/13_20Hz"
                output(i,10) = SF_cor;

            elseif FC_freq == "img_COH/20_32Hz"
                output(i,11) = SF_cor;

            elseif FC_freq == "PLV/1_4Hz"
                output(i,12) = SF_cor;

            elseif FC_freq == "PLV/4_8Hz"
                output(i,13) = SF_cor;

            elseif FC_freq == "PLV/8_13Hz"
                output(i,14) = SF_cor;

            elseif FC_freq == "PLV/13_20Hz"
                output(i,15) = SF_cor;

            elseif FC_freq == "PLV/20_32Hz"
                output(i,16) = SF_cor;

            elseif FC_freq == "wPLI/1_4Hz"
                output(i,17) = SF_cor;

            elseif FC_freq == "wPLI/4_8Hz"
                output(i,18) = SF_cor;

            elseif FC_freq == "wPLI/8_13Hz"
                output(i,19) = SF_cor;

            elseif FC_freq == "wPLI/13_20Hz"
                output(i,20) = SF_cor;

            elseif FC_freq == "wPLI/20_32Hz"
                output(i,21) = SF_cor;

            elseif FC_freq == "AEC/1_4Hz"
                output(i,22) = SF_cor;

            elseif FC_freq == "AEC/4_8Hz"
                output(i,23) = SF_cor;

            elseif FC_freq == "AEC/8_13Hz"
                output(i,24) = SF_cor;

            elseif FC_freq == "AEC/13_20Hz"
                output(i,25) = SF_cor;

            elseif FC_freq == "AEC/20_32Hz"
                output(i,26) = SF_cor;

            elseif FC_freq == "orth_AEC/1_4Hz"
                output(i,27) = SF_cor;

            elseif FC_freq == "orth_AEC/4_8Hz"
                output(i,28) = SF_cor;

            elseif FC_freq == "orth_AEC/8_13Hz"
                output(i,29) = SF_cor;

            elseif FC_freq == "orth_AEC/13_20Hz"
                output(i,30) = SF_cor;

            elseif FC_freq == "orth_AEC/20_32Hz"
                output(i,31) = SF_cor;

            end

        end
    end
end

%%
%display data as a table
Summary_data = array2table(output);