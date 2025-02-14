% this script produces streamline count, proxy for WM signal transmission time,
% penalisability of phase delay (penalisation scale), in preparation for post hoc mediation
% analysis

% requirements: DWI matrices, phase lag matrices
% output: table with structural connectivity streamline count, proxy for WM signal transmission time,
% penalisability of phase delay (penalisation scale)

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript:
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------

%% clear variables and add paths as required
clear
clc

addpath('...') %add path to scripts

%LAG PATHS AND VALUES
path2lag = '...EEG_networks_DK_1200_methodspaper/pos_neg_lag';
lag_files = fullfile(path2lag, '**/', '*.mat');
lag_list = dir(lag_files);

% make an ID variable. I don't know how to do this in a structure, so I
% convert structure --> table and then table --> structure
lag_list = struct2table(lag_list);
lag_list.ID = lag_list.name;
lag_list.ID = erase(string(lag_list.ID), '.mat');
lag_list.ID = erase(string(lag_list.ID), 'nets_ROIdata ');
lag_list = sortrows(lag_list, {'name', 'folder'});
lag_list = table2struct(lag_list);

%DIFFUSION paths
MainDiffusion = '.../detailed_structural_metrics';
Diffusion_list = dir(fullfile(MainDiffusion,'**/', '*.mat'));

% make an ID variable. I don't know how to do this in a structure, so I
% convert structure --> table and then table --> structure
Diffusion_list = struct2table(Diffusion_list);
Diffusion_list.ID = Diffusion_list.name;
Diffusion_list.ID = erase(string(Diffusion_list.ID), '_hardiSD_ABS0.002_ang45_type2_pconn_d3.00_s0.00.mat');
Diffusion_list = sortrows(Diffusion_list, {'name'});
Diffusion_list = table2struct(Diffusion_list);

%%
path2output = '...';

% Diffusion imaging participants only
%%
output = [];

for i = 1:length(Diffusion_list) %go to the first lag network
    disp(i)
    %%
    %make an ID number
    Diffusion_ID = Diffusion_list(i).ID;

    % load diffusion image
    cd(Diffusion_list(i).folder)
    load(Diffusion_list(i).name);

    %% DIFFUSION IMAGING METRIC OF INTEREST

    % Xpln = Log-transformed, then normalised streamline count
    % Lp = Each edge represents the median length (in mm) of the streamlines connecting parcels i and j.
    % Hp = Each edge represents the median HMOA of the streamlines connecting parcels i and j.
    % (TO USE TRANSMISSION TIME, MUST SELECT HP)

    streamlines = Xp(triu(true(size(Xp)) ,1)); %this uses
    %non-normalised, non-log-transformed streamline count

    HMOA = Hp(triu(true(size(Hp)) ,1));
    distance = Lp(triu(true(size(Lp)) ,1));

    temp_per_participant = [];

end

for j = 1:length(lag_list)

    if eq(lag_list(j).ID,Diffusion_ID)

        cd(lag_list(j).folder);
        data = load(lag_list(j).name);
        fnames = fieldnames(data);
        lag_network = data.(fnames{1});

        %get rid of the temporarily needed variables
        clear data
        clear fnames

        lag_network = lag_network(:,:,1); %ensure you only have one lag matrix per participant

        %select only the TOP triange for each
        lag_vector = lag_network(triu(true(size(lag_network)) ,1));

        % Penalisation Scale
        % here, the most penalised (+-pi and 0) lags are given values of
        % 1, while lags of +- pi/2 are given values of 0.

        penalisation_value = (abs(abs(lag_vector) - pi/2))/(pi/2);

        % remove lag values where there is not a direct connection
        penalisation_value(isnan(HMOA)) = [];

        % sort by frequency band

        freq = erase(string(lag_list(j).folder), ".../EEG_networks_DK_1200_methodspaper/pos_neg_lag/");

        if freq == "1_4Hz"
            temp_per_participant(:,4) = penalisation_value;
        elseif freq == "4_8Hz"
            temp_per_participant(:,5) = penalisation_value;
        elseif freq == "8_13Hz"
            temp_per_participant(:,6) = penalisation_value;
        elseif freq == "13_20Hz"
            temp_per_participant(:,7) = penalisation_value;
        elseif freq == "20_32Hz"
            temp_per_participant(:,8) = penalisation_value;
        end

    end

end

streamlines(isnan(HMOA)) = [];
distance(isnan(HMOA)) = [];
HMOA(isnan(HMOA)) = [];

%TIME = DISTANCE / SPEED
time_cor = distance./HMOA;

temp_per_participant(:,1) = streamlines;
temp_per_participant(:,2) = time_cor;
temp_per_participant(:,3) = HMOA;
output = [output; temp_per_participant];

end

%display data as a table
Summary_data = array2table(output);
Summary_data.Properties.VariableNames = {'streamline_strength', 'time_proxy',	'HMOA',	'penalised_1_4', 'penalised_4_8', 'penalised_8_13',	'penalised_13_20',	'penalised_20_32'};

cd(path2output)
save_name = ...;
writetable(Summary_data,save_name)