% produces normalised, weighted shortest pathlength from adjacency matrices
% for each FC method, in 4-8 Hz

% this code uses weight_conversion.m, distance_wei.m, charpath.m,
% null_model_und_sign.m from the brain connectivity toolbox - https://sites.google.com/site/bctnet/

% input: adjacency matrices with SWM specific ROIs, in 4-8 Hz
% output: table with normalised, weighted shortest pathlength for all EEG
% FC methods

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript: 
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------

%% clear variables and add paths as required
clear
clc

%set up paths DEPENDING ON FREQUENCY BAND
path2nets = ('...'); %load networks
addpath('...') %add path to scripts etc
path2output = '...';

%load data
network_files = fullfile(path2nets, '**', 'ROIdata*'); %input adjacency matrices for each FC method, 4-8Hz
network_file_list = dir(network_files); %directory lists files and folders in the current folder

output = nan*ones(length(network_file_list),7); % 7 columns: 6 methods + ID
cd(path2nets)

%%
%run loop
for i = 1:length(network_file_list)
    %%
    disp(i)

    data = load(network_file_list(i).name);
    fnames = fieldnames(data);
    EEG_network = data.(fnames{1});

    %get rid of the temporarily needed variables
    clear data
    clear fnames

    %Ensure only the first network is chosen.
    connectivity_matrix = EEG_network(:,:,1); %only keep the first, get rid of the rest

    %path length code
    temp_matrix = weight_conversion(connectivity_matrix, 'lengths'); %weights into distances by inversion method
    [D,~] = distance_wei(temp_matrix); %calculates the distance between all pairs of nodes
    [lambda,~,~,~,~] = charpath(D); %calculates the average distance between all pairs of nodes

    % calculate CC and PL for surrogate networks
    surrogate_matrix = nan*ones(500,1);

    for k = 1:500

        %make null models
        [W0,~] = null_model_und_sign(connectivity_matrix,5,0.5); %Random graphs with preserved weight, degree and
        %strength distributions

        %path length stuff
        temp_matrix_sug = weight_conversion(W0, 'lengths'); %weights into distances by inversion method
        [D_sug,~] = distance_wei(temp_matrix_sug); %calculates the distance between all pairs of nodes
        [lambda_sug,~,~,~,~] = charpath(D_sug,0,0); %calculates the average distance between all pairs of nodes

        % save surrogate PL into surrogate matrix
        surrogate_matrix(k,1) = lambda_sug;

    end

    % weighted normalised pathlength = pathlength / pathlength of null models

    nPL = (lambda/nanmean(surrogate_matrix(:,1)));

    %store based on FC method and frequency bin
    %outputs per frequency per FC method

    FC_freq = erase(string(netfiles(i).folder), ".../SWM_matrices/");

    if FC_freq == "COH/4_8Hz"
        output(1,2) = nPL;

    elseif FC_freq == "img_COH/4_8Hz"
        output(1,3) = nPL;

    elseif FC_freq == "PLV/4_8Hz"
        output(1,4) = nPL;

    elseif FC_freq == "wPLI/4_8Hz"
        output(1,5) = nPL;

    elseif FC_freq == "AEC/4_8Hz"
        output(1,6) = nPL;

    elseif FC_freq == "Orth_AEC/4_8Hz"
        output(1,7) = nPL;

    end

    %make the file name into a number
    pre_ID = erase(string(network_file_list(i).name), ".mat");
    ID = str2double(erase(string(pre_ID), "ROIdata "));
    output(i,1) = ID;

end

