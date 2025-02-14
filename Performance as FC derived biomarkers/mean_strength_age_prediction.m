% produces table of mean strength for each FC method, in each frequency band

% this code uses strengths_und.m from the brain connectivity toolbox - https://sites.google.com/site/bctnet/

% input: FC adjacency matrices for all FC methods, in each freq band
% output: table with mean strength for all FC methods, in each freq band

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
MainEEG = '...EEG_networks_DK_1200_methodspaper/full_nets'; % define a main directory
EEG_list = dir(fullfile(MainEEG,'**/nets_ROIdata*'));

output = nan*ones(length(network_file_list),31); % 31 columns: 6 methods x 5 frequency bands + ID
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

    %mean strength code
    [str] = strengths_und(connectivity_matrix);
    ms = mean(str);

    %store based on FC method and frequency bin
    %outputs per frequency per FC method

    FC_freq = erase(string(EEG_list(j).folder), .../EEG_networks_DK_1200_methodspaper/full_nets/");

    if FC_freq == "COH/1_4Hz"
        output(1, 2) = ms;

    elseif FC_freq == "COH/4_8Hz"
        output(1,3) = ms;

    elseif FC_freq == "COH/8_13Hz"
        output(1,4) = ms;

    elseif FC_freq == "COH/13_20Hz"
        output(1,5) = ms;

    elseif FC_freq == "COH/20_32Hz"
        output(1,6) = ms;

    elseif FC_freq == "img_COH/1_4Hz"
        output(1,7) = ms;

    elseif FC_freq == "img_COH/4_8Hz"
        output(1,8) = ms;

    elseif FC_freq == "img_COH/8_13Hz"
        output(1,9) = ms;

    elseif FC_freq == "img_COH/13_20Hz"
        output(1,10) = ms;

    elseif FC_freq == "img_COH/20_32Hz"
        output(1,11) = ms;

    elseif FC_freq == "PLV/1_4Hz"
        output(1,12) = ms;

    elseif FC_freq == "PLV/4_8Hz"
        output(1,13) = ms;

    elseif FC_freq == "PLV/8_13Hz"
        output(1,14) = ms;

    elseif FC_freq == "PLV/13_20Hz"
        output(1,15) = ms;

    elseif FC_freq == "PLV/20_32Hz"
        output(1,16) = ms;

    elseif FC_freq == "wPLI/1_4Hz"
        output(1,17) = ms;

    elseif FC_freq == "wPLI/4_8Hz"
        output(1,18) = ms;

    elseif FC_freq == "wPLI/8_13Hz"
        output(1,19) = ms;

    elseif FC_freq == "wPLI/13_20Hz"
        output(1,20) = ms;

    elseif FC_freq == "wPLI/20_32Hz"
        output(1,21) = ms;

    elseif FC_freq == "AEC/1_4Hz"
        output(1,22) = ms;

    elseif FC_freq == "AEC/4_8Hz"
        output(1,23) = ms;

    elseif FC_freq == "AEC/8_13Hz"
        output(1,24) = ms;

    elseif FC_freq == "AEC/13_20Hz"
        output(1,25) = ms;

    elseif FC_freq == "AEC/20_32Hz"
        output(1,26) = ms;

    elseif FC_freq == "Orth_AEC/1_4Hz"
        output(1,27) = ms;

    elseif FC_freq == "Orth_AEC/4_8Hz"
        output(1,28) = ms;

    elseif FC_freq == "Orth_AEC/8_13Hz"
        output(1,29) = ms;

    elseif FC_freq == "Orth_AEC/13_20Hz"
        output(1,30) = ms;

    elseif FC_freq == "Orth_AEC/20_32Hz"
        output(1,31) = ms;

    end

    %make the file name into a number
    pre_ID = erase(string(network_file_list(i).name), ".mat");
    ID = str2double(erase(string(pre_ID), "ROIdata "));
    output(i,1) = ID;

end

