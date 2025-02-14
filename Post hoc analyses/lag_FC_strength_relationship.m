%% script that looks at the relationship between phase lag and strength of functional connections

% produces graph of FC strength vs phase lag
% performs statistical analysis of FC strength of near-zero vs
% non-near-zero interactions between neural populations

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript:
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------

%% clear variables and add paths as required
clear
clc

addpath('...') %add path to scripts

%lag paths
path2lag = '.../EEG_networks_DK_1200_methodspaper/pos_neg_lag';
lag_files = fullfile(path2lag, '**/', '*.mat');
lag_list = dir(lag_files);

% make an ID variable. I don't know how to do this in a structure, so I
% convert structure --> table and then table --> structure
lag_list = struct2table(lag_list);
lag_list.ID = lag_list.name;
lag_list.ID = erase(string(lag_list.ID), '.mat');
lag_list.ID = convertStringsToChars(erase(string(lag_list.ID), 'nets_ROIdata '));
lag_list.freq = lag_list.folder;
lag_list.freq = erase(string(lag_list.freq), '.../EEG_networks_DK_1200_methodspaper/pos_neg_lag/');
lag_list = sortrows(lag_list, {'name', 'folder'});
lag_list = table2struct(lag_list);

%EEGs paths in the 'grand structure'
MainEEG = '.../EEG_networks_DK_1200_methodspaper/full_nets'; % define a main directory
EEG_list = dir(fullfile(MainEEG,'**/nets_ROIdata*'));

%sort rows of EEG
EEG_list = struct2table(EEG_list); % convert the struct array to a table
EEG_list.ID = EEG_list.name;
EEG_list.ID = erase(string(EEG_list.ID), ".mat");
EEG_list.ID = convertStringsToChars(erase(string(EEG_list.ID), "nets_ROIdata "));
EEG_list.freq = EEG_list.folder;
EEG_list.freq = erase(string(EEG_list.freq), '.../EEG_networks_DK_1200_methodspaper/full_nets/');
EEG_list.freq = extractAfter(EEG_list.freq,"/");

EEG_list = sortrows(EEG_list, {'name', 'folder'}); % sort by name followed by folder type
EEG_list = table2struct(EEG_list); % change it back to struct array if necessary

%% limit by Euclidian distance

%distance paths
euclidean = '.../EEG_networks_DK_1200_methodspaper/euclidean_distance/';
euclidean_list = dir(fullfile(euclidean,'*_distMat.mat'));
euclidean_list = struct2table(euclidean_list);
euclidean_list.euc_ID = erase(euclidean_list.name,'_distMat.mat');
euclidean_list.euc_ID = convertCharsToStrings(euclidean_list.euc_ID);
euclidean_list = table2struct(euclidean_list);

%% MAKE AVERAGE EUCLIDEAN DISTANCE MATRIX
%I have decided to do just one across all age groups, as the mean euclidean
%distance between ROIs only changes from 68.34mm from the youngest
%participant, to 72.4850 as the maxiumum value.

cd(euclidean)

euclidean_all = nan*ones(68,68,length(euclidean_list));

for w = 1:length(euclidean_list) %bring up each set of networks one at a time
    load(euclidean_list(w).name);
    euclidean_all(:,:,w) = D;
end

mean_euclidean = mean(euclidean_all,3);
euc_vector = mean_euclidean(triu(true(size(mean_euclidean)) ,1));

%% STATE WHICH PARTICIPANTS YOU WANT TO SELECT
selected_participants = '...';

%% output matrices
coh_1_4 = [];
coh_4_8 = [];
coh_8_13 = [];
coh_13_20 = [];
coh_20_32 = [];

img_coh_1_4 = [];
img_coh_4_8 = [];
img_coh_8_13 = [];
img_coh_13_20 = [];
img_coh_20_32 = [];

plv_1_4 = [];
plv_4_8 = [];
plv_8_13 = [];
plv_13_20 = [];
plv_20_32 = [];

wpli_1_4 = [];
wpli_4_8 = [];
wpli_8_13 = [];
wpli_13_20 = [];
wpli_20_32 = [];

aec_1_4 = [];
aec_4_8 = [];
aec_8_13 = [];
aec_13_20 = [];
aec_20_32 = [];

orth_aec_1_4 = [];
orth_aec_4_8 = [];
orth_aec_8_13 = [];
orth_aec_13_20 = [];
orth_aec_20_32 = [];

%% go through each EEG at a time
for i = 1:length(EEG_list) %go to the first lag network
    disp(i)
    %%
    cd(EEG_list(i).folder);
    data = load(EEG_list(i).name);
    fnames = fieldnames(data);
    EEG_network = data.(fnames{1});
    EEG_network = EEG_network(:,:,1); %only keep the first, get rid of the rest

    EEG_cor = EEG_network(triu(true(size(EEG_network)) ,1));
    EEG_cor(euc_vector < 80) = []; %remove region-pairs (FC strength) that are closer than 80mm Euclidean distance

    %get rid of the temporarily needed variables
    clear data
    clear fnames

    EEG_ID = EEG_list(i).ID;
    EEG_freq = EEG_list(i).freq;

    for j = 1:length(lag_list)

        % proceed if correct ID match and frequency match
        if eq(lag_list(j).ID,EEG_ID)
            if eq(lag_list(j).freq,EEG_freq)

                cd(lag_list(j).folder);
                data = load(lag_list(j).name);
                fnames = fieldnames(data);
                lag_network = data.(fnames{1});
                lag_network = lag_network(:,:,1); %only keep the first

                lag_cor = lag_network(triu(true(size(lag_network)) ,1));
                lag_cor = abs(lag_cor);
                lag_cor(euc_vector < 80) = []; %remove region-pairs (lags) that are closer than 80mm Euclidean distance

                %get rid of the temporarily needed variables
                clear data
                clear fnames

            end
        end
    end

    FC_lag_list = [EEG_cor lag_cor];

    FC_freq = erase(string(EEG_list(i).folder), ".../EEG_networks_DK_1200_methodspaper/full_nets/");

    if FC_freq == "COH/1_4Hz"
        coh_1_4 = [coh_1_4;FC_lag_list];

    elseif FC_freq == "COH/4_8Hz"
        coh_4_8 = [coh_4_8;FC_lag_list];

    elseif FC_freq == "COH/8_13Hz"
        coh_8_13 = [coh_8_13;FC_lag_list];

    elseif FC_freq == "COH/13_20Hz"
        coh_13_20 = [coh_13_20;FC_lag_list];

    elseif FC_freq == "COH/20_32Hz"
        coh_20_32 = [coh_20_32;FC_lag_list];

    elseif FC_freq == "img_COH/1_4Hz"
        img_coh_1_4 = [img_coh_1_4;FC_lag_list];

    elseif FC_freq == "img_COH/4_8Hz"
        img_coh_4_8 = [img_coh_4_8;FC_lag_list];

    elseif FC_freq == "img_COH/8_13Hz"
        img_coh_8_13 = [img_coh_8_13;FC_lag_list];

    elseif FC_freq == "img_COH/13_20Hz"
        img_coh_13_20 = [img_coh_13_20;FC_lag_list];

    elseif FC_freq == "img_COH/20_32Hz"
        img_coh_20_32 = [img_coh_20_32;FC_lag_list];

    elseif FC_freq == "PLV/1_4Hz"
        plv_1_4 = [plv_1_4;FC_lag_list];

    elseif FC_freq == "PLV/4_8Hz"
        plv_4_8 = [plv_4_8;FC_lag_list];

    elseif FC_freq == "PLV/8_13Hz"
        plv_8_13 = [plv_8_13;FC_lag_list];

    elseif FC_freq == "PLV/13_20Hz"
        plv_13_20 = [plv_13_20;FC_lag_list];

    elseif FC_freq == "PLV/20_32Hz"
        plv_20_32 = [plv_20_32;FC_lag_list];

    elseif FC_freq == "wPLI/1_4Hz"
        wpli_1_4 = [wpli_1_4;FC_lag_list];

    elseif FC_freq == "wPLI/4_8Hz"
        wpli_4_8 = [wpli_4_8;FC_lag_list];

    elseif FC_freq == "wPLI/8_13Hz"
        wpli_8_13 = [wpli_8_13;FC_lag_list];

    elseif FC_freq == "wPLI/13_20Hz"
        wpli_13_20 = [wpli_13_20;FC_lag_list];

    elseif FC_freq == "wPLI/20_32Hz"
        wpli_20_32 = [wpli_20_32;FC_lag_list];

    elseif FC_freq == "AEC/1_4Hz"
        aec_1_4 = [aec_1_4;FC_lag_list];

    elseif FC_freq == "AEC/4_8Hz"
        aec_4_8 = [aec_4_8;FC_lag_list];

    elseif FC_freq == "AEC/8_13Hz"
        aec_8_13 = [aec_8_13;FC_lag_list];

    elseif FC_freq == "AEC/13_20Hz"
        aec_13_20 = [aec_13_20;FC_lag_list];

    elseif FC_freq == "AEC/20_32Hz"
        aec_20_32 = [aec_20_32;FC_lag_list];

    elseif FC_freq == "orth_AEC/1_4Hz"
        orth_aec_1_4 = [orth_aec_1_4;FC_lag_list];

    elseif FC_freq == "orth_AEC/4_8Hz"
        orth_aec_4_8 = [orth_aec_4_8;FC_lag_list];

    elseif FC_freq == "orth_AEC/8_13Hz"
        orth_aec_8_13 = [orth_aec_8_13;FC_lag_list];

    elseif FC_freq == "orth_AEC/13_20Hz"
        orth_aec_13_20 = [orth_aec_13_20;FC_lag_list];

    elseif FC_freq == "orth_AEC/20_32Hz"
        orth_aec_20_32 = [orth_aec_20_32;FC_lag_list];

    end

end

%%
cor_matrices = {coh_1_4,	coh_4_8,	coh_8_13,	coh_13_20,	coh_20_32,	img_coh_1_4,	img_coh_4_8,	img_coh_8_13,	img_coh_13_20,	img_coh_20_32,	plv_1_4,	plv_4_8,	plv_8_13,	plv_13_20,	plv_20_32,	wpli_1_4,	wpli_4_8,	wpli_8_13,	wpli_13_20,	wpli_20_32,	aec_1_4,	aec_4_8,	aec_8_13,	aec_13_20,	aec_20_32,	orth_aec_1_4,	orth_aec_4_8,	orth_aec_8_13,	orth_aec_13_20,	orth_aec_20_32};
there_da_names = {	'coh_1_4',	'coh_4_8',	'coh_8_13',	'coh_13_20',	'coh_20_32',	'img_coh_1_4',	'img_coh_4_8',	'img_coh_8_13',	'img_coh_13_20',	'img_coh_20_32',	'plv_1_4',	'plv_4_8',	'plv_8_13',	'plv_13_20',	'plv_20_32',	'wpli_1_4',	'wpli_4_8',	'wpli_8_13',	'wpli_13_20',	'wpli_20_32',	'aec_1_4',	'aec_4_8',	'aec_8_13',	'aec_13_20',	'aec_20_32',	'orth_aec_1_4',	'orth_aec_4_8',	'orth_aec_8_13',	'orth_aec_13_20',	'orth_aec_20_32'	};
EEG_FC_name = {'coh',	'coh',	'coh',	'coh',	'coh',	'img_coh',	'img_coh',	'img_coh',	'img_coh',	'img_coh',	'plv',	'plv',	'plv',	'plv',	'plv',	'wpli',	'wpli',	'wpli',	'wpli',	'wpli',	'aec',	'aec',	'aec',	'aec',	'aec',	'orth_aec',	'orth_aec',	'orth_aec',	'orth_aec',	'orth_aec'	};
freq = {'A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E'};
lag_EEG_structure = struct('cor_matrices', cor_matrices, 'method_freq',there_da_names, 'freq_category', freq, 'method',EEG_FC_name);

%% PLOT ALL OF THE GRAPHS
T = struct2table(lag_EEG_structure); % convert the struct array to a table
sortedT = sortrows(T, {'freq_category'}); % sort by name followed by folder type
lag_EEG_structure = table2struct(sortedT); % change it back to struct array if necessary
clear T
clear sortedT

t = tiledlayout(5,6);
t.Padding = 'none';
t.TileSpacing = 'compact';

for z = 1:length(lag_EEG_structure);
    disp(z)
    graph = lag_EEG_structure(z).cor_matrices;
    absolute_lag = graph(:,2);
    EEG_FC_value = graph(:,1);

    nexttile

    C = ksdensity([absolute_lag, EEG_FC_value],[absolute_lag, EEG_FC_value]);
    scatter(absolute_lag, EEG_FC_value,10,C, 'filled');
    colorbar()
    ylim([0 1]);
    % colormap(parula)
    % title(lag_EEG_structure(z).method_freq);
    % xlabel('Absolute Phase Lag');
    % ylabel('EEG FC Value');
    % caxis([0 37]);

end


%%  PLOT ONLY ALPHA SO THAT THEY ARE MORE CLEARLY VISIBLE

% rearrange structure so that it's plotted properly
T = struct2table(lag_EEG_structure); % convert the struct array to a table
sortedT = sortrows(T, {'freq_category'}); % sort by name followed by folder type
lag_EEG_structure = table2struct(sortedT); % change it back to struct array if necessary
clear T
clear sortedT

alpha_lag_structure = lag_EEG_structure;

T = struct2table(alpha_lag_structure); % convert the struct array to a table
T = T(13:18,:);
reorder = {'A','D','B','E','C','F'}; %order so that the dirty methods are on top and the clean methods are below

for c = 1:length(reorder)
    T.freq_category(c) = reorder(c);
end

T = sortrows(T, {'freq_category'});
alpha_lag_structure = table2struct(T); % change it back to struct array if necessary
clear T

t = tiledlayout(2,3);
t.Padding = 'none';
t.TileSpacing = 'compact';

alpha_lag_structure(1).method = 'COH';
alpha_lag_structure(2).method = 'PLV';
alpha_lag_structure(3).method = 'AEC';
alpha_lag_structure(4).method = 'Img COH';
alpha_lag_structure(5).method = 'wPLI';
alpha_lag_structure(6).method = 'Orth AEC';

for z = 1:length(alpha_lag_structure);
    disp(z)
    graph = alpha_lag_structure(z).cor_matrices;
    absolute_lag = graph(:,2);
    EEG_FC_value = graph(:,1);

    nexttile

    C = ksdensity([absolute_lag, EEG_FC_value],[absolute_lag, EEG_FC_value]);
    scatter(absolute_lag, EEG_FC_value,12,C, 'filled')
    colorbar('FontSize',16)
    %title(alpha_lag_structure(z).method,'FontSize',32,'Interpreter','none')


    %xlabel('Absolute Phase Lag','FontSize',16);
    xlim([0 pi]);
    ylim([0 1])
    %ylabel('EEG FC Value','FontSize',16);
    set(gca,'FontSize',16)


end

%% compare connectivity strength for near-zero versus non-near-zero phase-delay connectivity

% assign each phase value to near zero or not near zero.
% near zero group = 1, not-near zero group = 0;

nz = 0.3; %CAN TOGGLE DEFINITION

T = struct2table(lag_EEG_structure); % convert the struct array to a table
T.near_zero_connectivity_vs_not = T.cor_matrices;
lag_EEG_structure = table2struct(T);
clear T

for i = 1:length(lag_EEG_structure)
    X = lag_EEG_structure(i).near_zero_connectivity_vs_not;
    X = array2table(X, 'VariableNames',{'strength','phase'});
    mask = X.phase > - nz & X.phase < nz | X.phase > (pi-nz) | X.phase<(-pi+nz);
    X.phase_group(mask) = 1;
    lag_EEG_structure(i).near_zero_connectivity_vs_not = X;
end

clear X

% mean connectivity strength and standard deviation followed by t test
% you can report the results in another table

for i = 1:length(lag_EEG_structure)
    X = lag_EEG_structure(i).near_zero_connectivity_vs_not;
    thing = groupsummary(X,"phase_group",["max", "mean", "std","median", "min"], "strength");
    near_zero = X.strength(X.phase_group==1);
    not_near_zero = X.strength(X.phase_group==0);
    [p,h] = ranksum(near_zero,not_near_zero);
    thing(:,8) = {p*30};
    thing(:,9) = {lag_EEG_structure(i).method_freq};

    summary_t_tests(2*i-1:2*i,:) = thing;

end

path2output = '...';
cd(path2output)
save_name = 't_tests_penalisable.xlsx';
writetable(summary_t_tests,save_name)