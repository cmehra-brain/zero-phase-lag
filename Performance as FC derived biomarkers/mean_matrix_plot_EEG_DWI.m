% plots mean EEG and DWI matrices across participants

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

%%
coh_1_4 = nan*ones(68,68,length(EEG_list));
coh_4_8 = nan*ones(68,68,length(EEG_list));
coh_8_13 = nan*ones(68,68,length(EEG_list));
coh_13_20 = nan*ones(68,68,length(EEG_list));
coh_20_32 = nan*ones(68,68,length(EEG_list));

img_coh_1_4 = nan*ones(68,68,length(EEG_list));
img_coh_4_8 = nan*ones(68,68,length(EEG_list));
img_coh_8_13 = nan*ones(68,68,length(EEG_list));
img_coh_13_20 = nan*ones(68,68,length(EEG_list));
img_coh_20_32 = nan*ones(68,68,length(EEG_list));

PLV_1_4 = nan*ones(68,68,length(EEG_list));
PLV_4_8 = nan*ones(68,68,length(EEG_list));
PLV_8_13 = nan*ones(68,68,length(EEG_list));
PLV_13_20 = nan*ones(68,68,length(EEG_list));
PLV_20_32 = nan*ones(68,68,length(EEG_list));

wPLI_1_4 = nan*ones(68,68,length(EEG_list));
wPLI_4_8 = nan*ones(68,68,length(EEG_list));
wPLI_8_13 = nan*ones(68,68,length(EEG_list));
wPLI_13_20 = nan*ones(68,68,length(EEG_list));
wPLI_20_32 = nan*ones(68,68,length(EEG_list));

AEC_1_4 = nan*ones(68,68,length(EEG_list));
AEC_4_8 = nan*ones(68,68,length(EEG_list));
AEC_8_13 = nan*ones(68,68,length(EEG_list));
AEC_13_20 = nan*ones(68,68,length(EEG_list));
AEC_20_32 = nan*ones(68,68,length(EEG_list));

orth_AEC_1_4 = nan*ones(68,68,length(EEG_list));
orth_AEC_4_8 = nan*ones(68,68,length(EEG_list));
orth_AEC_8_13 = nan*ones(68,68,length(EEG_list));
orth_AEC_13_20 = nan*ones(68,68,length(EEG_list));
orth_AEC_20_32 = nan*ones(68,68,length(EEG_list));

%%
    for j = 1:length(EEG_list)

    cd(EEG_list(j).folder);
            data = load(EEG_list(j).name);
            fnames = fieldnames(data);
            EEG_network = data.(fnames{1});

            %get rid of the temporarily needed variables
            clear data
            clear fnames

            %Ensure only the first network is chosen.
            EEG_network = EEG_network(:,:,1); %only keep the first, get rid of the rest

            FC_freq = erase(string(EEG_list(j).folder), .../EEG_networks_DK_1200_methodspaper/full_nets/");

            %Note: These matrices will mainly consist of Nans
            if FC_freq == "COH/1_4Hz"
                coh_1_4(:,:,i) = EEG_network;

            elseif FC_freq == "COH/4_8Hz"
                coh_4_8(:,:,i) = EEG_network;

            elseif FC_freq == "COH/8_13Hz"
                coh_8_13(:,:,i) = EEG_network;

            elseif FC_freq == "COH/13_20Hz"
                coh_13_20(:,:,i) = EEG_network;

            elseif FC_freq == "COH/20_32Hz"
                coh_20_32(:,:,i) = EEG_network;

            elseif FC_freq == "img_COH/1_4Hz"
                img_coh_1_4(:,:,i) = EEG_network;

            elseif FC_freq == "img_COH/4_8Hz"
                img_coh_4_8(:,:,i) = EEG_network;

            elseif FC_freq == "img_COH/8_13Hz"
                img_coh_8_13(:,:,i) = EEG_network;

            elseif FC_freq == "img_COH/13_20Hz"
                img_coh_13_20(:,:,i) = EEG_network;

            elseif FC_freq == "img_COH/20_32Hz"
                img_coh_20_32(:,:,i) = EEG_network;

            elseif FC_freq == "PLV/1_4Hz"
                PLV_1_4(:,:,i) = EEG_network;

            elseif FC_freq == "PLV/4_8Hz"
                PLV_4_8(:,:,i) = EEG_network;

            elseif FC_freq == "PLV/8_13Hz"
                PLV_8_13(:,:,i) = EEG_network;

            elseif FC_freq == "PLV/13_20Hz"
                PLV_13_20(:,:,i) = EEG_network;

            elseif FC_freq == "PLV/20_32Hz"
                PLV_20_32(:,:,i) = EEG_network;

            elseif FC_freq == "wPLI/1_4Hz"
                wPLI_1_4(:,:,i) = EEG_network;

            elseif FC_freq == "wPLI/4_8Hz"
                wPLI_4_8(:,:,i) = EEG_network;

            elseif FC_freq == "wPLI/8_13Hz"
                wPLI_8_13(:,:,i) = EEG_network;

            elseif FC_freq == "wPLI/13_20Hz"
                wPLI_13_20(:,:,i) = EEG_network;

            elseif FC_freq == "wPLI/20_32Hz"
                wPLI_20_32(:,:,i) = EEG_network;

            elseif FC_freq == "AEC/1_4Hz"
                AEC_1_4(:,:,i) = EEG_network;

            elseif FC_freq == "AEC/4_8Hz"
                AEC_4_8(:,:,i) = EEG_network;

            elseif FC_freq == "AEC/8_13Hz"
                AEC_8_13(:,:,i) = EEG_network;

            elseif FC_freq == "AEC/13_20Hz"
                AEC_13_20(:,:,i) = EEG_network;

            elseif FC_freq == "AEC/20_32Hz"
                AEC_20_32(:,:,i) = EEG_network;

            elseif FC_freq == "orth_AEC/1_4Hz"
                orth_AEC_1_4(:,:,i) = EEG_network;

            elseif FC_freq == "orth_AEC/4_8Hz"
                orth_AEC_4_8(:,:,i) = EEG_network;

            elseif FC_freq == "orth_AEC/8_13Hz"
                orth_AEC_8_13(:,:,i) = EEG_network;

            elseif FC_freq == "orth_AEC/13_20Hz"
                orth_AEC_13_20(:,:,i) = EEG_network;

            elseif FC_freq == "orth_AEC/20_32Hz"
                orth_AEC_20_32(:,:,i) = EEG_network;

            end

        end

    end
end

%%
%put them all in a structure

mats = {coh_1_4,coh_4_8,coh_8_13,coh_13_20,	coh_20_32,	img_coh_1_4,img_coh_4_8,img_coh_8_13,	img_coh_13_20,	img_coh_20_32,	PLV_1_4,	PLV_4_8,	PLV_8_13,	PLV_13_20,	PLV_20_32,	wPLI_1_4,	wPLI_4_8,	wPLI_8_13,	wPLI_13_20,	wPLI_20_32, AEC_1_4,	AEC_4_8,	AEC_8_13,	AEC_13_20,	AEC_20_32,	orth_AEC_1_4,	orth_AEC_4_8,	orth_AEC_8_13,	orth_AEC_13_20,	orth_AEC_20_32};
names = {	'coh_1_4',	'coh_4_8',	'coh_8_13',	'coh_13_20',	'coh_20_32',	'img_coh_1_4',	'img_coh_4_8',	'img_coh_8_13',	'img_coh_13_20',	'img_coh_20_32',	'PLV_1_4',	'PLV_4_8',	'PLV_8_13',	'PLV_13_20',	'PLV_20_32',	'wPLI_1_4',	'wPLI_4_8',	'wPLI_8_13',	'wPLI_13_20',	'wPLI_20_32', 'AEC_1_4',	'AEC_4_8',	'AEC_8_13',	'AEC_13_20',	'AEC_20_32',	'orth_AEC_1_4',	'orth_AEC_4_8',	'orth_AEC_8_13','orth_AEC_13_20', 'orth_AEC_20_32'};
freq = {'A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E'};
dummies = {nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
structure_mean_mat = struct('matrices', mats, 'method_freq',names, 'mean_network' , dummies, 'freq_category', freq);

for k = 1:length(structure_mean_mat)
    mat_method_freq = structure_mean_mat(k).matrices;
    mean_net = mean(mat_method_freq,3);
    structure_mean_mat(k).mean_network = mean_net;
end

%%
%rearrange structure so that it's plotted properly
T = struct2table(structure_mean_mat); % convert the struct array to a table
sortedT = sortrows(T, 'freq_category'); % sort by name followed by folder type
structure_mean_mat = table2struct(sortedT); % change it back to struct array if necessary
clear T
clear sortedT

%% reorder matrices based on proximity
anatomical_order = [31	13	11	18	17	19	26	3	27	16	23	21	28	30	7	24	20	12	4	10	33	29	1	14	8	6	5	15	32	34	25	2	22	9	65	47	45	52	51	53	60	37	61	50	57	55	62	64	41	58	54	46	38	44	67	63	35	48	42	40	39	49	66	68	59	36	56	43];
for k = 1:length(structure_mean_mat)
    hi = structure_mean_mat(k).mean_network;
    hi = hi(anatomical_order,anatomical_order);
    structure_mean_mat(k).mean_network_reorder = hi;
end

%% produce all the adjacency matrices 

t = tiledlayout(5,6);
t.Padding = 'none';
t.TileSpacing = 'compact';

% for z = 1:length(structure_mean_mat)
% 
%     nexttile
%     imagesc(structure_mean_mat(z).mean_network_reorder, [0,1])
%     colorbar()
%     axis off
%     axis square
% 
%     % cm = flip([linspace(1,0,20)',repmat(linspace(0,1,20)',1,2); repmat(linspace(1,0,20)',1,2),linspace(0,1,20)']);
% 
% end

for z = 1:5

       nexttile
    imagesc(structure_mean_mat(6*z - 5).mean_network_reorder, [0,1])
    colorbar()
    axis off
    axis square

        nexttile
    imagesc(structure_mean_mat(6*z - 4).mean_network_reorder, [0,0.2])
    colorbar()
    axis off
    axis square

     nexttile
    imagesc(structure_mean_mat(6*z - 3).mean_network_reorder, [0,1])
    colorbar()
    axis off
    axis square

     nexttile
    imagesc(structure_mean_mat(6*z - 2).mean_network_reorder, [0,0.25])
    colorbar()
    axis off
    axis square

     nexttile
    imagesc(structure_mean_mat(6*z - 1).mean_network_reorder, [0,1])
    colorbar()
    axis off
    axis square

     nexttile
    imagesc(structure_mean_mat(6*z).mean_network_reorder, [0,0.14])
    colorbar()
    axis off
    axis square

end



%% Diffusion MEAN IMAGE
diffusion_files = fullfile(path2DWI, '/**', '*.mat');
diffusion_list = dir(diffusion_files);

DWI_all = nan*ones(68,68,length(diffusion_list));

% for streamline count
for w = 1:length(diffusion_list)
    cd(diffusion_list(w).folder)
    load(diffusion_list(w).name)
    DWI_all(:,:,w) = Xpln;

end 


mean_DWI = nanmean(DWI_all,3);
            DWI_vector = mean_DWI(triu(true(size(mean_DWI)) ,1));
            nanmedian(DWI_vector)

anatomical_order = [31	13	11	18	17	19	26	3	27	16	23	21	28	30	7	24	20	12	4	10	33	29	1	14	8	6	5	15	32	34	25	2	22	9	65	47	45	52	51	53	60	37	61	50	57	55	62	64	41	58	54	46	38	44	67	63	35	48	42	40	39	49	66	68	59	36	56	43];
mean_DWI = mean_DWI(anatomical_order,anatomical_order);

imagesc(mean_DWI)
axis square
axis off

%  % code to label :
%     tickLabelsP = {'F_L',	'F_L',	'F_L',	'F_L',	'F_L',	'F_L',	'F_L',	'F_L',	'F_L',	'F_L',	'F_L',	'P_L',	'P_L',	'P_L',	'P_L',	'P_L',	'O_L',	'O_L',	'O_L',	'O_L',	'T_L',	'T_L',	'T_L',	'T_L',	'T_L',	'T_L',	'T_L',	'T_L',	'T_L',	'I_L',	'C_L',	'C_L',	'C_L',	'C_L',...
%         'F_R',	'F_R',	'F_R',	'F_R',	'F_R',	'F_R',	'F_R',	'F_R',	'F_R',	'F_R',	'F_R',	'P_R',	'P_R',	'P_R',	'P_R',	'P_R',	'O_R',	'O_R',	'O_R',	'O_R',	'T_R',	'T_R',	'T_R',	'T_R',	'T_R',	'T_R',	'T_R',	'T_R',	'T_R',	'I_R',	'C_R',	'C_R',	'C_R',	'C_R'};
% 
%     uniqueLabels = unique(tickLabelsP,'stable');
%     numLabels    = numel(uniqueLabels);
%     labelBounds  = zeros(numLabels,2);
%     for L = 1:numLabels
%         parcelsThisLabel = find(strcmp(tickLabelsP,uniqueLabels{L}));
%         labelBounds(L,1) = parcelsThisLabel(1)-0.5;
%         labelBounds(L,2) = parcelsThisLabel(end)+0.5;
%     end
% 
%     labelCentre  = mean(labelBounds,2);
%     tickValues   = unique([labelBounds(:); labelCentre]);
% 
%     tickLabelsL  = cell(numel(tickValues),1);
%     for L = 1:numLabels
%         tickLabelsL(tickValues==labelCentre(L)) = uniqueLabels(L);
%     end
% 
%     xticks(tickValues)
%     yticks(tickValues)
%     xticklabels(tickLabelsL)
%     yticklabels(tickLabelsL)
% 
% 
%     set(gca,'TickLength',[0 0],'FontSize',12) %if you want to change the font size of the axes labels, do through here
%     for L = 2:numLabels
%         xline(labelBounds(L,1),'-k','LineWidth',1.5);
%         yline(labelBounds(L,1),'-k','LineWidth',1.5);
%     end
% %
% % make an ID variable. I don't know how to do this in a structure, so I
% % convert structure --> table and then table --> structure
% DWI_list = struct2table(DWI_list);
% DWI_list.ID = DWI_list.name;
% DWI_list.ID = erase(string(DWI_list.ID), '_hardiSD_ABS0.002_ang45_type2_pconn_d3.00_s0.00.mat');
% DWI_list = table2struct(DWI_list);

