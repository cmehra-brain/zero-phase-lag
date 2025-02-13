%% calculating edge-wise ICC
% This script produces a figure of test-retest reliability for each method
% in each freq band

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript: 
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------

%% clear variables and add paths as required
clear
clc

addpath('...') %add path to scripts

MainDirectory = '...EEG_networks_DK_1200_methodspaper/full_nets'; % define a main directory. The folder full_nets contains subfolders methods --> frequency bands 
netfiles = dir(fullfile(MainDirectory,'**/nets_ROIdata*'));

T = struct2table(netfiles); % convert the struct array to a table
sortedT = sortrows(T, {'name', 'folder'}); % sort by name followed by folder type
netfiles = table2struct(sortedT); % change it back to struct array if necessary
clear T
clear sortedT

path2output = '...';

%% Define participants

all_ages = {...}; %cell containing participant IDs included in analysis

%% STATE WHICH PARTICIPANTS YOU WANT TO SELECT FOR THIS ANALYSIS
selected_participants = all_ages;

%% make output matrices
coh_1_4 = nan*ones(68,68,2,length(netfiles));
coh_4_8 = nan*ones(68,68,2,length(netfiles));
coh_8_13 = nan*ones(68,68,2,length(netfiles));
coh_13_20 = nan*ones(68,68,2,length(netfiles));
coh_20_32 = nan*ones(68,68,2,length(netfiles));

img_coh_1_4 = nan*ones(68,68,2,length(netfiles));
img_coh_4_8 = nan*ones(68,68,2,length(netfiles));
img_coh_8_13 = nan*ones(68,68,2,length(netfiles));
img_coh_13_20 = nan*ones(68,68,2,length(netfiles));
img_coh_20_32 = nan*ones(68,68,2,length(netfiles));

PLV_1_4 = nan*ones(68,68,2,length(netfiles));
PLV_4_8 = nan*ones(68,68,2,length(netfiles));
PLV_8_13 = nan*ones(68,68,2,length(netfiles));
PLV_13_20 = nan*ones(68,68,2,length(netfiles));
PLV_20_32 = nan*ones(68,68,2,length(netfiles));

wPLI_1_4 = nan*ones(68,68,2,length(netfiles));
wPLI_4_8 = nan*ones(68,68,2,length(netfiles));
wPLI_8_13 = nan*ones(68,68,2,length(netfiles));
wPLI_13_20 = nan*ones(68,68,2,length(netfiles));
wPLI_20_32 = nan*ones(68,68,2,length(netfiles));

AEC_1_4 = nan*ones(68,68,2,length(netfiles));
AEC_4_8 = nan*ones(68,68,2,length(netfiles));
AEC_8_13 = nan*ones(68,68,2,length(netfiles));
AEC_13_20 = nan*ones(68,68,2,length(netfiles));
AEC_20_32 = nan*ones(68,68,2,length(netfiles));

orth_AEC_1_4 = nan*ones(68,68,2,length(netfiles));
orth_AEC_4_8 = nan*ones(68,68,2,length(netfiles));
orth_AEC_8_13 = nan*ones(68,68,2,length(netfiles));
orth_AEC_13_20 = nan*ones(68,68,2,length(netfiles));
orth_AEC_20_32 = nan*ones(68,68,2,length(netfiles));

%% run through for each adjacency matrix in each frequency band / FC method

for i = 1:length(netfiles) %bring up each set of networks one at a time
    disp(i)
    %figure out what the name is
    pre_ID = erase(string(netfiles(i).name), ".mat");
    ID = str2double(erase(string(pre_ID), "nets_ROIdata "));
    FC_freq = erase(string(netfiles(i).folder), "...EEG_networks_DK_1200_methodspaper/full_nets/");

    if any(strcmp(all_participants,num2str(ID))) %this equals one if Subject ID is in the list above

        %load the file
        cd(netfiles(i).folder);
        data = load(netfiles(i).name);
        fnames = fieldnames(data);
        temp_networks = data.(fnames{1});

        %get rid of the temporarily needed variables
        clear data
        clear fnames

        %only act if there are 2 networks or more - other participants are completely ignored
        if  size(temp_networks,3) > 1
            networks = temp_networks(:,:,1:2); %only keep 2, get rid of the rest

            % now, go through the whole network structure. See what FC
            % method and frequency each network was made in. Then, place
            % them each in their own matrix.

            %Note: These matrices will mainly consist of Nans
            if FC_freq == "COH/1_4Hz"
                coh_1_4(:,:,:,i) = networks;

            elseif FC_freq == "COH/4_8Hz"
                coh_4_8(:,:,:,i) = networks;

            elseif FC_freq == "COH/8_13Hz"
                coh_8_13(:,:,:,i) = networks;

            elseif FC_freq == "COH/13_20Hz"
                coh_13_20(:,:,:,i) = networks;

            elseif FC_freq == "COH/20_32Hz"
                coh_20_32(:,:,:,i) = networks;

            elseif FC_freq == "img_COH/1_4Hz"
                img_coh_1_4(:,:,:,i) = networks;

            elseif FC_freq == "img_COH/4_8Hz"
                img_coh_4_8(:,:,:,i) = networks;

            elseif FC_freq == "img_COH/8_13Hz"
                img_coh_8_13(:,:,:,i) = networks;

            elseif FC_freq == "img_COH/13_20Hz"
                img_coh_13_20(:,:,:,i) = networks;

            elseif FC_freq == "img_COH/20_32Hz"
                img_coh_20_32(:,:,:,i) = networks;

            elseif FC_freq == "PLV/1_4Hz"
                PLV_1_4(:,:,:,i) = networks;

            elseif FC_freq == "PLV/4_8Hz"
                PLV_4_8(:,:,:,i) = networks;

            elseif FC_freq == "PLV/8_13Hz"
                PLV_8_13(:,:,:,i) = networks;

            elseif FC_freq == "PLV/13_20Hz"
                PLV_13_20(:,:,:,i) = networks;

            elseif FC_freq == "PLV/20_32Hz"
                PLV_20_32(:,:,:,i) = networks;

            elseif FC_freq == "wPLI/1_4Hz"
                wPLI_1_4(:,:,:,i) = networks;

            elseif FC_freq == "wPLI/4_8Hz"
                wPLI_4_8(:,:,:,i) = networks;

            elseif FC_freq == "wPLI/8_13Hz"
                wPLI_8_13(:,:,:,i) = networks;

            elseif FC_freq == "wPLI/13_20Hz"
                wPLI_13_20(:,:,:,i) = networks;

            elseif FC_freq == "wPLI/20_32Hz"
                wPLI_20_32(:,:,:,i) = networks;

            elseif FC_freq == "AEC/1_4Hz"
                AEC_1_4(:,:,:,i) = networks;

            elseif FC_freq == "AEC/4_8Hz"
                AEC_4_8(:,:,:,i) = networks;

            elseif FC_freq == "AEC/8_13Hz"
                AEC_8_13(:,:,:,i) = networks;

            elseif FC_freq == "AEC/13_20Hz"
                AEC_13_20(:,:,:,i) = networks;

            elseif FC_freq == "AEC/20_32Hz"
                AEC_20_32(:,:,:,i) = networks;

            elseif FC_freq == "orth_AEC/1_4Hz"
                orth_AEC_1_4(:,:,:,i) = networks;

            elseif FC_freq == "orth_AEC/4_8Hz"
                orth_AEC_4_8(:,:,:,i) = networks;

            elseif FC_freq == "orth_AEC/8_13Hz"
                orth_AEC_8_13(:,:,:,i) = networks;

            elseif FC_freq == "orth_AEC/13_20Hz"
                orth_AEC_13_20(:,:,:,i) = networks;

            elseif FC_freq == "orth_AEC/20_32Hz"
                orth_AEC_20_32(:,:,:,i) = networks;

            end

        end

    end

end

%% get rid of nan 3D matrices

%now that you have 20 matrices (4 FC x 5 freq band), which are mostly nan,
%get rid of all the 3D matrices with only nans:

%create all the matrices in a structure
%then go through each of them and perform the same action

trt_mats = {coh_1_4,	coh_4_8,	coh_8_13,	coh_13_20,	coh_20_32,	img_coh_1_4,	img_coh_4_8,	img_coh_8_13,	img_coh_13_20,	img_coh_20_32,	PLV_1_4,	PLV_4_8,	PLV_8_13,	PLV_13_20,	PLV_20_32,	wPLI_1_4,	wPLI_4_8,	wPLI_8_13,	wPLI_13_20,	wPLI_20_32,	AEC_1_4,	AEC_4_8,	AEC_8_13,	AEC_13_20,	AEC_20_32,	orth_AEC_1_4,	orth_AEC_4_8,	orth_AEC_8_13,	orth_AEC_13_20,	orth_AEC_20_32};
names = {	'coh_1_4',	'coh_4_8',	'coh_8_13',	'coh_13_20',	'coh_20_32',	'img_coh_1_4',	'img_coh_4_8',	'img_coh_8_13',	'img_coh_13_20',	'img_coh_20_32',	'PLV_1_4',	'PLV_4_8',	'PLV_8_13',	'PLV_13_20',	'PLV_20_32',	'wPLI_1_4',	'wPLI_4_8',	'wPLI_8_13',	'wPLI_13_20',	'wPLI_20_32',	'AEC_1_4',	'AEC_4_8',	'AEC_8_13',	'AEC_13_20',	'AEC_20_32',	'orth_AEC_1_4',	'orth_AEC_4_8',	'orth_AEC_8_13',	'orth_AEC_13_20',	'orth_AEC_20_32'	};
EEG_FC_name = {'coh',	'coh',	'coh',	'coh',	'coh',	'img_coh',	'img_coh',	'img_coh',	'img_coh',	'img_coh',	'PLV',	'PLV',	'PLV',	'PLV',	'PLV',	'wPLI',	'wPLI',	'wPLI',	'wPLI',	'wPLI',	'AEC',	'AEC',	'AEC',	'AEC',	'AEC',	'orth_AEC',	'orth_AEC',	'orth_AEC',	'orth_AEC',	'orth_AEC'	};
freq = {'A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E','A', 'B','C','D','E'};
dummies = {nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
structure_trt = struct('matrices', trt_mats, 'method_freq',names, 'ICC_result' , dummies, 'single_ICC', dummies, 'freq_category', freq, 'method',EEG_FC_name);

for j = 1:length(structure_trt)

    trt_mat_1 = structure_trt(j).matrices;

    toRemove = false(size(trt_mat_1,4),1); %this creates an array (one column wide) of logical zeros. One for each 'page' in the 4th dimension. In this array, you will 'bookmark' pages that you want removed.
    for i = 1:size(trt_mat_1,4)
        if all(isnan(trt_mat_1(:,:,:,i))) %goes through each page of your 4th dimension - if there are all nans, then
            toRemove(i) = true; %applies the 'true' bookmark.
        end
    end

    trt_mat_1(:,:,:,toRemove) = []; %deletes bookmarked 'pages'

    structure_trt(j).matrices = trt_mat_1;
end

%% Calcualte ICC per edge
% you now have a structure with 30 rows. Each row contains a matrix of a
% given f with a given method
% you now want to create 30 ICC matrices.

for x = 1:length(structure_trt) %load each matrix one by one
    unitICC = zeros(68,68);
    ICC_mat_1 = structure_trt(x).matrices;

    for k = 1:size(ICC_mat_1,1) %68
        for m = 1:size(ICC_mat_1,1)
            pre_ICC = (squeeze(ICC_mat_1(k,m,:,:)))';
            unitICC(k,m) = ICC(pre_ICC, 'A-1'); %still some negatives
            %the intraclass correlation will be negative whenever the variability within groups exceeds the variability across groups. This means that scores in a group ``diverge'' relative to the noise present in the individuals.

            %this gets each edge through ALL PARTICIPANTS, THROUGH ALL
            %REPEATS (I.E. 2 REPEATS)
        end
    end

    structure_trt(x).ICC_result = unitICC;

    %then calculate the single median for that frequency range with that FC
    %method:

    %make each TOP HALF OF TRIANGLE, except the diagonal, into a column vector
    a = unitICC(triu(true(size(unitICC)) ,1));
    single_ICC = median(a);
    ICC_first_quantile = quantile(a,[0.25]);
    ICC_third_quantile = quantile(a,[0.75]);
    structure_trt(x).median_ICC = single_ICC;        
    structure_trt(x).IQR_twentyfive = ICC_first_quantile;
    structure_trt(x).IQR_seventyfive = ICC_third_quantile;

end

%%
%rearrange structure so that it's plotted properly
T = struct2table(structure_trt); % convert the struct array to a table
sortedT = sortrows(T, {'freq_category'}); % sort by name followed by folder type
structure_trt = table2struct(sortedT); % change it back to struct array if necessary
clear T
clear sortedT

%% reorder matrices based on proximity
anatomical_order = [31	13	11	18	17	19	26	3	27	16	23	21	28	30	7	24	20	12	4	10	33	29	1	14	8	6	5	15	32	34	25	2	22	9	65	47	45	52	51	53	60	37	61	50	57	55	62	64	41	58	54	46	38	44	67	63	35	48	42	40	39	49	66	68	59	36	56	43];
for k = 1:length(structure_trt)
    hi = structure_trt(k).ICC_result;
    hi = hi(anatomical_order,anatomical_order);
    structure_trt(k).ICC_result_reorder = hi;
end

%% produce all the adjacency matrices for TRT reliabiltiy

t = tiledlayout(5,6)
t.Padding = 'none';
t.TileSpacing = 'compact';

for z = 1:length(structure_trt)

    nexttile
    imagesc(structure_trt(z).ICC_result_reorder, [0 1])
    axis off
    axis square
    % cm = flip([linspace(1,0,20)',repmat(linspace(0,1,20)',1,2); repmat(linspace(1,0,20)',1,2),linspace(0,1,20)']);

end

%% focusing on only alpha
% alpha_trt = structure_trt;
% 
% T = struct2table(alpha_trt); % convert the struct array to a table
% T = T(13:18,:);
% reorder = {'A','D','B','E','C','F'}; %order so that the dirty methods are on top and the clean methods are below
% 
% for c = 1:length(reorder)
% T.freq_category(c) = reorder(c);
% end

%% calculating the variances of ICC values
% var_output = nan*zeros(6,2);
% for v = 1:length(alpha_trt)
% var_matrix = structure_trt(v).ICC_result;
% variances_trt = var_matrix(triu(true(size(var_matrix)) ,1));
% var_1 = var(variances_trt);
% var_output(v,2) = var_1;
% var_output(v,1) = convertCharsToStrings(alpha_trt(v).freq_category);
% end

%% formatting for one matrix with zoomed in ROIs

% %there are lines separating the lobes.
% 
% t = tiledlayout(2,3)
% t.Padding = 'none';
% t.TileSpacing = 'compact';
% 
% for z = 1:length(alpha_trt);
%     nexttile
%     imagesc(alpha_trt(z).ICC_result_reorder, [-0.6 1])
%     title(alpha_trt(z).method,'FontSize',30,'Interpreter','none') 
%     %colorbar('FontSize',18)
%     axis square
% 
%     % labels:
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
% end