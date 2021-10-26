% generate_paper_suppfigures.m
% 
% Matlab code to generate the supplemetary figures of the paper

%% LOAD CONNECTOME MATRICES AND REGION LOCATIONS/NAMES

load('data/connectome_human.mat')
load('data/connectome_chimp.mat')
load('data/connectome_macaque.mat')
load('data/connectome_marmoset.mat')
load('data/region_props_chimp.mat', 'region_centroids', 'region_names') 
load('data/RSN.mat', 'RSN_indices', 'RSN_names')
load('data/macaque_networks.mat', 'networks_indices', 'networks_names')

%% DEFINE FIGURE-RELATED PROPERTIES

Yeo_7net_colormap = [120, 18, 134;
                     70, 130, 180;
                     0, 118, 14;
                     196, 68, 250;
                     220, 248, 164;
                     230, 148, 34;
                     205, 62, 78;
                     200, 200, 200]/255;
Yeo_names = {'VIS', 'SM', 'DA', 'VA', 'LIM', 'FP', 'DM'};
macaque_network_names = {'FRONT', 'TEMP', 'PAR', 'OCC', 'LIM'};

cb = cbrewer('qual', 'Set1', 8, 'pchip');
color_brown = cb(7,:);
color_green = cb(3,:);
color_blue = cb(2,:);
color_purple = cb(4,:);
color_orange = cb(5,:);
cb = cbrewer('qual', 'Accent', 8, 'pchip');
color_brown_light = cb(7,:);
cb = cbrewer('qual', 'Dark2', 8, 'pchip');
color_green_light = cb(1,:);
color_brown_lighter = cb(7,:);
cb = cbrewer('qual', 'Set2', 8, 'pchip');
color_green_lighter = cb(1,:);

% artworks
human_female = imread('artworks/human_female_2.png');
human_male = imread('artworks/human_male_2.png');
chimpanzee = imread('artworks/chimpanzee.png');
macaque = imread('artworks/macaque.png');
marmoset = imread('artworks/marmoset.png');

fontsize_axis = 10;
fontsize_label = 12;

data_foldername = 'data/data_figures';

%% SUPPLEMENTARY FIGURE 1

% load all data relevant to Supplementary Figure 1
data_SuppFigure1 = load(sprintf('%s/SuppFigure1.mat', data_foldername));

titles = {'human', 'macaque'};

cmap = [color_brown; color_blue];
          
N = 114;
slice = 'axial';   
node_interest = 86;
frac = 0.05;
markersize = 100;
node_interests_colors = brighten(lines(3),0.3);
node_interests_colors = node_interests_colors(2,:);
node_cmap = [1, 1, 1; node_interests_colors];
node_cval = ones(N,1);
node_cval(node_interest) = 2;
nodeLocations = region_centroids;
nodeSizes = (markersize/10)*ones(N,1);
nodeSizes(node_interest) = markersize;
edge_width = 1;
image_width = 0.08;

fig = figure('Position', [200 100 700 500]);
% =========================================================================
% A left: brain network
% =========================================================================
ax1 = axes('Position', [0.03 0.7 0.25 0.25]);
 
W = connectome_human;
[edge_X, edge_Y, edge_Z] = adjacency_plot_und(threshold_proportional(W, frac), nodeLocations);  % get all the edges
[nodeLocations, edges] = utils.extract_scatterBrain_locs_edges(nodeLocations, edge_X, edge_Y, edge_Z, slice);
 
hold on;
if strcmpi(slice, 'axial')
    scatter3(ax1, nodeLocations(:,1), nodeLocations(:,2), max(nodeLocations(:,3))*ones(N,1), ...
             nodeSizes, node_cval, 'filled', 'markeredgecolor', 'k');
else
    scatter3(ax1, nodeLocations(:,1), nodeLocations(:,2), nodeLocations(:,3), ...
         nodeSizes, node_cval, 'filled', 'markeredgecolor', 'k');
end
plot3(edges(:,1), edges(:,2), edges(:,3), 'color', 'k', 'linewidth', edge_width);
hold off;
text(nodeLocations(node_interest,1)+8, nodeLocations(node_interest,2), max(nodeLocations(:,3)), 'region $i$', ...
    'horizontalalignment', 'left', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(-40, -46, 'LH', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, ...
    'fontweight', 'b')
text(40, -46, 'RH', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, ...
    'fontweight', 'b')
colormap(ax1, node_cmap)
view(ax1, 2)
axis(ax1, 'equal')
set(ax1, 'visible', 'off')

% =========================================================================
% A middle top: simulated neural activity
% =========================================================================
ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)+0.1 ax1.Position(2)+ax1.Position(4)*0.65 0.2 0.05]);
T = data_SuppFigure1.timeseries_neural.time;
plot(T, data_SuppFigure1.timeseries_neural.S(node_interest,:), 'color', node_interests_colors(1,:), 'linewidth', 1.5)
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [min(T) max(T)], 'xtick', [], 'ytick', [], 'ycolor', 'none', 'box', 'off');
axis off
title('simulated neural activity', 'fontsize', fontsize_axis)

% =========================================================================
% A middle bottom: simulated fMRI signal
% ========================================================================= 
ax3 = axes('Position', [ax2.Position(1) ax1.Position(2) ax2.Position(3) ax2.Position(4)]);
T = data_SuppFigure1.timeseries_fMRI.time;
plot(T, data_SuppFigure1.timeseries_fMRI.BOLD(node_interest,:), 'color', node_interests_colors(1,:), 'linewidth', 1.5)
set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [min(T) max(T)], 'xtick', [], 'ytick', [], 'ycolor', 'none', 'box', 'off');
axis off
title('simulated fMRI signal', 'fontsize', fontsize_axis)

% =========================================================================
% A right: simulated functional connectivity
% =========================================================================  
ax4 = axes('Position', [ax2.Position(1)+ax2.Position(3)+0.12 ax3.Position(2) 0.3 0.25]);
FC = corr(data_SuppFigure1.timeseries_fMRI.BOLD');
imagesc(FC)
caxis([-1 1])
colormap(ax4, utils.bluewhitered)
cbar = colorbar;
set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xtick', [], 'ytick', []);
set(cbar, 'fontsize', fontsize_axis-2)
xlabel('region', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('region', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel(cbar, '$FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
axis square
title('simulated functional connectivity', 'fontsize', fontsize_axis)

% =========================================================================
% B left: simulated vs empirical human functional connectivity
% =========================================================================  
ax5 = axes('Position', [ax1.Position(1)+0.02 0.1 0.48 0.48]);
FC_emp = data_SuppFigure1.FC_human.emp;
FC_sim = data_SuppFigure1.FC_human.sim;
tri_ind = find(triu(ones(size(FC_emp,1)),1));
data_to_plot_x = FC_emp(tri_ind);
data_to_plot_y = FC_sim(tri_ind);
[rho_temp,pval_temp] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
limits_min = min([data_to_plot_x; data_to_plot_y]);
 
hold on;
plot(data_to_plot_x, data_to_plot_y, '.', 'markersize', 8, 'color', cmap(1,:))
plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'k-', 'linewidth', 2);
hold off;
set(ax5, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [limits_min 1], 'ylim', [limits_min 1], 'xtick', -1:0.2:1, 'ytick', -1:0.2:1)
text(max(get(gca, 'xlim')), min(get(gca, 'ylim'))*0.9+0.1, sprintf('r = %.2g', rho_temp), ...
        'fontsize', 10, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
text(max(get(gca, 'xlim')), min(get(gca, 'ylim'))*0.9, sprintf('p = %.2g', pval_temp), ...
    'fontsize', 10, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
xlabel('empirical $FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('simulated $FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
axis square
title(titles{1}, 'fontsize', fontsize_label)
 
image_to_plot = human_female;
bw = image_to_plot>0;
ax5_image = axes('Position', [ax5.Position(1)+ax5.Position(3)*0.15 ax5.Position(2)+ax5.Position(4)*0.9 image_width image_width]);
imshow(bw);
colormap(ax5_image, flipud(gray))

% =========================================================================
% B right: simulated vs empirical macaque functional connectivity
% =========================================================================
ax6 = axes('Position', [ax5.Position(1)+ax5.Position(3) ax5.Position(2) ax5.Position(3) ax5.Position(4)]);
FC_emp = data_SuppFigure1.FC_macaque.emp;
FC_sim = data_SuppFigure1.FC_macaque.sim;
tri_ind = find(triu(ones(size(FC_emp,1)),1));
data_to_plot_x = FC_emp(tri_ind);
data_to_plot_y = FC_sim(tri_ind);
[rho_temp,pval_temp] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
limits_min = min([data_to_plot_x; data_to_plot_y]);
 
hold on;
plot(data_to_plot_x, data_to_plot_y, '.', 'markersize', 8, 'color', cmap(2,:))
plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'k-', 'linewidth', 2);
hold off;
set(ax6, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [limits_min 1], 'ylim', [limits_min 1], 'xtick', -1:0.2:1, 'ytick', -1:0.2:1)
text(max(get(gca, 'xlim')), min(get(gca, 'ylim'))*0.9+0.1, sprintf('r = %.2g', rho_temp), ...
        'fontsize', 10, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
text(max(get(gca, 'xlim')), min(get(gca, 'ylim'))*0.9, sprintf('p = %.2g', pval_temp), ...
    'fontsize', 10, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
xlabel('empirical $FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('simulated $FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
axis square
title(titles{2}, 'fontsize', fontsize_label)
 
image_to_plot = macaque;
bw = image_to_plot>0;
ax6_image = axes('Position', [ax6.Position(1)+ax6.Position(3)*0.15 ax6.Position(2)+ax6.Position(4)*0.9 image_width image_width]);
imshow(bw);
colormap(ax6_image, flipud(gray))
 
%%% arrows and text
annotation('arrow', [ax2.Position(1)-0.08 ax2.Position(1)-0.02], [ax2.Position(2)+ax2.Position(4)*0.5 ax2.Position(2)+ax2.Position(4)*0.5])
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*0.5 ax2.Position(1)+ax2.Position(3)*0.5], [ax2.Position(2) ax3.Position(2)+ax3.Position(4)*1.85])
annotation('arrow', [ax3.Position(1)+ax3.Position(3)+0.02 ax3.Position(1)+ax3.Position(3)+0.08], [ax3.Position(2)+ax3.Position(4)*0.5 ax3.Position(2)+ax3.Position(4)*0.5])
annotation('textbox', [ax2.Position(1)+ax2.Position(3)*0.75, ax3.Position(2)+ax3.Position(4)*1.5, 0.05, 0.1], 'string', 'hemodynamic model', 'edgecolor', 'none', ...
            'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
 
%%% panel letters
annotation(fig, 'textbox', [0.05, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.05, 0.65, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
  
%% SUPPLEMENTARY FIGURE 2

% load all data relevant to Supplementary Figure 2
data_SuppFigure2 = load(sprintf('%s/SuppFigure2.mat', data_foldername));

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];
S_transition = 0.3;
image_width = 0.05;
 
fig = figure('Position', [400 200 600 800]);
% =========================================================================
% A: tuning curve
% =========================================================================
for type_ind = 1:length(types)
    type = types{type_ind};
    
    stats_temp = data_SuppFigure2.tuningcurve.stats{type_ind};
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    gain_cmap = repmat(cmap(type_ind,:),length(stats_temp.max),1);
    brighten_vec = linspace(0.8,0,length(stats_temp.max));
    
    [~,xval_transition_ind] = sort(stats_temp.xval_transition, 'ascend');
    
    %%% A: dynamic range sigmoids
    subplot(3,2,type_ind)
    ax1 = gca;
    hold on;
    for ii=1:length(stats_temp.max)
        plot(data_SuppFigure2.tuningcurve.w, data_SuppFigure2.tuningcurve.S{type_ind}(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
    end
    plot([min(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition, S_transition], 'k-', 'linewidth', 2)
    plot([min(stats_temp.xval_transition), min(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    plot([max(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    hold off;
    text(max(stats_temp.xval_transition)*1.04, S_transition, ...
        ['width = ', sprintf('%.2f', stats_temp.dx_transition)])
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel('global recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax1_image= axes('Position', [ax1.Position(1)+0.01 ax1.Position(2)+ax1.Position(4)*0.75 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
end

% =========================================================================
% B: dynamic range std
% =========================================================================
DR_std = {data_SuppFigure2.subject_stats.DR_std{1}; data_SuppFigure2.subject_stats.DR_std{2}};
grp = cell2mat(arrayfun(@(i){i*ones(numel(DR_std{i}),1)},(1:numel(DR_std))')); 
data_violin = utils.padconcatenation(data_SuppFigure2.subject_stats.DR_std{1}, data_SuppFigure2.subject_stats.DR_std{2}, 2);
  
subplot(3,2,3)
set(gca, 'Position', get(gca, 'Position')+[0.21 0 0 0 ])
ax2 = gca;
violins = utils.violinplot(data_violin, {});
for type_ind = 1:length(types)
    violins(type_ind).MedianPlot.SizeData = 30;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 15;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.02;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
end
 
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.01, 0.01], 'xticklabel', titles)
set(gca, 'ylim', get(gca, 'ylim').*[0 1])
ylabel({'dynamic range $\sigma$'; 'per subject'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
[~,pval] = ttest2(DR_std{1}, DR_std{2});
text(0.6, max(get(gca, 'ylim'))*0.9, sprintf('p = %.2g', pval), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'middle', 'horizontalalignment', 'left');

% =========================================================================
% C: dynamic range std versus brain volume
% =========================================================================
for type_ind=1:length(types)
    type = types{type_ind};
    
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    
    subplot(3,2,4+type_ind)
    ax3 = gca;
    data_to_plot_x = data_SuppFigure2.subject_stats.volume{type_ind};
    data_to_plot_y = data_SuppFigure2.subject_stats.DR_std{type_ind};
 
    hold on;
    plot(data_to_plot_x, data_to_plot_y, '.', 'color', cmap(type_ind,:), 'markersize', 15)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
            'k-', 'linewidth', 2);
    hold off;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', get(ax2, 'ylim'));
    xlabel('brain volume (mm$^3$)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'dynamic range $\sigma$'; 'per subject'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    text(max(get(gca,'xlim')), max(get(gca, 'ylim'))*0.9, sprintf('r = %.2g', rho), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
    text(max(get(gca,'xlim')), max(get(gca, 'ylim'))*0.8, sprintf('p = %.2g', pval), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
    
    ax3_image= axes('Position', [ax3.Position(1)+0.01 ax3.Position(2)+ax3.Position(4)*0.75 image_width image_width]);
    imshow(bw);
    colormap(ax3_image, flipud(gray))
end

%%% panel letters
annotation(fig, 'textbox', [0.05, 0.96, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.26, 0.64, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.05, 0.36, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
 
%% SUPPLEMENTARY FIGURE 3

% load all data relevant to Supplementary Figure 3
data_SuppFigure3 = load(sprintf('%s/SuppFigure3.mat', data_foldername));

types = {'human', 'chimp', 'human_equidensity'};
titles = {'human', 'chimpanzee', 'human (equal density)'};

cmap = [color_brown; color_green; color_brown_lighter];
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';
image_width = 0.07;

fig = figure('Position', [200 200 800 500]);
% =========================================================================
% A: tuning curve
% =========================================================================
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    elseif type_ind==3
        image_to_plot = human_female;
    end
    bw = image_to_plot>0;
    gain_cmap = repmat(cmap(type_ind,:),length(data_SuppFigure3.tuningcurve.stats{type_ind}.max),1);
    brighten_vec = linspace(0.8,0,length(data_SuppFigure3.tuningcurve.stats{type_ind}.max));

    [~,xval_transition_ind] = sort(data_SuppFigure3.tuningcurve.stats{type_ind}.xval_transition, 'ascend');
    
    subplot(2,length(types),type_ind)
    ax1 = gca;
    hold on;
    for ii=1:length(data_SuppFigure3.tuningcurve.stats{type_ind}.max)
        plot(data_SuppFigure3.tuningcurve.w, data_SuppFigure3.tuningcurve.S{type_ind}(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
    end
    plot([min(data_SuppFigure3.tuningcurve.stats{type_ind}.xval_transition), max(data_SuppFigure3.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition, S_transition], 'k-', 'linewidth', 2)
    plot([min(data_SuppFigure3.tuningcurve.stats{type_ind}.xval_transition), min(data_SuppFigure3.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    plot([max(data_SuppFigure3.tuningcurve.stats{type_ind}.xval_transition), max(data_SuppFigure3.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    hold off;
    text(max(data_SuppFigure3.tuningcurve.stats{type_ind}.xval_transition)*1.05, S_transition, ...
        ['width = ', sprintf('%.2f', data_SuppFigure3.tuningcurve.stats{type_ind}.dx_transition)])
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel('global recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax1_image= axes('Position', [ax1.Position(1)-0.01 ax1.Position(2)+ax1.Position(4)*0.75 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
end

% =========================================================================
% B: dynamic range distribution
% =========================================================================
subplot(2,length(types),(length(types)+1):length(types)*2)
ax2 = gca;

data_violin = [];
for type_ind = 1:length(types)
    data_violin = utils.padconcatenation(data_violin, data_SuppFigure3.tuningcurve.stats{type_ind}.(stats_to_plot)-mean(data_SuppFigure3.tuningcurve.stats{type_ind}.(stats_to_plot)), 2);
end
    
violins = utils.violinplot(data_violin, {});
for type_ind = 1:length(types)
    violins(type_ind).MedianPlot.SizeData = 50;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 10;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.03;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
    text(type_ind, 1.3, titles{type_ind}, 'color', cmap(type_ind,:), ...
    'fontsize', fontsize_axis, 'fontweight', 'b', 'verticalalignment', 'bottom')
    text(type_ind, 1.3, ['(\sigma = ', sprintf('%.2f)', std(data_SuppFigure3.tuningcurve.stats{type_ind}.(stats_to_plot)-mean(data_SuppFigure3.tuningcurve.stats{type_ind}.(stats_to_plot))))], 'color', cmap(type_ind,:), ...
        'fontsize', fontsize_axis, 'fontweight', 'b', 'horizontalalignment', 'left', 'verticalalignment', 'top')
end
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0.5, length(types)+0.5], 'ylim', [-1.8 1.8], 'xtick', []);
ylabel([stats_to_plot_name, ' (mean-subtracted)'], 'fontsize', fontsize_label, 'interpreter', 'latex')
view(90, 90);
box off
ax2.XAxis.Visible = 'off';

for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    elseif type_ind==3
        image_to_plot = human_female;
    end
    bw = image_to_plot>0;
    
    ax2_image= axes('Position', [ax2.Position(1)+ax2.Position(3)*0.78 ax2.Position(2)+ax2.Position(4)*0.73-0.115*(type_ind-1) image_width image_width]);
    imshow(bw);
    colormap(ax2_image, flipud(gray))
end

%%% panel letters
annotation(fig, 'textbox', [0.07, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.07, 0.49, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% SUPPLEMENTARY FIGURE 4

% load all data relevant to Supplementary Figure 4
data_SuppFigure4 = load(sprintf('%s/SuppFigure4.mat', data_foldername));

types = {'human', 'chimp', 'human_noisy', 'chimp_noisy'};
titles = {'human', 'chimpanzee', 'human (rescaled)', 'chimpanzee (rescaled)'};

cmap = [color_brown; color_green; color_brown_lighter; color_green_lighter];
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';
image_width = 0.06;

fig = figure('Position', [200 200 800 500]);
% =========================================================================
% A: tuning curve
% =========================================================================
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    elseif type_ind==3
        image_to_plot = human_female;
    elseif type_ind==4
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    gain_cmap = repmat(cmap(type_ind,:),length(data_SuppFigure4.tuningcurve.stats{type_ind}.max),1);
    brighten_vec = linspace(0.8,0,length(data_SuppFigure4.tuningcurve.stats{type_ind}.max));

    [~,xval_transition_ind] = sort(data_SuppFigure4.tuningcurve.stats{type_ind}.xval_transition, 'ascend');
    
    subplot(2,length(types),type_ind)
    ax1 = gca;
    hold on;
    for ii=1:length(data_SuppFigure4.tuningcurve.stats{type_ind}.max)
        plot(data_SuppFigure4.tuningcurve.w, data_SuppFigure4.tuningcurve.S{type_ind}(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
    end
    plot([min(data_SuppFigure4.tuningcurve.stats{type_ind}.xval_transition), max(data_SuppFigure4.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition, S_transition], 'k-', 'linewidth', 2)
    plot([min(data_SuppFigure4.tuningcurve.stats{type_ind}.xval_transition), min(data_SuppFigure4.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    plot([max(data_SuppFigure4.tuningcurve.stats{type_ind}.xval_transition), max(data_SuppFigure4.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    hold off;
    text(max(data_SuppFigure4.tuningcurve.stats{type_ind}.xval_transition)*1.05, S_transition, ...
        ['width = ', sprintf('%.2f', data_SuppFigure4.tuningcurve.stats{type_ind}.dx_transition)])
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel({'global recurrent'; 'strength, $w$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax1_image= axes('Position', [ax1.Position(1)-0.01 ax1.Position(2)+ax1.Position(4)*0.78 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
end

% =========================================================================
% B: dynamic range distribution
% =========================================================================
subplot(2,length(types),(length(types)+1):length(types)*2)
ax2 = gca;

data_violin = [];
for type_ind = 1:length(types)
    data_violin = utils.padconcatenation(data_violin, data_SuppFigure4.tuningcurve.stats{type_ind}.(stats_to_plot)-mean(data_SuppFigure4.tuningcurve.stats{type_ind}.(stats_to_plot)), 2);
end
    
violins = utils.violinplot(data_violin, {});
for type_ind = 1:length(types)
    violins(type_ind).MedianPlot.SizeData = 50;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 10;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.03;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
    text(type_ind, 1.3, titles{type_ind}, 'color', cmap(type_ind,:), ...
    'fontsize', fontsize_axis, 'fontweight', 'b', 'verticalalignment', 'bottom')
    text(type_ind, 1.3, ['(\sigma = ', sprintf('%.2f)', std(data_SuppFigure4.tuningcurve.stats{type_ind}.(stats_to_plot)-mean(data_SuppFigure4.tuningcurve.stats{type_ind}.(stats_to_plot))))], 'color', cmap(type_ind,:), ...
        'fontsize', fontsize_axis, 'fontweight', 'b', 'horizontalalignment', 'left', 'verticalalignment', 'top')
end
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0.5, length(types)+0.5], 'ylim', [-1.8 1.8], 'xtick', []);
ylabel([stats_to_plot_name, ' (mean-subtracted)'], 'fontsize', fontsize_label, 'interpreter', 'latex')
view(90, 90);
box off
ax2.XAxis.Visible = 'off';

for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    elseif type_ind==3
        image_to_plot = human_female;
    elseif type_ind==4
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    
    ax2_image= axes('Position', [ax2.Position(1)+ax2.Position(3)*0.79 ax2.Position(2)+ax2.Position(4)*0.8-0.088*(type_ind-1) image_width image_width]);
    imshow(bw);
    colormap(ax2_image, flipud(gray))
end

%%% panel letters
annotation(fig, 'textbox', [0.07, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.07, 0.49, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% SUPPLEMENTARY FIGURE 5

% load all data relevant to Supplementary Figure 5
data_SuppFigure5 = load(sprintf('%s/SuppFigure5.mat', data_foldername)); 

types = {'human', 'chimp', 'human_equisample'};
titles = {'human', 'chimpanzee', 'human (equal sample)'};

cmap = [color_brown; color_green; color_brown_lighter];
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';
image_width = 0.05;
 
fig = figure('Position', [200 100 800 800]);
% =========================================================================
% A: tuning curve
% =========================================================================
for type_ind = 1:length(types)
    S_temp = data_SuppFigure5.tuningcurve.S{type_ind};
    if type_ind==1 || type_ind==2
        stats_temp = data_SuppFigure5.tuningcurve.stats{type_ind};
    else
        stats_temp = data_SuppFigure5.tuningcurve.stats{type_ind}{1};
    end
 
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    elseif type_ind==3
        image_to_plot = human_female;
    end
    bw = image_to_plot>0;
    gain_cmap = repmat(cmap(type_ind,:),length(stats_temp.max),1);
    brighten_vec = linspace(0.8,0,length(stats_temp.max));
 
    [~,xval_transition_ind] = sort(stats_temp.xval_transition, 'ascend');
    
    %%% A: dynamic range sigmoids
    subplot(3,length(types),type_ind)
    ax1 = gca;
    hold on;
    for ii=1:length(stats_temp.max)
        plot(data_SuppFigure5.tuningcurve.w, S_temp(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
    end
    plot([min(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition, S_transition], 'k-', 'linewidth', 2)
    plot([min(stats_temp.xval_transition), min(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    plot([max(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    hold off;
    text(max(stats_temp.xval_transition)*1.05, S_transition, ...
        ['width = ', sprintf('%.2f', stats_temp.dx_transition)])
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel('global recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax1_image= axes('Position', [ax1.Position(1) ax1.Position(2)+ax1.Position(4)*0.75 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
end
 
% =========================================================================
% B: dynamic range distribution
% =========================================================================
subplot(3,length(types),(length(types)+1):length(types)*2)
ax2 = gca;
data_violin = [];
for type_ind = 1:length(types)
    S_temp = data_SuppFigure5.tuningcurve.S{type_ind};
    if type_ind==1 || type_ind==2
        stats_temp = data_SuppFigure5.tuningcurve.stats{type_ind};
    else
        stats_temp = data_SuppFigure5.tuningcurve.stats{type_ind}{1};
    end
    
    data_violin = utils.padconcatenation(data_violin, stats_temp.(stats_to_plot)-mean(stats_temp.(stats_to_plot)), 2);
end
    
violins = utils.violinplot(data_violin, {});
for type_ind = 1:length(types)
    S_temp = data_SuppFigure5.tuningcurve.S{type_ind};
    if type_ind==1 || type_ind==2
        stats_temp = data_SuppFigure5.tuningcurve.stats{type_ind};
    else
        stats_temp = data_SuppFigure5.tuningcurve.stats{type_ind}{1};
    end
    
    violins(type_ind).MedianPlot.SizeData = 50;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 10;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.03;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
    text(type_ind, 1.3, titles{type_ind}, 'color', cmap(type_ind,:), ...
    'fontsize', fontsize_axis, 'fontweight', 'b', 'verticalalignment', 'bottom')
    text(type_ind, 1.3, ['(\sigma = ', sprintf('%.2f)', std(stats_temp.(stats_to_plot)-mean(stats_temp.(stats_to_plot))))], 'color', cmap(type_ind,:), ...
        'fontsize', fontsize_axis, 'fontweight', 'b', 'horizontalalignment', 'left', 'verticalalignment', 'top')
end
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0.5, length(types)+0.5], 'ylim', [-1.8 1.8], 'xtick', []);
ylabel([stats_to_plot_name, ' (mean-subtracted)'], 'fontsize', fontsize_label, 'interpreter', 'latex')
view(90, 90);
box off
ax2.XAxis.Visible = 'off';
 
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    elseif type_ind==3
        image_to_plot = human_female;
    elseif type_ind==4
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    
    ax2_image= axes('Position', [ax2.Position(1)+ax2.Position(3)*0.79 ax2.Position(2)+ax2.Position(4)*0.73-0.075*(type_ind-1) image_width image_width]);
    imshow(bw);
    colormap(ax2_image, flipud(gray))
end
 
% =========================================================================
% C: dynamic range std
% =========================================================================
subplot(3,length(types),2*length(types)+(1:length(types)))
stats_std_2 = zeros(1,length(data_SuppFigure5.tuningcurve.stats{3}));
for ii=1:length(data_SuppFigure5.tuningcurve.stats{3})
    stats_std_2(ii) = std(data_SuppFigure5.tuningcurve.stats{3}{ii}.(stats_to_plot));
end
hold on;
plot(1:ii, std(data_SuppFigure5.tuningcurve.stats{1}.(stats_to_plot))*ones(1,length(stats_std_2)), '-', 'color', cmap(1,:), 'linewidth', 2);
plot(1:ii, std(data_SuppFigure5.tuningcurve.stats{2}.(stats_to_plot))*ones(1,length(stats_std_2)), '-', 'color', cmap(2,:), 'linewidth', 2);
plot(1:ii, stats_std_2, 'k.', 'color', cmap(3,:), 'markersize', 20)
hold off;
legend(titles, 'location', 'best')
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02]);
xlabel('trial', 'fontsize', fontsize_label)
ylabel('dynamic range \sigma', 'fontsize', fontsize_label)
 
%%% panel letters
annotation(fig, 'textbox', [0.07, 0.97, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.07, 0.65, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.07, 0.36, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% SUPPLEMENTARY FIGURE 6

% load all data relevant to Supplementary Figure 6
data_SuppFigure6 = load(sprintf('%s/SuppFigure6.mat', data_foldername));

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';
image_width = 0.05;

fig = figure('Position', [200 100 800 800]);
% =========================================================================
% A: time delay distribution
% =========================================================================
subplot(3,length(types),1)
ax1 = gca;
set(ax1, 'Position', ax1.Position+[0.22 0 0 0])
violins = utils.violinplot(data_SuppFigure6.time_delay, {});
for type_ind = 1:length(types)
    violins(type_ind).MedianPlot.SizeData = 30;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 1;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.01;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
end
set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xticklabel', titles, 'ylim', [0 1.3*max(data_SuppFigure6.time_delay(:))])
ylabel('time delay, $t^d$ (ms)', 'fontsize', fontsize_label, 'interpreter', 'latex')

for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;

    ax1_image = axes('Position', [ax1.Position(1)+0.059+(type_ind-1)*0.17 ax1.Position(2)+ax1.Position(4)*0.8 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
end

% =========================================================================
% B: tuning curve
% =========================================================================
for type_ind = 1:length(types)
    S_temp = data_SuppFigure6.tuningcurve.S{type_ind};
    stats_temp = data_SuppFigure6.tuningcurve.stats{type_ind};
 
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    gain_cmap = repmat(cmap(type_ind,:),length(stats_temp.max),1);
    brighten_vec = linspace(0.8,0,length(stats_temp.max));
 
    [~,xval_transition_ind] = sort(stats_temp.xval_transition, 'ascend');
    
    %%% B: dynamic range sigmoids
    subplot(3,length(types),length(types)+type_ind)
    ax2 = gca;
    hold on;
    for ii=1:length(stats_temp.max)
        plot(data_SuppFigure6.tuningcurve.w, S_temp(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
    end
    plot([min(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition, S_transition], 'k-', 'linewidth', 2)
    plot([min(stats_temp.xval_transition), min(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    plot([max(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    hold off;
    text(max(stats_temp.xval_transition)*1.05, S_transition, ...
        ['width = ', sprintf('%.2f', stats_temp.dx_transition)])
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel('global recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    
    ax2_image= axes('Position', [ax2.Position(1) ax2.Position(2)+ax2.Position(4)*0.75 image_width image_width]);
    imshow(bw);
    colormap(ax2_image, flipud(gray))
end
 
% =========================================================================
% C: dynamic range distribution
% =========================================================================
subplot(3,length(types),(2*length(types)+1):length(types)*3)
ax3 = gca;
data_violin = [];
for type_ind = 1:length(types)
    S_temp = data_SuppFigure6.tuningcurve.S{type_ind};
    stats_temp = data_SuppFigure6.tuningcurve.stats{type_ind};
        
    data_violin = utils.padconcatenation(data_violin, stats_temp.(stats_to_plot)-mean(stats_temp.(stats_to_plot)), 2);
end
    
violins = utils.violinplot(data_violin, {});
for type_ind = 1:length(types)
    S_temp = data_SuppFigure6.tuningcurve.S{type_ind};
    stats_temp = data_SuppFigure6.tuningcurve.stats{type_ind};
    
    violins(type_ind).MedianPlot.SizeData = 50;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 10;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.02;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
    text(type_ind, 1.3, titles{type_ind}, 'color', cmap(type_ind,:), ...
    'fontsize', fontsize_axis, 'fontweight', 'b', 'verticalalignment', 'bottom')
    text(type_ind, 1.3, ['(\sigma = ', sprintf('%.2f)', std(stats_temp.(stats_to_plot)-mean(stats_temp.(stats_to_plot))))], 'color', cmap(type_ind,:), ...
        'fontsize', fontsize_axis, 'fontweight', 'b', 'horizontalalignment', 'left', 'verticalalignment', 'top')
end
set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0.5, length(types)+0.5], 'ylim', [-1.8 1.8], 'xtick', []);
ylabel([stats_to_plot_name, ' (mean-subtracted)'], 'fontsize', fontsize_label, 'interpreter', 'latex')
view(90, 90);
box off
ax3.XAxis.Visible = 'off';
 
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    
    ax3_image= axes('Position', [ax3.Position(1)+ax3.Position(3)*0.79 ax3.Position(2)+ax3.Position(4)*0.65-0.11*(type_ind-1) image_width image_width]);
    imshow(bw);
    colormap(ax3_image, flipud(gray))
end
 
%%% panel letters
annotation(fig, 'textbox', [0.29, 0.95, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.07, 0.66, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.07, 0.35, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% SUPPLEMENTARY FIGURE 7

% load all data relevant to Supplementary Figure 7
data_SuppFigure7 = load(sprintf('%s/SuppFigure7.mat', data_foldername)); 

types = {'human', 'HCP_100_cortex_thres'};
titles = {'human', 'human (HCP)'};

cmap = [color_brown; color_brown_lighter];
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';
 
fig = figure('Position', [200 100 500 400]);
% =========================================================================
% A: tuning curve
% =========================================================================
for type_ind = 1:length(types)
    S_temp = data_SuppFigure7.tuningcurve.S{type_ind};
    stats_temp = data_SuppFigure7.tuningcurve.stats{type_ind};
    
    gain_cmap = repmat(cmap(type_ind,:),length(stats_temp.max),1);
    brighten_vec = linspace(0.8,0,length(stats_temp.max));
 
    [~,xval_transition_ind] = sort(stats_temp.xval_transition, 'ascend');
    
    subplot(2,length(types),type_ind)
    ax1 = gca;
    hold on;
    for ii=1:length(stats_temp.max)
        plot(data_SuppFigure7.tuningcurve.w, S_temp(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
    end
    plot([min(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition, S_transition], 'k-', 'linewidth', 2)
    plot([min(stats_temp.xval_transition), min(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    plot([max(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    hold off;
    text(max(stats_temp.xval_transition)*1.05, S_transition, ...
        ['width = ', sprintf('%.2f', stats_temp.dx_transition)])
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel('global recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
end
 
% =========================================================================
% B: dynamic range distribution
% =========================================================================
subplot(2,length(types),(length(types)+1):length(types)*2)
ax2 = gca;
data_violin = [];
for type_ind = 1:length(types)
    data_violin = utils.padconcatenation(data_violin, data_SuppFigure7.tuningcurve.stats{type_ind}.(stats_to_plot)-mean(data_SuppFigure7.tuningcurve.stats{type_ind}.(stats_to_plot)), 2);
end
    
violins = utils.violinplot(data_violin, {});
for type_ind = 1:length(types)
    violins(type_ind).MedianPlot.SizeData = 50;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 15;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.03;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
    text(type_ind, 1.3, titles{type_ind}, 'color', cmap(type_ind,:), ...
    'fontsize', fontsize_axis, 'fontweight', 'b', 'verticalalignment', 'bottom')
    text(type_ind, 1.3, ['(\sigma = ', sprintf('%.2f)', std(data_SuppFigure7.tuningcurve.stats{type_ind}.(stats_to_plot)-mean(data_SuppFigure7.tuningcurve.stats{type_ind}.(stats_to_plot))))], 'color', cmap(type_ind,:), ...
        'fontsize', fontsize_axis, 'fontweight', 'b', 'horizontalalignment', 'left', 'verticalalignment', 'top')
end
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0.5, length(types)+0.5], 'ylim', [-1.8 1.8], 'xtick', []);
ylabel([stats_to_plot_name, ' (mean-subtracted)'], 'fontsize', fontsize_label, 'interpreter', 'latex')
view(90, 90);
box off
ax2.XAxis.Visible = 'off';
 
%%% panel letters
annotation(fig, 'textbox', [0.07, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.07, 0.49, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% SUPPLEMENTARY FIGURE 8

% load all data relevant to Supplementary Figure 8
data_SuppFigure8 = load(sprintf('%s/SuppFigure8.mat', data_foldername));

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];
S_transition = 0.15;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';
image_width = 0.05;

N = 114;
slice = 'axial';   
node_interests = [86,84];
frac = 0.05;
markersize = 100;
node_interests_colors = brighten(lines(6),0.3);
node_interests_colors = node_interests_colors([4,6],:);
node_cmap = [1, 1, 1; node_interests_colors];
node_cval = ones(N,1);
node_cval(node_interests(1)) = 2;
node_cval(node_interests(2)) = 3;
nodeSizes = (markersize/10)*ones(N,1);
nodeSizes(node_interests) = markersize;
nodeLocations = region_centroids;
positions = [0.1, 0.53, 0.42, 0.42; ...
             0.1, 0.05, 0.4, 0.4];
edge_alpha = 0.8;
edge_width = 1;
edge_colors = [cmap(1,:), edge_alpha;
               cmap(2,:), edge_alpha];
 
fig = figure('Position', [200 100 650 700]);
 
% =========================================================================
% A: brain network
% =========================================================================
ax1 = axes('Position', [0.05 0.7 0.25 0.25]);
 
W = connectome_human;
[edge_X, edge_Y, edge_Z] = adjacency_plot_und(threshold_proportional(W, frac), nodeLocations);  % get all the edges
[nodeLocations, edges] = utils.extract_scatterBrain_locs_edges(nodeLocations, edge_X, edge_Y, edge_Z, slice);
 
hold on;
if strcmpi(slice, 'axial')
    scatter3(ax1, nodeLocations(:,1), nodeLocations(:,2), max(nodeLocations(:,3))*ones(N,1), ...
             nodeSizes, node_cval, 'filled', 'markeredgecolor', 'k');
else
    scatter3(ax1, nodeLocations(:,1), nodeLocations(:,2), nodeLocations(:,3), ...
         nodeSizes, node_cval, 'filled', 'markeredgecolor', 'k');
end
plot3(edges(:,1), edges(:,2), edges(:,3), 'color', 'k', 'linewidth', edge_width);
hold off;
text(nodeLocations(node_interests(1),1)+6, nodeLocations(node_interests(1),2), max(nodeLocations(:,3)), 'region $i$', ...
    'horizontalalignment', 'left', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(nodeLocations(node_interests(2),1)+6, nodeLocations(node_interests(2),2), max(nodeLocations(:,3)), 'region $j$', ...
    'horizontalalignment', 'left', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(-38, -46, 'LH', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, ...
    'fontweight', 'b')
text(38, -46, 'RH', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, ...
    'fontweight', 'b')
colormap(ax1, node_cmap)
view(ax1, 2)
axis(ax1, 'equal')
set(ax1, 'visible', 'off')

% =========================================================================
% B: model schematic
% =========================================================================
ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)*1.5 ax1.Position(2) 0.6 ax1.Position(4)]);
circle_locs_x = [0.25, 0.6; 
                 0.25, 0.6];
circle_locs_y = [0.72, 0.72; 
                 0.3, 0.3];
             
% annotation boxes
annotation('rectangle', [ax2.Position(1)+ax2.Position(3)*0.02 ax2.Position(2)-0.012 0.21 ax2.Position(4)*1.13], ...
    'linewidth', 2, 'color', node_interests_colors(1,:))
annotation('rectangle', [ax2.Position(1)+ax2.Position(3)*0.48 ax2.Position(2)-0.012 0.21 ax2.Position(4)*1.13], ...
    'linewidth', 2, 'color', node_interests_colors(2,:))
 
% model circles and lines
hold on;
plot([circle_locs_x(1,1),circle_locs_x(1,2)], [circle_locs_y(1,1),circle_locs_y(1,2)], 'k-', 'linewidth', 3)
plot([1,1]*circle_locs_x(1,1)-0.02, [circle_locs_y(1,1),circle_locs_y(2,1)], 'k-', 'linewidth', 3)
plot([1,1]*circle_locs_x(1,1)+0.02, [circle_locs_y(1,1),circle_locs_y(2,1)], 'k-', 'linewidth', 3)
plot([1,1]*circle_locs_x(1,2)-0.02, [circle_locs_y(1,2),circle_locs_y(2,2)], 'k-', 'linewidth', 3)
plot([1,1]*circle_locs_x(1,2)+0.02, [circle_locs_y(1,2),circle_locs_y(2,2)], 'k-', 'linewidth', 3)
plot(circle_locs_x(1,1), circle_locs_y(1,1)*1.19, 'ko', 'markersize', 25, 'linewidth', 3, 'markerfacecolor', 'none')
plot(circle_locs_x(1,2), circle_locs_y(1,2)*1.19, 'ko', 'markersize', 25, 'linewidth', 3, 'markerfacecolor', 'none')
plot(circle_locs_x(2,1), circle_locs_y(2,1)*0.57, 'ko', 'markersize', 25, 'linewidth', 3, 'markerfacecolor', 'none')
plot(circle_locs_x(2,2), circle_locs_y(2,2)*0.57, 'ko', 'markersize', 25, 'linewidth', 3, 'markerfacecolor', 'none')
plot(circle_locs_x(1,1), circle_locs_y(1,1), 'ko', 'markersize', 35, 'markerfacecolor', 'w')
plot(circle_locs_x(2,1), circle_locs_y(2,1), 'ko', 'markersize', 35, 'markerfacecolor', 'w')
plot(circle_locs_x(1,2), circle_locs_y(1,2), 'ko', 'markersize', 35, 'markerfacecolor', 'w')
plot(circle_locs_x(2,2), circle_locs_y(2,2), 'ko', 'markersize', 35, 'markerfacecolor', 'w')
hold off;
set(ax2, 'xtick', [], 'ytick', [], 'xlim', [0 1], 'ylim', [0 1])
 
% annotation lines
% wEE
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,1)*1.165 ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,1)*1.155], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,1)*1.22 ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,1)*1.13], ...
    'linestyle', 'none', 'headlength', 10, 'headwidth', 12)
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,2)*1.07 ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,2)*1.065], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,2)*1.22 ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,2)*1.13], ...
    'linestyle', 'none', 'headlength', 10, 'headwidth', 12)
% wII
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,1)*0.83 ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,1)*0.825], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,1)*0.62 ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,1)*0.69], ...
    'linestyle', 'none', 'headlength', 10, 'headwidth', 12)
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,2)*0.931 ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,2)*0.929], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,2)*0.62 ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,2)*0.69], ...
    'linestyle', 'none', 'headlength', 10, 'headwidth', 12)
% wEI 
annotation('arrow', ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,1)*0.917*[1,1], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,1)*0.8 ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,1)*0.83], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
annotation('arrow', ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,2)*0.965*[1,1], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,2)*0.8 ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,2)*0.83], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
% wIE 
annotation('arrow', ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,1)*1.08*[1,1], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,1)*1.42 ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,1)*1.39], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
annotation('arrow', ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,2)*1.032*[1,1], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,2)*1.42 ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,2)*1.39], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
% PE
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,1)*0.42 ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,1)*0.75], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,1)*1.1 ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,1)*1.04], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,2)*1.24 ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,2)*1.1], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,2)*1.1 ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,2)*1.04], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
% sigma_E
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,1)*0.42 ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,1)*0.75], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,1)*0.9 ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,1)*0.96], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,2)*1.24 ax2.Position(1)+ax2.Position(3)*circle_locs_x(1,2)*1.1], [ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,2)*0.9 ax2.Position(2)+ax2.Position(4)*circle_locs_y(1,2)*0.96], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
% sigma_I
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,1)*0.42 ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,1)*0.75], ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,1)*[1,1], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
annotation('arrow', [ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,2)*1.24 ax2.Position(1)+ax2.Position(3)*circle_locs_x(2,2)*1.1], ax2.Position(2)+ax2.Position(4)*circle_locs_y(2,2)*[1,1], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12)
 
% texts
text(0.1, 1.15, 'region $i$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(0.75, 1.15, 'region $j$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(1,1), circle_locs_y(1,1), '$E$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(1,2), circle_locs_y(1,2), '$E$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(2,1), circle_locs_y(2,1), '$I$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(2,2), circle_locs_y(2,2), '$I$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(1,1), circle_locs_y(1,1)+0.3, '$w_{EE}$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(1,2), circle_locs_y(1,2)+0.3, '$w_{EE}$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(2,1), circle_locs_y(2,1)-0.28, '$w_{II}$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(2,2), circle_locs_y(2,2)-0.28, '$w_{II}$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(1,1)-0.08, mean(circle_locs_y(:,1)), '$w_{EI}$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(1,2)-0.08, mean(circle_locs_y(:,2)), '$w_{EI}$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(1,1)+0.08, mean(circle_locs_y(:,1)), '$w_{IE}$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(1,2)+0.08, mean(circle_locs_y(:,2)), '$w_{IE}$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(1,1)-0.18, circle_locs_y(1,1)+0.1, '$P_E$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(1,2)+0.18, circle_locs_y(1,2)+0.1, '$P_E$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(1,1)-0.18, circle_locs_y(1,1)-0.1, '$D_E$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(1,2)+0.18, circle_locs_y(1,2)-0.1, '$D_E$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(2,1)-0.18, circle_locs_y(2,1), '$D_I$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(circle_locs_x(2,2)+0.18, circle_locs_y(2,2), '$D_I$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(mean(circle_locs_x(1,:)), circle_locs_y(1,2)-0.07, '$A_{ij}$', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, 'interpreter', 'latex')
axis off

% =========================================================================
% C: tuning curve
% =========================================================================
for type_ind = 1:length(types)
    S_temp = data_SuppFigure8.tuningcurve.SE{type_ind};
    stats_temp = data_SuppFigure8.tuningcurve.stats_SE{type_ind};
    
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    gain_cmap = repmat(cmap(type_ind,:),length(stats_temp.max),1);
    brighten_vec = linspace(0.8,0,length(stats_temp.max));
 
    [~,xval_transition_ind] = sort(stats_temp.xval_transition, 'ascend');
    
    subplot(3,length(types),2+type_ind)
    ax3 = gca;
    hold on;
    for ii=1:length(stats_temp.max)
        plot(data_SuppFigure8.tuningcurve.wEE, S_temp(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
    end
    plot([min(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition, S_transition], 'k-', 'linewidth', 2)
    plot([min(stats_temp.xval_transition), min(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    plot([max(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    hold off;
    text(max(stats_temp.xval_transition)*1.05, S_transition, ...
        ['width = ', sprintf('%.2f', stats_temp.dx_transition)])
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 0.5], ...
                'xtick', 0:5:30, 'xlim', [0 30]);
    xlabel({'global excitatory'; 'recurrent strength, $w_{EE}$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'excitatory'; 'firing rate, $S_E$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax3_image= axes('Position', [ax3.Position(1)+0.01 ax3.Position(2)+ax3.Position(4)*0.75 image_width image_width]);
    imshow(bw);
    colormap(ax3_image, flipud(gray))
end
 
% =========================================================================
% D: dynamic range distribution
% =========================================================================
subplot(3,length(types),2+((length(types)+1):length(types)*2))
ax4 = gca;
data_violin = [];
for type_ind = 1:length(types)
    data_violin = utils.padconcatenation(data_violin, data_SuppFigure8.tuningcurve.stats_SE{type_ind}.(stats_to_plot)-mean(data_SuppFigure8.tuningcurve.stats_SE{type_ind}.(stats_to_plot)), 2);
end
    
violins = utils.violinplot(data_violin, {});
for type_ind = 1:length(types)
    violins(type_ind).MedianPlot.SizeData = 50;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 10;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.03;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
    text(type_ind, 0.55, titles{type_ind}, 'color', cmap(type_ind,:), ...
    'fontsize', fontsize_axis, 'fontweight', 'b', 'verticalalignment', 'bottom')
    text(type_ind, 0.55, ['(\sigma = ', sprintf('%.2f)', std(data_SuppFigure8.tuningcurve.stats_SE{type_ind}.(stats_to_plot)-mean(data_SuppFigure8.tuningcurve.stats_SE{type_ind}.(stats_to_plot))))], 'color', cmap(type_ind,:), ...
        'fontsize', fontsize_axis, 'fontweight', 'b', 'horizontalalignment', 'left', 'verticalalignment', 'top')
end
set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0.5, length(types)+0.5], 'ylim', [-0.7 0.7], 'xtick', []);
ylabel([stats_to_plot_name, ' (mean-subtracted)'], 'fontsize', fontsize_label, 'interpreter', 'latex')
view(90, 90);
box off
ax4.XAxis.Visible = 'off';
 
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    
    ax4_image= axes('Position', [ax4.Position(1)+ax4.Position(3)*0.81 ax4.Position(2)+ax4.Position(4)*0.63-0.11*(type_ind-1) image_width image_width]);
    imshow(bw);
    colormap(ax4_image, flipud(gray))
end
 
%%% panel letters
annotation(fig, 'textbox', [0.07, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.07, 0.67, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.07, 0.35, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% SUPPLEMENTARY FIGURE 9

% load all data relevant to Supplementary Figure 9
data_SuppFigure9 = load(sprintf('%s/SuppFigure9.mat', data_foldername)); 

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];
S_transition = 0.15;
stats_to_plot = 'dynamic_range';
image_width = 0.15;
 
fig = figure('Position', [200 200 800 300]);
for type_ind = 1:length(types)
    type = types{type_ind};
    nodeLocations = data_SuppFigure9.coordinates{type_ind};
    
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;
    
    subplot(1,2,type_ind)
    ax1 = gca;
    data_to_plot_x = nodeLocations(:,2);
    data_to_plot_y = zscore(data_SuppFigure9.tuningcurve.stats{type_ind}.(stats_to_plot));
    hold on;
    plot(data_to_plot_x, zscore(data_to_plot_y), '.', 'markersize', 20, 'color', cmap(type_ind,:))
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,zscore(data_to_plot_y),1), data_to_plot_x), ...
        'k-', 'linewidth', 2)
    hold off;
    xline(0, 'k:');
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xtick', [], 'xticklabel', {});
    ylabel('dynamic range (z score)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    text(max(get(gca, 'xlim')), min(get(gca, 'ylim'))*0.8, sprintf('r = %.2g', rho), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
    text(max(get(gca, 'xlim')), min(get(gca, 'ylim'))*0.95, sprintf('p = %.2g', pval), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
    title(titles{type_ind}, 'fontsize', fontsize_label)
    annot = annotation('textarrow', [ax1.Position(1)+ax1.Position(3)*0.65 ax1.Position(1)+ax1.Position(3)*0.35], (ax1.Position(2)-0.04)*ones(1,2), 'String', 'anterior', ...
               'fontsize', fontsize_label, 'horizontalalignment', 'center', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5);
    annotation('textarrow', [annot.X(2) annot.X(1)], [annot.Y(1) annot.Y(2)], 'String', 'posterior', ...
               'fontsize', fontsize_label, 'horizontalalignment', 'center', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5)
           
           
    ax1_image= axes('Position', [ax1.Position(1)-0.04 ax1.Position(2)+ax1.Position(4)*0.8 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
end

%% SUPPLEMENTARY FIGURE 10

markersize = 10;
node_interest_markersize = 50;
nodeLocations = region_centroids;
 
fig = figure('Position', [200 200 800 400]);
% =========================================================================
% axial view
% =========================================================================
slice = 'axial';
for net = 1:7
    if net==1
        ax1 = axes('Position', [0.04 0.6 0.15 0.23]);
    else
        ax1 = axes('Position', [ax1.Position(1)+0.135 ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
    end
    cvals = double(RSN_indices==net);
    net_ind = find(cvals);
    node_interest = net_ind;
    node_interest_data = net*ones(length(net_ind),1);
    node_interest_cmap = Yeo_7net_colormap(net,:);
    ax1 = utils.draw_brain_nodeHighlight(ax1, node_interest, nodeLocations, ...
                                 markersize, node_interest_data, node_interest_cmap, ...
                                 node_interest_markersize, slice);
    if net==1
        annotation('textbox', [ax1.Position(1)-0.03, ax1.Position(2)-0.05, 0.1, 0.1], 'string', 'LH', ...
           'fontsize', 6, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
           'verticalalignment', 'middle', 'fontweight', 'b')
        annotation('textbox', [ax1.Position(1)+ax1.Position(3)-0.07, ax1.Position(2)-0.05, 0.1, 0.1], 'string', 'RH', ...
           'fontsize', 6, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
           'verticalalignment', 'middle', 'fontweight', 'b')
    end
    text(0, 70, Yeo_names{net}, 'horizontalalignment', 'center', 'fontsize', fontsize_label, ...
        'fontweight', 'b')
end
 
% =========================================================================
% sagittal left view
% =========================================================================
slice = 'sagittal_left';
for net = 1:7
    if net==1
        ax2 = axes('Position', [0.06 ax1.Position(2)-0.26 0.11 0.2]);
    else
        ax2 = axes('Position', [ax2.Position(1)+0.135 ax2.Position(2) ax2.Position(3) ax2.Position(4)]);
    end
    cvals = double(RSN_indices==net);
    net_ind = find(cvals);
    node_interest = net_ind;
    node_interest_data = net*ones(length(net_ind),1);
    node_interest_cmap = Yeo_7net_colormap(net,:);
    ax2 = utils.draw_brain_nodeHighlight(ax2, node_interest, nodeLocations, ...
                                 markersize, node_interest_data, node_interest_cmap, ...
                                 node_interest_markersize, slice);
    if net==1
        annotation('textbox', [ax2.Position(1)-0.04, ax2.Position(2)-0.02, 0.1, 0.1], 'string', 'LH', ...
           'fontsize', 6, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
           'verticalalignment', 'middle', 'fontweight', 'b')
    end
end
 
% =========================================================================
% sagittal right view
% =========================================================================
slice = 'sagittal_right';
for net = 1:7
    if net==1
        ax3 = axes('Position', [0.06 ax2.Position(2)-0.24 ax2.Position(3) ax2.Position(4)]);
    else
        ax3 = axes('Position', [ax3.Position(1)+0.135 ax3.Position(2) ax3.Position(3) ax3.Position(4)]);
    end
    cvals = double(RSN_indices==net);
    net_ind = find(cvals);
    node_interest = net_ind;
    node_interest_data = net*ones(length(net_ind),1);
    node_interest_cmap = Yeo_7net_colormap(net,:);
    ax3 = utils.draw_brain_nodeHighlight(ax3, node_interest, nodeLocations, ...
                                 markersize, node_interest_data, node_interest_cmap, ...
                                 node_interest_markersize, slice);
    if net==1
        annotation('textbox', [ax3.Position(1)+ax3.Position(3)-0.06, ax3.Position(2)-0.02, 0.1, 0.1], 'string', 'RH', ...
           'fontsize', 6, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
           'verticalalignment', 'middle', 'fontweight', 'b')
    end
end

%% SUPPLEMENTARY FIGURE 11

% load all data relevant to Supplementary Figure 11
data_SuppFigure11 = load(sprintf('%s/SuppFigure11.mat', data_foldername)); 

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];
nodeLocations = region_centroids;
image_width = 0.05;
N = 114;

fig = figure('Position', [400 100 800 750]);
for type_ind = 1:length(types)
    
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
 
    bw = image_to_plot>0;
    
    ax_positions = [0.15+0.43*(type_ind-1) 0.7 0.35 0.25];
    
    % =====================================================================
    % A: path length distribution
    % =====================================================================
    if type_ind==1
        ax1 = axes('Position', [ax_positions(1)+ax_positions(3)*0.6 ax_positions(2) ax_positions(3) ax_positions(4)]);
 
        data_violin = utils.padconcatenation(data_SuppFigure11.tuningcurve.stats{1}.pl, data_SuppFigure11.tuningcurve.stats{2}.pl, 2);
        violins = utils.violinplot(data_violin, {});
        for type_ind_2 = 1:length(types)
            violins(type_ind_2).MedianPlot.SizeData = 30;
            violins(type_ind_2).ViolinColor = cmap(type_ind_2,:);
            violins(type_ind_2).ViolinAlpha = 0.2;
            violins(type_ind_2).ScatterPlot.SizeData = 10;
            violins(type_ind_2).ScatterPlot.MarkerFaceAlpha = 1;
            violins(type_ind_2).BoxColor = [0 0 0];
            violins(type_ind_2).BoxWidth = 0.01;
            violins(type_ind_2).WhiskerPlot.LineStyle = 'none';
 
%             text(type_ind_2, max(data_violin(:,type_ind_2))*1.01, ['\sigma = ', sprintf('%.2f', std(data_violin(:,type_ind_2)))], 'color', cmap(type_ind_2,:), ...
%             'fontsize', fontsize_label, 'fontweight', 'b', 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
        end
 
        set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xticklabel', titles)
        ylabel('regional path length', 'fontsize', fontsize_label, 'interpreter', 'latex')
        
        for type_ind_2 = 1:length(types)
            if type_ind_2==1
                image_to_plot = human_female;
                yloc = 0.35;
            elseif type_ind_2==2
                image_to_plot = chimpanzee;
                yloc = 0.9;
            end
            bw2 = image_to_plot>0;
 
            ax1_image = axes('Position', [ax1.Position(1)+0.02+(type_ind_2-1)*0.17 ax1.Position(2)+ax1.Position(4)*yloc image_width image_width]);
            imshow(bw2);
            colormap(ax1_image, flipud(gray))
        end
    end
    
    stats_temp = data_SuppFigure11.tuningcurve.stats{type_ind};
    % =====================================================================
    % B: path length vs dynamic range   
    % =====================================================================
    ax2 = axes('Position', [ax_positions(1) ax1.Position(2)-ax1.Position(4)*1.3 ax_positions(3) ax_positions(4)]);
    data_to_plot_x = zscore(stats_temp.dynamic_range);
    data_to_plot_y = stats_temp.pl;
    hold on;
    plot(data_to_plot_x, data_to_plot_y, '.', 'color', cmap(type_ind,:), 'markersize', 12)
%     plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
%         'k-', 'linewidth', 2)
    hold off;
    set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'box', 'off', 'ylim', [15 40])
    xlabel('dynamic range (z score)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    ylabel('regional path length', 'fontsize', fontsize_label, 'interpreter', 'latex')
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'spearman');
    text(max(data_to_plot_x)*1, 39, ['\rho ', sprintf('= %.2g', rho)], ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
    text(max(data_to_plot_x)*1, 36, sprintf('p = %.2g', pval), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
 
    ax2_image = axes('Position', [ax2.Position(1)+0.01 ax2.Position(2)+ax2.Position(4)*0.05 image_width image_width]);
    imshow(bw);
    colormap(ax2_image, flipud(gray)) 
    
    % =====================================================================
    % C: path length spatial organization
    % =====================================================================
    pl_cmap = cbrewer('seq', 'YlOrBr', 125, 'pchip');
    pl_cmap = flipud(pl_cmap(end-N+1:end,:));
    
    ax3 = axes('Position', [ax_positions(1) ax2.Position(2)-ax2.Position(4)*1.12 ax_positions(3) ax_positions(4)]);
    ax3_positions = ax3.Position;
    set(ax3, 'visible', 'off')
    
    ax3_1 = axes('Position', [ax3_positions(1)+0.02 ax3_positions(2)-0.05 0.15 0.25]);
    [ax3_1, ~] = utils.draw_scatterBrain(ax3_1, nodeLocations, stats_temp.pl, 30, 'axial');
    colormap(ax3_1, pl_cmap)
    annotation('textbox', [ax3_1.Position(1)-0.04, ax3_1.Position(2)+ax3_1.Position(4)*0.0, 0.1, 0.1], 'string', 'LH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annotation('textbox', [ax3_1.Position(1)+ax3_1.Position(3)-0.06, ax3_1.Position(2)+ax3_1.Position(4)*0.0, 0.1, 0.1], 'string', 'RH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annot = annotation('textarrow', ax3_1.Position(1)*ones(1,2), [ax3_1.Position(2)+ax3_1.Position(4)*0.55 ax3_1.Position(2)+ax3_1.Position(4)*0.4], 'String', 'a', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5);
    annotation('textarrow', [annot.X(1) annot.X(2)], [annot.Y(2) annot.Y(1)], 'String', 'p', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5)
   
    ax3_1_image = axes('Position', [ax3_1.Position(1)-0.02 ax3_1.Position(2)+ax3_1.Position(4)*0.7 image_width image_width]);
    imshow(bw);
    colormap(ax3_1_image, flipud(gray))
 
    ax3_2 = axes('Position', [ax3_1.Position(1)+ax3_1.Position(3)*1.2 ax3_1.Position(2)+ax3_1.Position(4)*0.22 0.10 0.25]);
    [ax3_2, ~] = utils.draw_scatterBrain(ax3_2, nodeLocations, stats_temp.pl, 20, 'sagittal_left');
    colormap(ax3_2, pl_cmap)
    annotation('textbox', [ax3_2.Position(1)-0.04, ax3_2.Position(2)+ax3_2.Position(4)*0.18, 0.1, 0.1], 'string', 'LH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annot = annotation('textarrow', [ax3_2.Position(1)+ax3_2.Position(3)*0.7 ax3_2.Position(1)+ax3_2.Position(3)*0.95], (ax3_2.Position(2)+ax3_2.Position(4)*0.33)*ones(1,2), 'String', 'a', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5);
    annotation('textarrow', [annot.X(2) annot.X(1)], [annot.Y(1) annot.Y(2)], 'String', 'p', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5)
           
    ax3_3 = axes('Position', [ax3_2.Position(1) ax3_2.Position(2)-0.1 ax3_2.Position(3) ax3_2.Position(4)]);
    [ax3_3, ~] = utils.draw_scatterBrain(ax3_3, nodeLocations, stats_temp.pl, 20, 'sagittal_right');
    colormap(ax3_3, pl_cmap)
    annotation('textbox', [ax3_3.Position(1)+ax3_3.Position(3)*0.45, ax3_3.Position(2)+ax3_3.Position(4)*0.18, 0.1, 0.1], 'string', 'RH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annot = annotation('textarrow', [ax3_3.Position(1)+ax3_3.Position(3)*0.3 ax3_3.Position(1)+ax3_3.Position(3)*0.05], (ax3_3.Position(2)+ax3_3.Position(4)*0.33)*ones(1,2), 'String', 'a', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5);
    annotation('textarrow', [annot.X(2) annot.X(1)], [annot.Y(1) annot.Y(2)], 'String', 'p', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5)
           
    cbar = colorbar(ax3_3,'southoutside');
    ylabel(cbar, 'regional path length', 'fontsize', fontsize_label, 'interpreter', 'latex')
    set(cbar, 'fontsize', fontsize_axis, 'ticklength', 0.01, ...
        'position', [ax3_positions(1)+ax3_positions(3)*0.36, ax3_3.Position(2)-0.05, ax3_3.Position(3)*0.9, 0.01], ...
        'ytick', [])
    annotation(fig, 'textbox', [cbar.Position(1)-0.045, cbar.Position(2)*1, 0.04, 0.01], 'string', 'LOW', 'edgecolor', 'none', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    annotation(fig, 'textbox', [cbar.Position(1)+cbar.Position(3)+0.01, cbar.Position(2)*1, 0.04, 0.01], 'string', 'HIGH', 'edgecolor', 'none', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')        
end

%%% panel letters
annotation(fig, 'textbox', [0.3, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.1, 0.66, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.1, 0.29, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
    
%% SUPPLEMENTARY FIGURE 12

% load all data relevant to Supplementary Figure 12
data_SuppFigure12 = load(sprintf('%s/SuppFigure12.mat', data_foldername));

cmap = color_blue;
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';
image_width = 0.1;
 
fig = figure('Position', [200 200 400 450]);

stats_temp = data_SuppFigure12.tuningcurve.stats;
image_to_plot = macaque;
bw = image_to_plot>0;
gain_cmap = repmat(cmap,length(stats_temp.max),1);
brighten_vec = linspace(0.8,0,length(stats_temp.max));
 
[~,xval_transition_ind] = sort(stats_temp.xval_transition, 'ascend');

% =========================================================================
% A: tuning curve
% =========================================================================
subplot(2,1,1)
ax1 = gca;
set(ax1, 'Position', [ax1.Position(1)+0.035 ax1.Position(2)-0.14 ax1.Position(3)+0.01 ax1.Position(4)+0.14])
hold on;
for ii=1:length(stats_temp.max)
    plot(data_SuppFigure12.tuningcurve.w, data_SuppFigure12.tuningcurve.S(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
end
plot([min(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
     [S_transition, S_transition], 'k-', 'linewidth', 2)
plot([min(stats_temp.xval_transition), min(stats_temp.xval_transition)], ...
     [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
plot([max(stats_temp.xval_transition), max(stats_temp.xval_transition)], ...
     [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
hold off;
text(max(stats_temp.xval_transition)*1.05, S_transition, ...
    ['width = ', sprintf('%.2f', stats_temp.dx_transition)])
set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
            'xtick', 0:0.5:2, 'xlim', [0 2]);
xlabel('global recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
title('macaque', 'fontsize', fontsize_label)
 
ax1_image= axes('Position', [ax1.Position(1)+0.02 ax1.Position(2)+ax1.Position(4)*0.795 image_width image_width]);
imshow(bw);
colormap(ax1_image, flipud(gray))
    
box_dodge_amounts = linspace(0.2,2.5,2);
dot_dodge_amounts = box_dodge_amounts+0.2;

% =========================================================================
% B: dynamic range distribution
% =========================================================================
ax2 = axes('Position', [ax1.Position(1) 0.12 ax1.Position(3) 0.15]);
 
data_violin = stats_temp.(stats_to_plot);
violins = utils.violinplot(data_violin, {});
violins.MedianPlot.SizeData = 50;
violins.ViolinColor = cmap;
violins.ViolinAlpha = 0.2;
violins.ScatterPlot.SizeData = 10;
violins.ScatterPlot.MarkerFaceAlpha = 1;
violins.BoxColor = [0 0 0];
violins.BoxWidth = 0.02;
violins.WhiskerPlot.LineStyle = 'none';
text(1.25, min(data_violin), ['\sigma = ', sprintf('%.2f', std(data_violin))], 'color', 'k', ...
    'fontsize', fontsize_axis, 'fontweight', 'b', 'horizontalalignment', 'left')
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0.5, 1.5], 'xtick', []);
ylabel(stats_to_plot_name, 'fontsize', fontsize_label, 'interpreter', 'latex')
view(90, 90);
box off
ax2.XAxis.Visible = 'off';
 
%%% panel letters
annotation(fig, 'textbox', [0.05, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.05, 0.32, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% SUPPLEMENTARY FIGURE 13

% load all data relevant to Supplementary Figure 13
data_SuppFigure13 = load(sprintf('%s/SuppFigure13.mat', data_foldername)); 

fig = figure('Position', [200 200 700 300]);

% =========================================================================
% left: time series
% =========================================================================
ax1_left = axes('Position', [0.07 0.15 0.38 0.8]);
plot(data_SuppFigure13.timeseries.time, data_SuppFigure13.timeseries.S, 'k', 'linewidth', 1.5)
set(ax1_left, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ytick', []);
xlabel('time (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('synaptic gating, $S$', 'fontsize', fontsize_label, 'interpreter', 'latex')
box off;

% =========================================================================
% right: autocorrelation fitting method
% =========================================================================
ax1_right = axes('Position', [ax1_left.Position(1)+ax1_left.Position(3)+0.15 ax1_left.Position(2) ax1_left.Position(3) ax1_left.Position(4)]);
hold on;
plot(data_SuppFigure13.autocorrelation.lags, data_SuppFigure13.autocorrelation.acf_data/max(data_SuppFigure13.autocorrelation.acf_data), 'k', 'linewidth', 1.5)
plot(data_SuppFigure13.autocorrelation.lags, data_SuppFigure13.autocorrelation.acf_fit/max(data_SuppFigure13.autocorrelation.acf_fit), 'r--', 'linewidth', 2)
hold off;
set(ax1_right, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 5], 'ytick', []);
xlabel('time lag, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('autocorrelation', 'fontsize', fontsize_label, 'interpreter', 'latex')
leg = legend('{\rm data}', 'fit = $c_1e^{-t/\tau} + c_2$', 'location', 'northeast', 'interpreter', 'latex');
set(leg, 'fontsize', fontsize_label, 'box', 'off', 'numcolumns', 1)
box off;
annotation('arrow', [ax1_left.Position(1)+ax1_left.Position(3)+0.02 ax1_right.Position(1)-0.05], [ax1_left.Position(2)+ax1_left.Position(4)*0.5 ax1_left.Position(2)+ax1_left.Position(4)*0.5])

%% SUPPLEMENTARY FIGURE 14

% load all data relevant to Supplementary Figure 14
data_SuppFigure14 = load(sprintf('%s/SuppFigure14.mat', data_foldername));

types = {'human', 'chimp', 'macaque', 'marmoset'};
titles = {'human', 'chimpanzee', 'macaque', 'marmoset'};

cmap = [color_brown; color_green; color_blue; color_purple];

w_interest_list = [0.45, 0.45, 0.32, 0.45];
w_interest_ind = dsearchn(data_SuppFigure14.tuningcurve.w', w_interest_list');
 
image_width = 0.06;
 
fig = figure('Position', [200 200 600 600]);
 
ax0 = axes('Position', [0.13 0.6 0.35 0.3]);
ax0_position = ax0.Position;
set(ax0,'visible','off')

for type_ind=[3,4]
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    elseif type_ind==3
        image_to_plot = macaque;
    elseif type_ind==4
        image_to_plot = marmoset;
    end
    bw = image_to_plot>0;
    
    % =========================================================================
    % A: correlation of timescales and dynamic range
    % =========================================================================
    ax1 = axes('Position', [ax0_position(1)+(ax0_position(3)+0.1)*(type_ind-1-2) ax0_position(2) ax0_position(3) ax0_position(4)]);
    data_to_plot_x = tiedrank(zscore(data_SuppFigure14.tuningcurve.stats{type_ind}.dynamic_range));
    data_to_plot_y = tiedrank(data_SuppFigure14.tuningcurve.stats{type_ind}.tau_neural(:,w_interest_ind(type_ind)));
    hold on;
    plot(data_to_plot_x, data_to_plot_y, '.', 'color', cmap(type_ind,:), 'markersize', 20)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'k-', 'linewidth', 2);
    hold off;
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'box', 'off', 'xlim', [1 10*ceil(length(data_to_plot_x)/10)], 'ylim', [1 10*ceil(length(data_to_plot_x)/10)])
    xlabel('rank of dynamic range', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1+2
        ylabel('rank of timescale', 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    if type_ind==3
        text(max(get(ax1,'xlim'))+8, max(get(ax1,'ylim'))*0.38, sprintf('r = %.2g', rho), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
        text(max(get(ax1,'xlim'))+8, max(get(ax1,'ylim'))*0.30, sprintf('p = %.2g', pval), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
    else
        text(max(get(ax1,'xlim'))+8, max(get(ax1,'ylim'))*0.26, sprintf('r = %.2g', rho), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
        text(max(get(ax1,'xlim'))+8, max(get(ax1,'ylim'))*0.19, sprintf('p = %.2g', pval), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
    end

    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax1_image = axes('Position', [ax1.Position(1)+0.01 ax1.Position(2)+ax1.Position(4)*0.8 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))

    % =========================================================================
    % B: human - macaque/marmoset whole-brain accuracy
    % =========================================================================
    ax2 = axes('Position', [ax1.Position(1) 0.1 ax0_position(3) ax0_position(4)]);
    data_to_plot_x = data_SuppFigure14.decision.time;
    data_to_plot_y = mean(data_SuppFigure14.decision.accuracy{1},1)-mean(data_SuppFigure14.decision.accuracy{type_ind},1);
    [~, min_diff_ind] = min(data_to_plot_y);
    hold on;
    plot(data_to_plot_x, data_to_plot_y, 'k-', 'linewidth', 2)
    plot(data_to_plot_x, 0*(ones(size(data_to_plot_x))), 'k:');
    hold off;

    set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 1.5], 'ylim', [-8, 12])
    hold on;
    plot(data_to_plot_x(min_diff_ind)*ones(1,2), get(gca,'ylim'), 'k--');
    hold off;
    xlabel('time (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    text(data_to_plot_x(min_diff_ind), min(get(gca,'ylim')), '$t_{\rm min}$', 'interpreter', 'latex', ...
         'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontsize', fontsize_label)
    if type_ind==3
        ylabel({'difference in'; 'whole-brain accuracy ($\%$)'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(sprintf('%s – %s', titles{1}, titles{type_ind}), 'fontsize', fontsize_label)
end

%%% titles
annotation(fig, 'textbox', [0.1, 0.96, 0.8, 0.01], 'string', 'neural timescales', 'edgecolor', 'none', ...
        'fontsize', fontsize_label, 'fontweight', 'b', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
annotation(fig, 'textbox', [0.1, 0.47, 0.8, 0.01], 'string', 'computational capacity', 'edgecolor', 'none', ...
        'fontsize', fontsize_label, 'fontweight', 'b', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
      
%%% panel letters
annotation(fig, 'textbox', [0.03, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.03, 0.49, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
    
%% SUPPLEMENTARY FIGURE 15

% load all data relevant to Supplementary Figure 15
data_SuppFigure15 = load(sprintf('%s/SuppFigure15.mat', data_foldername));

nodeLocations = region_centroids;
markersize = 50;

FN = 7;
data_to_plot_x = data_SuppFigure15.decision.time;
data_to_plot_y = mean(data_SuppFigure15.decision.accuracy{1}(RSN_indices==FN,:),1)-mean(data_SuppFigure15.decision.accuracy{2}(RSN_indices==FN,:),1);

fig = figure;
% =========================================================================
% human - chimpanzee network accuracy
% =========================================================================
ax1 = axes('Position', [0.12 0.13 0.6 0.8]);
hold on;
plot(data_to_plot_x, data_to_plot_y, '-', 'linewidth', 2, 'color', Yeo_7net_colormap(FN,:))
hold off;
yline(0, 'k:');
set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 1.5])
xlabel('time (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel({'human -- chimpanzee'; 'network accuracy ($\%$)'}, 'fontsize', fontsize_label, 'interpreter', 'latex')

% =========================================================================
% axial view
% =========================================================================
ax1_1 = axes('Position', [ax1.Position(1)+ax1.Position(3)+0.04 ax1.Position(2)+ax1.Position(4)*0.56 0.2 0.4]);
cvals = double(RSN_indices==FN);
net_ind = find(cvals);
newnodeLocations = [nodeLocations(net_ind,1), nodeLocations(net_ind,2), max(nodeLocations(net_ind,3))*ones(length(net_ind),1)];
newcvals = ones(length(net_ind),1);
hold on;
[ax1_1, ~] = utils.draw_scatterBrain(ax1_1, nodeLocations, cvals, markersize, 'axial');
[ax1_1, ~] = utils.draw_scatterBrain(ax1_1, newnodeLocations, newcvals, markersize, 'axial');
hold off;
colormap(ax1_1, [0.85*ones(1,3); Yeo_7net_colormap(FN,:)])
annotation('textbox', [ax1_1.Position(1)-0.045, ax1_1.Position(2)+0.02, 0.1, 0.1], 'string', 'LH', ...
   'fontsize', fontsize_axis-2, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
   'verticalalignment', 'middle', 'fontweight', 'b')
annotation('textbox', [ax1_1.Position(1)+ax1_1.Position(3)-0.05, ax1_1.Position(2)+0.02, 0.1, 0.1], 'string', 'RH', ...
   'fontsize', fontsize_axis-2, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
   'verticalalignment', 'middle', 'fontweight', 'b')

% =========================================================================
% sagittal left view
% =========================================================================
ax1_2 = axes('Position', [ax1_1.Position(1) ax1_1.Position(2)-0.3 ax1_1.Position(3) ax1_1.Position(4)]);
cvals = double(RSN_indices==FN);
net_ind = find(cvals);
newnodeLocations = [min(nodeLocations(:,1))*ones(length(net_ind),1), nodeLocations(net_ind,2), nodeLocations(net_ind,3)];
newcvals = ones(length(net_ind),1);
hold on;
[ax1_2, ~] = utils.draw_scatterBrain(ax1_2, nodeLocations, cvals, markersize, 'sagittal_left');
[ax1_2, ~] = utils.draw_scatterBrain(ax1_2, newnodeLocations, newcvals, markersize, 'sagittal_left');
hold off;
colormap(ax1_2, get(ax1_1, 'Colormap'))
annotation('textbox', [ax1_2.Position(1)-0.045, ax1_2.Position(2)+0.07, 0.1, 0.1], 'string', 'LH', ...
   'fontsize', fontsize_axis-2, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
   'verticalalignment', 'middle', 'fontweight', 'b')

% =========================================================================
% sagittal right view
% =========================================================================
ax1_3 = axes('Position', [ax1_2.Position(1) ax1_2.Position(2)-0.25 ax1_2.Position(3) ax1_2.Position(4)]);
cvals = double(RSN_indices==FN);
net_ind = find(cvals);
newnodeLocations = [max(nodeLocations(:,1))*ones(length(net_ind),1), nodeLocations(net_ind,2), nodeLocations(net_ind,3)];
newcvals = ones(length(net_ind),1);
hold on;
[ax1_3, ~] = utils.draw_scatterBrain(ax1_3, nodeLocations, cvals, markersize, 'sagittal_right');
[ax1_3, ~] = utils.draw_scatterBrain(ax1_3, newnodeLocations, newcvals, markersize, 'sagittal_right');
hold off;
colormap(ax1_3, get(ax1_1, 'Colormap'))
annotation('textbox', [ax1_3.Position(1)+ax1_3.Position(3)-0.05, ax1_3.Position(2)+0.07, 0.1, 0.1], 'string', 'RH', ...
   'fontsize', fontsize_axis-2, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
   'verticalalignment', 'middle', 'fontweight', 'b')
