% generate_paper_figures_long.m
% 
% Matlab code to generate the main figures of the paper

%% LOAD CONNECTOME MATRICES AND REGION LOCATIONS/NAMES

load('data/connectome_human.mat')
load('data/connectome_chimp.mat')
load('data/connectome_macaque.mat')
load('data/connectome_marmoset.mat')
load('data/region_props_chimp.mat', 'region_centroids', 'region_names') 
load('data/RSN.mat', 'RSN_indices', 'RSN_names')

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

cb = cbrewer('qual', 'Set1', 8, 'pchip');
color_brown = cb(7,:);
color_green = cb(3,:);
color_blue = cb(2,:);
color_purple = cb(4,:);

% artworks
human_female = imread('artworks/human_female_2.png');
human_male = imread('artworks/human_male_2.png');
chimpanzee = imread('artworks/chimpanzee.png');
macaque = imread('artworks/macaque.png');
marmoset = imread('artworks/marmoset.png');

fontsize_axis = 10;
fontsize_label = 12;

data_foldername = 'data/data_figures_long';

%% FIGURE 1

% load all data relevant to Figure 1
data_Figure1 = load(sprintf('%s/Figure1.mat', data_foldername));

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];
matrix_cmap = [1,1,1; cbrewer('seq', 'Blues', 100, 'pchip')];

edge_alpha = 0.8;
edge_width = 1;
edge_colors = [cmap(1,:), edge_alpha;
               cmap(2,:), edge_alpha];
           
N = 114;
slice = 'axial';   
node_interests = [86,84];
frac = 0.1;
markersize = 40;
node_interests_colors = brighten(lines(3),0.3);
node_interests_colors = node_interests_colors(2:end,:);
node_cmap = [1, 1, 1; node_interests_colors];
node_cval = ones(N,1);
node_cval(node_interests(1)) = 2;
node_cval(node_interests(2)) = 3;
nodeLocations = region_centroids;
nodeSizes = (markersize/10)*ones(N,1);
nodeSizes(node_interests) = markersize;

fig = figure('Position', [200 200 800 400]);
% =========================================================================
% A: connectomes
% =========================================================================
for type_ind=1:length(types)
    type = types{type_ind};
    
    W = eval(sprintf('connectome_%s', type));
    [edge_X, edge_Y, edge_Z] = adjacency_plot_und(threshold_proportional(W, frac), nodeLocations);  % get all the edges
    
    [nodeLocations, edges] = utils.extract_scatterBrain_locs_edges(nodeLocations, edge_X, edge_Y, edge_Z, slice);
    
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;

    ax1 = axes('Position', [-0.05+0.26*(type_ind-1) 0.45 0.44 0.44]);
    hold on;
    if strcmpi(slice, 'axial')
        scatter3(ax1, nodeLocations(:,1), nodeLocations(:,2), max(nodeLocations(:,3))*ones(N,1), ...
                 nodeSizes, node_cval, 'filled', 'markeredgecolor', 'k');
    else
        scatter3(ax1, nodeLocations(:,1), nodeLocations(:,2), nodeLocations(:,3), ...
             	 nodeSizes, node_cval, 'filled', 'markeredgecolor', 'k');
    end
    plot3(edges(:,1), edges(:,2), edges(:,3), 'color', edge_colors(type_ind,:), 'linewidth', edge_width);
    hold off;
    text(0, 60, titles{type_ind}, 'horizontalalignment', 'center', 'fontsize', fontsize_label, ...
        'fontweight', 'b')
    text(nodeLocations(node_interests(1),1)+6, nodeLocations(node_interests(1),2), max(nodeLocations(:,3)), 'region $i$', ...
    'horizontalalignment', 'left', 'fontsize', fontsize_axis, 'interpreter', 'latex')
    text(nodeLocations(node_interests(2),1)+6, nodeLocations(node_interests(2),2), max(nodeLocations(:,3)), 'region $j$', ...
        'horizontalalignment', 'left', 'fontsize', fontsize_axis, 'interpreter', 'latex')
    text(-40, -41, 'LH', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, ...
        'fontweight', 'b')
    text(40, -41, 'RH', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, ...
        'fontweight', 'b')
    colormap(ax1, node_cmap)
    view(ax1, 2)
    axis(ax1, 'equal')
    set(ax1, 'visible', 'off')
    
    ax1_image = axes('Position', [ax1.Position(1)+0.08 ax1.Position(2)+ax1.Position(4)*0.8 0.1 0.1]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
end

% =========================================================================
% B left: model schematic
% =========================================================================
ax2 = axes('Position', [0.03 0.07 0.3 0.3]);
circle_locs_x = [0.25, 0.7];
circle_locs_y = [0.55, 0.55];
hold on;
plot(circle_locs_x, circle_locs_y, 'k-', 'linewidth', 3)
p = plot(circle_locs_x(1), circle_locs_y(1)+0.25, 'ko', 'markersize', 25, 'linewidth', 3, 'markerfacecolor', 'none');
plot(circle_locs_x(2), circle_locs_y(2)+0.25, 'ko', 'markersize', 25, 'linewidth', 3, 'markerfacecolor', 'none')
plot(circle_locs_x(1), circle_locs_y(1), 'ko', 'markersize', 40, 'markerfacecolor', node_interests_colors(1,:))
plot(circle_locs_x(2), circle_locs_y(2), 'ko', 'markersize', 40, 'markerfacecolor', node_interests_colors(2,:))
hold off;
set(ax2, 'xtick', [], 'ytick', [], 'xlim', [0 1], 'ylim', [0 1])

% annotation lines
utils.annotation_pinned('arrow', [circle_locs_x(1)+0.043 circle_locs_x(1)+0.025], [circle_locs_y(1)*1.37 circle_locs_y(1)*1.28], ...
    'linestyle', 'none', 'headlength', 10, 'headwidth', 12);
utils.annotation_pinned('arrow', [circle_locs_x(2)+0.043 circle_locs_x(2)+0.025], [circle_locs_y(2)*1.37 circle_locs_y(2)*1.28], ...
    'linestyle', 'none', 'headlength', 10, 'headwidth', 12);
utils.annotation_pinned('arrow', [circle_locs_x(1)-0.1 circle_locs_x(1)-0.02], [circle_locs_y(1)*0.35 circle_locs_y(1)*0.72], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12);
utils.annotation_pinned('arrow', [circle_locs_x(1)+0.1 circle_locs_x(1)+0.02], [circle_locs_y(1)*0.35 circle_locs_y(1)*0.72], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12);
utils.annotation_pinned('arrow', [circle_locs_x(2)-0.1 circle_locs_x(2)-0.02], [circle_locs_y(2)*0.35 circle_locs_y(2)*0.72], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12);
utils.annotation_pinned('arrow', [circle_locs_x(2)+0.1 circle_locs_x(2)+0.02], [circle_locs_y(2)*0.35 circle_locs_y(2)*0.72], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 12);

% texts
text(circle_locs_x(1), circle_locs_y(1), {'region'; '$i$'}, 'horizontalalignment', 'center', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(circle_locs_x(2), circle_locs_y(2), {'region'; '$j$'}, 'horizontalalignment', 'center', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(circle_locs_x(1), circle_locs_y(1)+0.48, '$w$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(2), circle_locs_y(1)+0.48, '$w$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(mean(circle_locs_x), circle_locs_y(1), '$A_{ij}$', 'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(1)+0.08, circle_locs_y(1)-0.36, '$D$', 'horizontalalignment', 'left', 'verticalalignment', 'top', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(2)+0.08, circle_locs_y(1)-0.36, '$D$', 'horizontalalignment', 'left', 'verticalalignment', 'top', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(1)-0.08, circle_locs_y(1)-0.36, '$I_0$', 'horizontalalignment', 'right', 'verticalalignment', 'top', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(2)-0.08, circle_locs_y(1)-0.36, '$I_0$', 'horizontalalignment', 'right', 'verticalalignment', 'top', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
axis off

% =========================================================================
% B right: time series
% =========================================================================
ax3 = axes('Position', [ax2.Position(1)+ax2.Position(3)*1.15 0.15 0.14 0.2]);
hold on;
plot(data_Figure1.timeseries.time, data_Figure1.timeseries.S(node_interests(1),:), 'color', node_interests_colors(1,:), 'linewidth', 1.5)
plot(data_Figure1.timeseries.time, data_Figure1.timeseries.S(node_interests(2),:), 'color', node_interests_colors(2,:), 'linewidth', 1.5)
hold off;
text(max(data_Figure1.timeseries.time), mean(data_Figure1.timeseries.S(node_interests(1),:)), 'region $i$', ...
    'horizontalalignment', 'left', 'fontsize', fontsize_axis, 'interpreter', 'latex')
text(max(data_Figure1.timeseries.time), mean(data_Figure1.timeseries.S(node_interests(2),:)), 'region $j$', ...
    'horizontalalignment', 'left', 'fontsize', fontsize_axis, 'interpreter', 'latex')
set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ytick', []);
xlabel('time (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')

% =========================================================================
% C: dynamic range calculation
% =========================================================================
node = 53;
ax4 = axes('Position', [ax3.Position(1)+ax3.Position(3)*2.1 0.15 0.3 0.77]);
hold on;
plot(data_Figure1.tuningcurve.w, data_Figure1.tuningcurve.S(node,:), 'k.');
plot(data_Figure1.tuningcurve.stats.xval_10(node)*ones(1,2), [0, 1], 'k:');
plot(data_Figure1.tuningcurve.stats.xval_90(node)*ones(1,2), [0, 1], 'k:');
plot(data_Figure1.tuningcurve.w, data_Figure1.tuningcurve.stats.Fval_10(node)*(ones(size(data_Figure1.tuningcurve.w))), 'k:');
plot(data_Figure1.tuningcurve.w, data_Figure1.tuningcurve.stats.Fval_90(node)*(ones(size(data_Figure1.tuningcurve.w))), 'k:');
hold off;
text(data_Figure1.tuningcurve.stats.xval_10(node), data_Figure1.tuningcurve.stats.Fval_10(node)-0.03, '($w_{10}$, $S_{10}$)', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(data_Figure1.tuningcurve.stats.xval_90(node), data_Figure1.tuningcurve.stats.Fval_90(node)-0.03, '($w_{90}$, $S_{90}$)', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(1.3, 0.4, {'dynamic range ='; '$10\log_{10}(w_{90}/w_{10})$'}, 'fontsize', fontsize_label, ...
     'fontweight', 'b', 'horizontalalignment', 'center', 'interpreter', 'latex')
set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
    'ylim', [0 1]);
xlabel('recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('regional synaptic gating, $S$', 'fontsize', fontsize_label, 'interpreter', 'latex')

%%% panel letters
annotation(fig, 'textbox', [0.04, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.04, 0.42, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.61, 0.98 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

    
%% FIGURE 2

% load all data relevant to Figure 2
data_Figure2 = load(sprintf('%s/Figure2.mat', data_foldername));

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';
nodeLocations = region_centroids;

ylim_min = min([data_Figure2.tuningcurve.stats{1}.(stats_to_plot); data_Figure2.tuningcurve.stats{2}.(stats_to_plot)]);
ylim_max = max([data_Figure2.tuningcurve.stats{1}.(stats_to_plot); data_Figure2.tuningcurve.stats{2}.(stats_to_plot)]);

fig = figure('Position', [400 100 800 1000]);
for type_ind = 1:length(types)
    type = types{type_ind};
    
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    gain_cmap = repmat(cmap(type_ind,:),N,1);
    brighten_vec = linspace(0.8,0,N);
    for ii=1:N
        gain_cmap(ii,:) = brighten(gain_cmap(ii,:), brighten_vec(ii));
    end

    bw = image_to_plot>0;
    
    [~,xval_transition_ind] = sort(data_Figure2.tuningcurve.stats{type_ind}.xval_transition, 'ascend');
    
    % =====================================================================
    % A: tuning curve
    % =====================================================================
    ax1 = axes('Position', [0.15+0.43*(type_ind-1) 0.78 0.35 0.18]);
    hold on;
    for ii=1:N
        plot(data_Figure2.tuningcurve.w, data_Figure2.tuningcurve.S{type_ind}(xval_transition_ind(ii),:), 'color', gain_cmap(ii,:))
    end
    plot([min(data_Figure2.tuningcurve.stats{type_ind}.xval_transition), max(data_Figure2.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition, S_transition], 'k-', 'linewidth', 2)
    plot([min(data_Figure2.tuningcurve.stats{type_ind}.xval_transition), min(data_Figure2.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    plot([max(data_Figure2.tuningcurve.stats{type_ind}.xval_transition), max(data_Figure2.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    hold off;
    text(max(data_Figure2.tuningcurve.stats{type_ind}.xval_transition)*1.03, S_transition, ...
        ['width = ', sprintf('%.2f', data_Figure2.tuningcurve.stats{type_ind}.dx_transition)])
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel('recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax1_image= axes('Position', [ax1.Position(1)+ax1.Position(3)*0.02 ax1.Position(2)+ax1.Position(4)*0.7 0.04 0.04]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
    
    ax1_1 = axes('Position', [ax1.Position(1)+ax1.Position(3)*0.78 ax1.Position(2)*0.95 0.08 0.18]);
    [ax1_1, ~] = utils.draw_scatterBrain(ax1_1, nodeLocations, data_Figure2.tuningcurve.stats{type_ind}.(stats_to_plot), 20, 'axial');
    colormap(ax1_1, gain_cmap)
    
    % =====================================================================
    % B: dynamic range distribution
    % =====================================================================
    if type_ind==1
        ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)*0.6 ax1.Position(2)-ax1.Position(4)*1.55 ax1.Position(3) ax1.Position(4)]);

        data_violin = utils.padconcatenation(data_Figure2.tuningcurve.stats{1}.(stats_to_plot), data_Figure2.tuningcurve.stats{2}.(stats_to_plot), 2);
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

            text(type_ind_2, max(data_violin(:,type_ind_2))*1.01, ['\sigma = ', sprintf('%.2f', std(data_violin(:,type_ind_2)))], 'color', cmap(type_ind_2,:), ...
            'fontsize', fontsize_label, 'fontweight', 'b', 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
        end

        set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xticklabel', titles)
        ylabel(stats_to_plot_name, 'fontsize', fontsize_label, 'interpreter', 'latex')
        
        for type_ind_2 = 1:length(types)
            if type_ind_2==1
                image_to_plot = human_female;
                yloc = 0.65;
            elseif type_ind_2==2
                image_to_plot = chimpanzee;
                yloc = 0.9;
            end
            bw2 = image_to_plot>0;

            ax2_image = axes('Position', [ax2.Position(1)+0.01+(type_ind_2-1)*0.17 ax2.Position(2)+ax2.Position(4)*yloc 0.04 0.04]);
            imshow(bw2);
            colormap(ax2_image, flipud(gray))
        end
    end
    
    % =====================================================================
    % C: dynamic range spatial arrangement
    % =====================================================================
    dynamic_range_cmap = cbrewer('seq', 'YlGnBu', 130, 'pchip');
    dynamic_range_cmap = flipud(dynamic_range_cmap(end-N+1:end,:));
    
    ax3 = axes('Position', [ax1.Position(1) ax2.Position(2)-ax2.Position(4)*1.2 ax1.Position(3) ax1.Position(4)]);
    ax3_positions = ax3.Position;
    set(ax3, 'visible', 'off')
    
    ax3_1 = axes('Position', [ax3_positions(1)+0.02 ax3_positions(2)-0.05 0.15 0.25]);
    [ax3_1, ~] = utils.draw_scatterBrain(ax3_1, nodeLocations, data_Figure2.tuningcurve.stats{type_ind}.(stats_to_plot), 30, 'axial');
    colormap(ax3_1, dynamic_range_cmap)
    annotation('textbox', [ax3_1.Position(1)-0.03, ax3_1.Position(2)+ax3_1.Position(4)*0.05, 0.1, 0.1], 'string', 'LH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annotation('textbox', [ax3_1.Position(1)+ax3_1.Position(3)-0.07, ax3_1.Position(2)+ax3_1.Position(4)*0.05, 0.1, 0.1], 'string', 'RH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annot = annotation('textarrow', ax3_1.Position(1)*ones(1,2), [ax3_1.Position(2)+ax3_1.Position(4)*0.55 ax3_1.Position(2)+ax3_1.Position(4)*0.4], 'String', 'a', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5);
    annotation('textarrow', [annot.X(1) annot.X(2)], [annot.Y(2) annot.Y(1)], 'String', 'p', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5)
   
    ax3_1_image = axes('Position', [ax3_1.Position(1) ax3_1.Position(2)+ax3_1.Position(4)*0.7 0.04 0.04]);
    imshow(bw);
    colormap(ax3_1_image, flipud(gray))

    ax3_2 = axes('Position', [ax3_1.Position(1)+ax3_1.Position(3)*1.2 ax3_1.Position(2)+ax3_1.Position(4)*0.15 0.10 0.25]);
    [ax3_2, ~] = utils.draw_scatterBrain(ax3_2, nodeLocations, data_Figure2.tuningcurve.stats{type_ind}.(stats_to_plot), 20, 'sagittal_left');
    colormap(ax3_2, dynamic_range_cmap)
    annotation('textbox', [ax3_2.Position(1)-0.04, ax3_2.Position(2)+ax3_2.Position(4)*0.22, 0.1, 0.1], 'string', 'LH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annot = annotation('textarrow', [ax3_2.Position(1)+ax3_2.Position(3)*0.7 ax3_2.Position(1)+ax3_2.Position(3)*0.95], (ax3_2.Position(2)+ax3_2.Position(4)*0.38)*ones(1,2), 'String', 'a', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5);
    annotation('textarrow', [annot.X(2) annot.X(1)], [annot.Y(1) annot.Y(2)], 'String', 'p', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5)
           
    ax3_3 = axes('Position', [ax3_2.Position(1) ax3_2.Position(2)-0.07 ax3_2.Position(3) ax3_2.Position(4)]);
    [ax3_3, ~] = utils.draw_scatterBrain(ax3_3, nodeLocations, data_Figure2.tuningcurve.stats{type_ind}.(stats_to_plot), 20, 'sagittal_right');
    colormap(ax3_3, dynamic_range_cmap)
    annotation('textbox', [ax3_3.Position(1)+ax3_3.Position(3)*0.45, ax3_3.Position(2)+ax3_3.Position(4)*0.22, 0.1, 0.1], 'string', 'RH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annot = annotation('textarrow', [ax3_3.Position(1)+ax3_3.Position(3)*0.3 ax3_3.Position(1)+ax3_3.Position(3)*0.05], (ax3_3.Position(2)+ax3_3.Position(4)*0.38)*ones(1,2), 'String', 'a', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5);
    annotation('textarrow', [annot.X(2) annot.X(1)], [annot.Y(1) annot.Y(2)], 'String', 'p', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'fontweight', 'b', 'headstyle', 'plain', 'headlength', 5, 'headwidth', 5)
           
    cbar = colorbar(ax3_3,'southoutside');
    ylabel(cbar, stats_to_plot_name, 'fontsize', fontsize_label, 'interpreter', 'latex')
    set(cbar, 'fontsize', fontsize_axis, 'ticklength', 0.01, ...
        'position', [ax3_positions(1)+ax3_positions(3)*0.36, ax3_3.Position(2)-0.012, ax3_3.Position(3)*0.9, 0.01], ...
        'ytick', [])
    annotation(fig, 'textbox', [cbar.Position(1)-0.055, cbar.Position(2)*1, 0.04, 0.01], 'string', 'SHARP', 'edgecolor', 'none', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    annotation(fig, 'textbox', [cbar.Position(1)+cbar.Position(3)+0.02, cbar.Position(2)*1, 0.04, 0.01], 'string', 'DIFFUSE', 'edgecolor', 'none', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
           
    % =====================================================================
    % D: dynamic range distribution of different networks
    % =====================================================================
    ax4 = axes('Position', [ax1.Position(1) ax3.Position(2)-ax3.Position(4)*1.4 ax1.Position(3) ax1.Position(4)]);
    data_violin = [];
    for FN = 1:7
        data_violin = cat(2, data_violin, data_Figure2.tuningcurve.stats{type_ind}.FN(:,FN));
    end

    violins = utils.violinplot(data_violin, {});

    for FN = 1:7
        violins(FN).MedianPlot.SizeData = 30;
        violins(FN).ViolinColor = Yeo_7net_colormap(FN,:);
        violins(FN).ViolinAlpha = 0.2;
        violins(FN).ScatterPlot.SizeData = 10;
        violins(FN).ScatterPlot.MarkerFaceAlpha = 1;
        violins(FN).BoxColor = [0 0 0];
        violins(FN).BoxWidth = 0.02;
        violins(FN).WhiskerPlot.LineStyle = 'none';
    end
    set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xticklabel', Yeo_names, ...
        'ylim', [ylim_min*0.95, ylim_max*1.05])
    ylabel(stats_to_plot_name, 'fontsize', fontsize_label, 'interpreter', 'latex')
    
    ax4_image= axes('Position', [ax4.Position(1)+ax4.Position(3)*0.02 ax4.Position(2)+ax4.Position(4)*0.7 0.04 0.04]);
    imshow(bw);
    colormap(ax4_image, flipud(gray))   
end

%%% panel letters
annotation(fig, 'textbox', [0.1, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.3, 0.71, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.1, 0.46, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.1, 0.25, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
    
    
%% FIGURE 3

% load all data relevant to Figure 3
data_Figure3 = load(sprintf('%s/Figure3.mat', data_foldername));

types = {'human', 'chimp', 'macaque', 'marmoset'};
titles = {'human', 'chimpanzee', 'macaque', 'marmoset'};

cmap = [color_brown; color_green; color_blue; color_purple];
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';

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
        image_to_plot = macaque;
    elseif type_ind==4
        image_to_plot = marmoset;
    end
    bw = image_to_plot>0;
    gain_cmap = repmat(cmap(type_ind,:),length(data_Figure3.tuningcurve.stats{type_ind}.max),1);
    brighten_vec = linspace(0.8,0,length(data_Figure3.tuningcurve.stats{type_ind}.max));

    [~,xval_transition_ind] = sort(data_Figure3.tuningcurve.stats{type_ind}.xval_transition, 'ascend');
    
    subplot(2,length(types),type_ind)
    ax1 = gca;
    hold on;
    for ii=1:length(data_Figure3.tuningcurve.stats{type_ind}.max)
        plot(data_Figure3.tuningcurve.w, data_Figure3.tuningcurve.S{type_ind}(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
    end
    plot([min(data_Figure3.tuningcurve.stats{type_ind}.xval_transition), max(data_Figure3.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition, S_transition], 'k-', 'linewidth', 2)
    plot([min(data_Figure3.tuningcurve.stats{type_ind}.xval_transition), min(data_Figure3.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    plot([max(data_Figure3.tuningcurve.stats{type_ind}.xval_transition), max(data_Figure3.tuningcurve.stats{type_ind}.xval_transition)], ...
         [S_transition*0.9, S_transition*1.1], 'k-', 'linewidth', 2)
    hold off;
    text(max(data_Figure3.tuningcurve.stats{type_ind}.xval_transition)*1.05, S_transition, ...
        ['width = ', sprintf('%.2f', data_Figure3.tuningcurve.stats{type_ind}.dx_transition)])
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel('recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax1_image= axes('Position', [ax1.Position(1)-0.01 ax1.Position(2)+ax1.Position(4)*0.795 0.07 0.07]);
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
    data_violin = utils.padconcatenation(data_violin, data_Figure3.tuningcurve.stats{type_ind}.(stats_to_plot)-mean(data_Figure3.tuningcurve.stats{type_ind}.(stats_to_plot)), 2);
end
    
violins = utils.violinplot(data_violin, {});
for type_ind = 1:length(types)
    violins(type_ind).MedianPlot.SizeData = 50;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 10;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.05;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
    text(type_ind-0.4, 1.3, titles{type_ind}, 'color', cmap(type_ind,:), ...
    'fontsize', fontsize_axis, 'fontweight', 'b')
    text(type_ind, 1.3, ['(\sigma = ', sprintf('%.2f)', std(data_Figure3.tuningcurve.stats{type_ind}.(stats_to_plot)-mean(data_Figure3.tuningcurve.stats{type_ind}.(stats_to_plot))))], 'color', cmap(type_ind,:), ...
        'fontsize', fontsize_axis, 'fontweight', 'b', 'horizontalalignment', 'left')
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
        image_to_plot = macaque;
    elseif type_ind==4
        image_to_plot = marmoset;
    end
    bw = image_to_plot>0;
    
    ax2_image= axes('Position', [ax2.Position(1)+ax2.Position(3)*0.78 ax2.Position(2)+ax2.Position(4)*0.8-0.085*(type_ind-1) 0.06 0.06]);
    imshow(bw);
    colormap(ax2_image, flipud(gray))
end

%%% panel letters
annotation(fig, 'textbox', [0.05, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.05, 0.5, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')


%% FIGURE 4

% load all data relevant to Figure 4
data_Figure4 = load(sprintf('%s/Figure4.mat', data_foldername));

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];

w_interest = 0.45;
w_interest_ind = dsearchn(data_Figure4.tuningcurve.w', w_interest);

fig = figure('Position', [200 200 700 600]);
% =========================================================================
% A: time series and autocorrelation fitting method
% =========================================================================
ax1_left = axes('Position', [0.1 0.8 0.35 0.15]);
plot(data_Figure4.timeseries.time, data_Figure4.timeseries.S, 'k', 'linewidth', 1.5)
set(ax1_left, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ytick', []);
xlabel('time (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel({'synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
box off;

ax1_right = axes('Position', [ax1_left.Position(1)+ax1_left.Position(3)+0.15 ax1_left.Position(2) ax1_left.Position(3) ax1_left.Position(4)]);
hold on;
plot(data_Figure4.autocorrelation.lags, data_Figure4.autocorrelation.acf_data, 'k', 'linewidth', 1.5)
plot(data_Figure4.autocorrelation.lags, data_Figure4.autocorrelation.acf_fit, 'r--', 'linewidth', 2)
hold off;
set(ax1_right, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 5], 'ytick', []);
xlabel('time lag, $t$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel('autocorrelation', 'fontsize', fontsize_label, 'interpreter', 'latex')
leg = legend('{\rm data}', 'fit = $c_1e^{-t/\tau} + c_2$', 'location', 'northeast', 'interpreter', 'latex');
set(leg, 'fontsize', fontsize_axis, 'box', 'off', 'numcolumns', 1)
box off;
annotation('arrow', [ax1_left.Position(1)+ax1_left.Position(3)+0.02 ax1_right.Position(1)-0.05], [ax1_left.Position(2)+ax1_left.Position(4)*0.5 ax1_left.Position(2)+ax1_left.Position(4)*0.5])

% =========================================================================
% B: timescales distribution
% =========================================================================
ax2 = axes('Position', [ax1_left.Position(1)+ax1_left.Position(3)*0.7 0.4 0.37 0.28]);

data_violin = utils.padconcatenation(data_Figure4.tuningcurve.stats{1}.tau_neural(:,w_interest_ind), data_Figure4.tuningcurve.stats{2}.tau_neural(:,w_interest_ind), 2);
violins = utils.violinplot(data_violin, {});
for type_ind = 1:length(types)
    violins(type_ind).MedianPlot.SizeData = 20;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 5;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.02;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
    
    text(type_ind, max(data_violin(:,type_ind))*1.02, ['\sigma = ', sprintf('%.2f', std(data_violin(:,type_ind)))], 'color', cmap(type_ind,:), ...
    'fontsize', fontsize_label, 'fontweight', 'b', 'horizontalalignment', 'center', 'verticalalignment', 'bottom')
end

set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.01, 0.01], 'xticklabel', titles)
ylabel('timescale, $\tau$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')

for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
        yloc = 0.25;
    elseif type_ind==2
        image_to_plot = chimpanzee;
        yloc = 0.9;
    end
    bw = image_to_plot>0;
    
    ax2_image = axes('Position', [ax2.Position(1)+(type_ind-1)*0.17 ax2.Position(2)+ax2.Position(4)*yloc 0.05 0.05]);
    imshow(bw);
    colormap(ax2_image, flipud(gray))
end

% =========================================================================
% C: correlation of timescales and dynamic range
% =========================================================================
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;

    ax3 = axes('Position', [ax1_left.Position(1)+(ax1_left.Position(3)+0.13)*(type_ind-1) 0.08 0.37 0.25]);
    data_to_plot_x = zscore(data_Figure4.tuningcurve.stats{type_ind}.dynamic_range);
    data_to_plot_y = data_Figure4.tuningcurve.stats{type_ind}.tau_neural(:,w_interest_ind);
    plot(data_to_plot_x, data_to_plot_y, '.', 'color', cmap(type_ind,:), 'markersize', 12)
    set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02])
    xlabel('dynamic range (z score)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    ylabel('timescale, $\tau$ (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'spearman');
    text(0, min(get(ax3, 'ylim'))+diff(get(ax3, 'ylim'))*0.95, ['\rho ', sprintf('= %.2g', rho)], ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'center');
    text(0, min(get(ax3, 'ylim'))+diff(get(ax3, 'ylim'))*0.85, sprintf('p = %.2g', pval), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'center');

    ax3_image= axes('Position', [ax3.Position(1)+0.02 ax3.Position(2)+ax3.Position(4)*0.7 0.05 0.05]);
    imshow(bw);
    colormap(ax3_image, flipud(gray))
    
    if type_ind==2
        ax3_inset = axes('Position', [ax3.Position(1)+0.07 ax3.Position(2)+0.08 ax3.Position(3)*0.25 ax3.Position(4)*0.25]);
        plot(data_to_plot_x, data_to_plot_y, '.', 'color', cmap(type_ind,:), 'markersize', 10)
        set(ax3_inset, 'fontsize', fontsize_axis-4, 'ticklength', [0.02, 0.02], ...
            'ytick', [0.1, 0.15, 0.2, 0.25], 'ylim', [0.1 0.25], 'box', 'on')
        xlabel('dynamic range', 'fontsize', fontsize_label-4, 'interpreter', 'latex')
        ylabel('$\tau$', 'fontsize', fontsize_label-4, 'interpreter', 'latex')
    end
end

%%% panel letters
annotation(fig, 'textbox', [0.03, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.27, 0.71, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.03, 0.36, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')


%% FIGURE 5

% load all data relevant to Figure 5
data_Figure5 = load(sprintf('%s/Figure5.mat', data_foldername));

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';
thres = 1;

edge_alpha = 0.8;
edge_width = 1;
           
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
nodeLocations = region_centroids;
nodeSizes = (markersize/10)*ones(N,1);
nodeSizes(node_interests) = markersize;

fig = figure('Position', [200 200 700 600]);
% =========================================================================
% A left: connectome
% =========================================================================
W = connectome_human;
ax1 = axes('Position', [0.02 0.63 0.3 0.3]);

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
text(-37, -45, 'LH', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, ...
    'fontweight', 'b')
text(37, -45, 'RH', 'horizontalalignment', 'center', 'fontsize', fontsize_axis, ...
    'fontweight', 'b')
colormap(ax1, node_cmap)
view(ax1, 2)
axis(ax1, 'equal')
set(ax1, 'visible', 'off')

% =========================================================================
% B left: model schematic
% =========================================================================
ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)*1.02 ax1.Position(2) 0.4 ax1.Position(4)]);
circle_locs_x = [0.25, 0.57];
circle_locs_y = [0.5, 0.5];
hold on;
plot(circle_locs_x, circle_locs_y, 'k-', 'linewidth', 3)
plot(circle_locs_x(1), circle_locs_y(1), 'ko', 'markersize', 45, 'markerfacecolor', node_interests_colors(1,:))
plot(circle_locs_x(2), circle_locs_y(2), 'ko', 'markersize', 45, 'markerfacecolor', node_interests_colors(2,:))
hold off;
set(ax2, 'xtick', [], 'ytick', [], 'xlim', [0 1], 'ylim', [0 1])

% annotation lines
utils.annotation_pinned('arrow', [circle_locs_x(1) circle_locs_x(1)], [circle_locs_y(1)*1.7 circle_locs_y(1)*1.25], 'linewidth', 3, 'headlength', 10, 'headwidth', 18);
utils.annotation_pinned('arrow', [circle_locs_x(2) circle_locs_x(2)], [circle_locs_y(1)*1.7 circle_locs_y(1)*1.25], 'linewidth', 3, 'headlength', 10, 'headwidth', 18);
utils.annotation_pinned('arrow', [circle_locs_x(1) circle_locs_x(1)], [circle_locs_y(1)*0.3 circle_locs_y(1)*0.75], 'linewidth', 3, 'headlength', 10, 'headwidth', 18);
utils.annotation_pinned('arrow', [circle_locs_x(2) circle_locs_x(2)], [circle_locs_y(1)*0.3 circle_locs_y(1)*0.75], 'linewidth', 3, 'headlength', 10, 'headwidth', 18);

% texts
text(circle_locs_x(1), circle_locs_y(1), {'region'; '$i$'}, 'horizontalalignment', 'center', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(circle_locs_x(2), circle_locs_y(2), {'region'; '$j$'}, 'horizontalalignment', 'center', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(circle_locs_x(1), circle_locs_y(1)+0.41, '$\beta$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(2), circle_locs_y(1)+0.41, '$\beta$', 'horizontalalignment', 'center', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(mean(circle_locs_x), circle_locs_y(1), '$L_{ij}$', 'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(1), circle_locs_y(1)-0.36, '$D$', 'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
text(circle_locs_x(2), circle_locs_y(1)-0.36, '$D$', 'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontsize', fontsize_label*1.2, 'interpreter', 'latex')
axis off

% =========================================================================
% B top right: time series
% =========================================================================
ax3 = axes('Position', [ax2.Position(1)+ax2.Position(3)*0.95 ax2.Position(2)+ax2.Position(4)*0.6 0.2 0.1]);
hold on;
plot(data_Figure5.timeseries.time, data_Figure5.timeseries.y(node_interests(1),:), 'color', node_interests_colors(1,:), 'linewidth', 1.5)
plot(data_Figure5.timeseries.time, data_Figure5.timeseries.y(node_interests(2),:), 'color', node_interests_colors(2,:), 'linewidth', 1.5)
plot(data_Figure5.timeseries.time, thres*ones(size(data_Figure5.timeseries.time)), 'k-', 'linewidth', 2);
plot(data_Figure5.timeseries.time, -thres*ones(size(data_Figure5.timeseries.time)), 'k-', 'linewidth', 2);
hold off;
text(4, thres, 'correct', ...
    'horizontalalignment', 'left', 'verticalalignment', 'middle','fontsize', fontsize_axis, 'fontweight', 'b', 'color', 'b')
text(4, -thres, 'incorrect', ...
    'horizontalalignment', 'left', 'verticalalignment', 'middle','fontsize', fontsize_axis, 'fontweight', 'b', 'color', 'r')
set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 4], 'ylim', [-thres,thres], 'ytick', []);
xlabel('time (s)', 'fontsize', fontsize_axis, 'interpreter', 'latex')
ylabel({'regional'; 'decision evidence'}, 'fontsize', fontsize_axis, 'interpreter', 'latex')

% =========================================================================
% B top bottom: regional accuracy
% =========================================================================
ax4 = axes('Position', [ax3.Position(1) ax3.Position(2)-0.18 ax3.Position(3) ax3.Position(4)]);
hold on;
plot(data_Figure5.timeseries.time, data_Figure5.timeseries.accuracy(node_interests(1),:), 'color', node_interests_colors(1,:), 'linewidth', 1.5)
plot(data_Figure5.timeseries.time, data_Figure5.timeseries.accuracy(node_interests(2),:), 'color', node_interests_colors(2,:), 'linewidth', 1.5)
hold off;
set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 1.5]);
xlabel('time (s)', 'fontsize', fontsize_axis, 'interpreter', 'latex')
ylabel({'regional'; 'accuracy ($\%$)'}, 'fontsize', fontsize_axis, 'interpreter', 'latex')

% =========================================================================
% C: human - chimpanzee whole-brain accuracy
% =========================================================================
ax5 = axes('Position', [0.13 0.1 0.35 0.4]);
data_to_plot_x = data_Figure5.decision.time;
data_to_plot_y = mean(data_Figure5.decision.accuracy{1},1)-mean(data_Figure5.decision.accuracy{2},1);
[min_diff, min_diff_ind] = min(data_to_plot_y);
hold on;
plot(data_to_plot_x, data_to_plot_y, 'k-', 'linewidth', 2)
plot(data_to_plot_x, 0*(ones(size(data_to_plot_x))), 'k:');
plot(data_to_plot_x(min_diff_ind)*ones(1,2), get(gca,'ylim'), 'k--');
hold off;
text(data_to_plot_x(min_diff_ind), min(get(gca,'ylim')), '$t_{\rm min}$', 'interpreter', 'latex', ...
     'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontsize', fontsize_label)
set(ax5, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 1.5])
xlabel('time (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel({'human -- chimpanzee'; 'whole-brain accuracy ($\%$)'}, 'fontsize', fontsize_label, 'interpreter', 'latex')

% =========================================================================
% D: correlation of accuracy at tmin with dynamic range
% =========================================================================
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;

    ax6 = axes('Position', [ax5.Position(1)+ax5.Position(3)+0.15 ax5.Position(2)+(ax5.Position(4)*0.6)*(2-type_ind) ax5.Position(3)*0.9 ax5.Position(4)*0.4]);
    data_to_plot_x = zscore(data_Figure5.tuningcurve.stats{type_ind}.dynamic_range);
    data_to_plot_y = squeeze(data_Figure5.decision.accuracy{type_ind}(:,min_diff_ind));
    hold on;
    plot(data_to_plot_x, data_to_plot_y, '.', 'color', cmap(type_ind,:), 'markersize', 12)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'k-', 'linewidth', 2);
    hold off;
    set(ax6, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [10 50])
    if type_ind==2
        xlabel('dynamic range (z score)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    if type_ind==2
        text(-2.6, min(data_to_plot_y)+5, {'regional accuracy at $t_{\rm min}$ ($\%$)'}, 'fontsize', fontsize_label, 'rotation', 90, 'interpreter', 'latex')
    end
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    text(max(data_to_plot_x), max(get(gca, 'ylim'))*0.9, sprintf('r = %.2g', rho), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
    text(max(data_to_plot_x), max(get(gca, 'ylim'))*0.78, sprintf('p = %.2g', pval), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');

    ax6_image = axes('Position', [ax6.Position(1)+0.01 ax6.Position(2)+ax6.Position(4)*0.75 0.06 0.06]);
    imshow(bw);
    colormap(ax6_image, flipud(gray))
end

%%% panel letters
annotation(fig, 'textbox', [0.06, 0.97, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.38, 0.97, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.06, 0.55, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.57, 0.55, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
    