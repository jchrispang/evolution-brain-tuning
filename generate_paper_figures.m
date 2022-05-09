% generate_paper_figures.m
% 
% Matlab code to generate the main figures of the paper

%% LOAD CONNECTOME MATRICES, REGION LOCATIONS/NAMES, CORTICAL SURFACES, PARCELLATION, MULTIPLE COMPARISONS P-VALUE

load('data/connectome_human.mat')
load('data/connectome_chimp.mat')
load('data/connectome_macaque.mat')
load('data/connectome_marmoset.mat')
load('data/region_props_chimp.mat', 'region_centroids', 'region_names') 
load('data/RSN.mat', 'RSN_indices', 'RSN_names')
load('data/macaque_networks.mat', 'networks_indices', 'networks_names')
load('data/data_figures/multiple_comparisons.mat', 'multiple_comparisons', 'pvals_adj')

human_surface = @(wb_hemi, surface_interest) sprintf('data/S1200.%s.%s_MSMAll.32k_fs_LR.surf.gii', wb_hemi, surface_interest);
chimp_surface = @(wb_hemi, surface_interest) sprintf('data/ChimpYerkes29_v1.2.%s.%s.32k_fs_LR.surf.gii', wb_hemi, surface_interest);

surface_interest = 'very_inflated';
surface_file{1}.lh = human_surface('L', surface_interest);
surface_file{1}.rh = human_surface('R', surface_interest);
surface_file{2}.lh = chimp_surface('L', surface_interest);
surface_file{2}.rh = chimp_surface('R', surface_interest);

parc_filename = @(hemisphere) sprintf('data/fsLR_32k_DK114-%s.txt', hemisphere);

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

show_padj = 1;

data_foldername = 'data/data_figures';

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
plot(data_Figure1.tuningcurve.stats.xval_10(node), data_Figure1.tuningcurve.stats.Fval_10(node), 'k.', 'markersize', 20)
plot(data_Figure1.tuningcurve.stats.xval_90(node), data_Figure1.tuningcurve.stats.Fval_90(node), 'k.', 'markersize', 20)
hold off;
text(data_Figure1.tuningcurve.stats.xval_10(node), data_Figure1.tuningcurve.stats.Fval_10(node)-0.04, '($w_{10}$, $S_{10}$)', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(data_Figure1.tuningcurve.stats.xval_90(node), data_Figure1.tuningcurve.stats.Fval_90(node)-0.04, '($w_{90}$, $S_{90}$)', 'fontsize', fontsize_label, 'interpreter', 'latex')
text(1.3, 0.4, {'dynamic range ='; '$10\log_{10}(w_{90}/w_{10})$'}, 'fontsize', fontsize_label, ...
     'fontweight', 'b', 'horizontalalignment', 'center', 'interpreter', 'latex')
set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
    'ylim', [0 1]);
xlabel('global recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
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

N = 114;
cmap = [color_brown; color_green];
S_transition = 0.3;
stats_to_plot = 'dynamic_range';
stats_to_plot_name = 'dynamic range';

parc_lh = dlmread(parc_filename('lh'));
parcels_lh = unique(parc_lh(parc_lh>0));
num_parcels_lh = length(parcels_lh);
parc_rh = dlmread(parc_filename('rh'));
parcels_rh = unique(parc_rh(parc_rh>0));
num_parcels_rh = length(parcels_rh);

ylim_min = min([data_Figure2.tuningcurve.stats{1}.(stats_to_plot); data_Figure2.tuningcurve.stats{2}.(stats_to_plot)]);
ylim_max = max([data_Figure2.tuningcurve.stats{1}.(stats_to_plot); data_Figure2.tuningcurve.stats{2}.(stats_to_plot)]);
image_width = 0.035;

fig = figure('Position', [400 100 700 1000]);
for type_ind = 1:length(types)
    type = types{type_ind};
    
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    gain_cmap = repmat(cmap(type_ind,:),N,1);
%     brighten_vec = linspace(0.8,0,N);
    brighten_vec = zeros(N,1);
    for ii=1:N
        gain_cmap(ii,:) = brighten(gain_cmap(ii,:), brighten_vec(ii));
    end

    bw = image_to_plot>0;
    
    [~,xval_transition_ind] = sort(data_Figure2.tuningcurve.stats{type_ind}.xval_transition, 'ascend');
    
    % =====================================================================
    % A: tuning curve
    % =====================================================================
    ax1 = axes('Position', [0.15+0.43*(type_ind-1) 0.84 0.35 0.13]);
    hold on;
    for ii=1:N
        plot(data_Figure2.tuningcurve.w, data_Figure2.tuningcurve.S{type_ind}(xval_transition_ind(ii),:), 'color', gain_cmap(ii,:))
    end
    hold off;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel('global recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax1_image= axes('Position', [ax1.Position(1)+ax1.Position(3)*0.02 ax1.Position(2)+ax1.Position(4)*0.7 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
    
    % =====================================================================
    % B: dynamic range distribution
    % =====================================================================
    if type_ind==1
        ax2 = axes('Position', [ax1.Position(1)+ax1.Position(3)*0.6 ax1.Position(2)-ax1.Position(4)*1.65 ax1.Position(3) ax1.Position(4)]);

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

            ax2_image = axes('Position', [ax2.Position(1)+0.01+(type_ind_2-1)*0.17 ax2.Position(2)+ax2.Position(4)*yloc image_width image_width]);
            imshow(bw2);
            colormap(ax2_image, flipud(gray))
        end
    end
    
    % =====================================================================
    % C: dynamic range spatial arrangement
    % =====================================================================
    dynamic_range_cmap = cbrewer('seq', 'YlGnBu', 130, 'pchip');
    dynamic_range_cmap = [0.5 0.5 0.5; flipud(dynamic_range_cmap(end-N+1:end,:))];
    
    surface_to_plot_lh = gifti(surface_file{type_ind}.lh);
    surface_to_plot_rh = gifti(surface_file{type_ind}.rh);
    
    boundary_method = 'midpoint';
    BOUNDARY_lh = findROIboundaries(surface_to_plot_lh.vertices,surface_to_plot_lh.faces,parc_lh,boundary_method);
    BOUNDARY_rh = findROIboundaries(surface_to_plot_rh.vertices,surface_to_plot_rh.faces,parc_rh,boundary_method);
    
    data_parcel = data_Figure2.tuningcurve.stats{type_ind}.(stats_to_plot);
    
    data_to_plot_lh = zeros(size(surface_to_plot_lh.vertices,1),1);
    data_to_plot_rh = zeros(size(surface_to_plot_rh.vertices,1),1);
    for ii=1:num_parcels_lh
        parcel = parcels_lh(ii);
        
        data_to_plot_lh(find(parc_lh==parcel)) = data_parcel(ii);
    end
    for ii=1:num_parcels_rh
        parcel = parcels_rh(ii);
        
        data_to_plot_rh(find(parc_rh==parcel)) = data_parcel(num_parcels_lh+ii);
    end
    
%     clims = [min(data_parcel(data_parcel>0)), max(data_parcel(data_parcel>0))];
    clims = prctile(data_parcel(data_parcel>0), [4 96]);
    
    ind_zeros = find(parc_lh==0);
    data_to_plot_lh(data_to_plot_lh<clims(1)) = clims(1);
    data_to_plot_lh(data_to_plot_lh>clims(2)) = clims(2);
    data_to_plot_lh(ind_zeros) = clims(1)-0.1;
    ind_zeros = find(parc_rh==0);
    data_to_plot_rh(data_to_plot_rh<clims(1)) = clims(1);
    data_to_plot_rh(data_to_plot_rh>clims(2)) = clims(2);
    data_to_plot_rh(ind_zeros) = clims(1)-0.1;
    clims(1) = clims(1)-0.1;
    
    ax3 = axes('Position', [ax1.Position(1) ax2.Position(2)-ax2.Position(4)*1.4 ax1.Position(3) ax1.Position(4)]);
    ax3_positions = ax3.Position;
    set(ax3, 'visible', 'off')
    
    ax3_1 = axes('Position', [ax3_positions(1)+0.01 ax3_positions(2)-0.03 0.16 0.26]);
    patch(ax3_1, 'Vertices', surface_to_plot_lh.vertices, 'Faces', surface_to_plot_lh.faces, 'FaceVertexCData', data_to_plot_lh, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for ii = 1:length(BOUNDARY_lh)
            plot3(BOUNDARY_lh{ii}(:,1), BOUNDARY_lh{ii}(:,2), BOUNDARY_lh{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end
    hold off;
    view([-90 0]);
    caxis(clims)
    material dull
    camlight('headlight');
    colormap(ax3_1, dynamic_range_cmap)
    axis off
    axis image
    
    ax3_1_image = axes('Position', [ax3_positions(1)-0.03 ax3_1.Position(2)+ax3_1.Position(4)*0.55 image_width image_width]);
    imshow(bw);
    colormap(ax3_1_image, flipud(gray))
    
    ax3_2 = axes('Position', [ax3_1.Position(1) ax3_1.Position(2)-ax3_1.Position(4)/3.1 ax3_1.Position(3) ax3_1.Position(4)]);
    patch(ax3_2, 'Vertices', surface_to_plot_lh.vertices, 'Faces', surface_to_plot_lh.faces, 'FaceVertexCData', data_to_plot_lh, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for ii = 1:length(BOUNDARY_lh)
            plot3(BOUNDARY_lh{ii}(:,1), BOUNDARY_lh{ii}(:,2), BOUNDARY_lh{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end
    hold off;
    view([90 0]);
    caxis(clims)
    material dull
    camlight('headlight');
    colormap(ax3_2, dynamic_range_cmap)
    axis off
    axis image
    
    ax3_3 = axes('Position', [ax3_1.Position(1)+ax3_1.Position(3)*1.1 ax3_1.Position(2) ax3_1.Position(3) ax3_1.Position(4)]);
    patch(ax3_3, 'Vertices', surface_to_plot_rh.vertices, 'Faces', surface_to_plot_rh.faces, 'FaceVertexCData', data_to_plot_rh, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for ii = 1:length(BOUNDARY_rh)
            plot3(BOUNDARY_rh{ii}(:,1), BOUNDARY_rh{ii}(:,2), BOUNDARY_rh{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end
    hold off;
    view([90 0]);
    caxis(clims)
    material dull
    camlight('headlight');
    colormap(ax3_3, dynamic_range_cmap)
    axis off
    axis image
    
    ax3_4 = axes('Position', [ax3_3.Position(1) ax3_2.Position(2) ax3_1.Position(3) ax3_1.Position(4)]);
    patch(ax3_4, 'Vertices', surface_to_plot_rh.vertices, 'Faces', surface_to_plot_rh.faces, 'FaceVertexCData', data_to_plot_rh, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for ii = 1:length(BOUNDARY_rh)
            plot3(BOUNDARY_rh{ii}(:,1), BOUNDARY_rh{ii}(:,2), BOUNDARY_rh{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end
    hold off;
    view([-90 0]);
    caxis(clims)
    material dull
    camlight('headlight');
    colormap(ax3_4, dynamic_range_cmap)
    axis off
    axis image
    
    cbar = colorbar(ax3_4,'southoutside');
    ylabel(cbar, stats_to_plot_name, 'fontsize', fontsize_label, 'interpreter', 'latex')
    set(cbar, 'fontsize', fontsize_axis, 'ticklength', 0.01, ...
        'position', [ax3_positions(1)+ax3_positions(3)*0.32, ax3_3.Position(2)-0.012, ax3_3.Position(3)*0.8, 0.01], ...
        'ytick', [])
    annotation(fig, 'textbox', [cbar.Position(1)-0.05, cbar.Position(2)*1, 0.04, 0.01], 'string', 'sharp', 'edgecolor', 'none', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    annotation(fig, 'textbox', [cbar.Position(1)+cbar.Position(3)+0.015, cbar.Position(2)*1, 0.04, 0.01], 'string', 'diffuse', 'edgecolor', 'none', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')    
   
    annotation('textbox', [ax3_1.Position(1)-0.04, ax3_1.Position(2)+0.038, 0.1, 0.1], 'string', 'LH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annotation('textbox', [ax3_1.Position(1)+ax3_1.Position(3)*1.73, ax3_1.Position(2)+0.038, 0.1, 0.1], 'string', 'RH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    
    % =====================================================================
    % D: dynamic range distribution of different networks
    % =====================================================================
    ax4 = axes('Position', [ax1.Position(1) ax3.Position(2)-ax3.Position(4)*1.65 ax1.Position(3) ax1.Position(4)]);
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
    
    ax4_image= axes('Position', [ax4.Position(1)+ax4.Position(3)*0.02 ax4.Position(2)+ax4.Position(4)*0.7 image_width image_width]);
    imshow(bw);
    colormap(ax4_image, flipud(gray))
    
    % =====================================================================
    % E: FC of different networks
    % =====================================================================
    ax5 = axes('Position', [ax1.Position(1) ax4.Position(2)-ax4.Position(4)*1.35 ax1.Position(3) ax1.Position(4)]);
    hold on;
    for FN=1:7
        plot(data_Figure2.tuningcurve.w, squeeze(data_Figure2.between_FC{type_ind}(FN,FN,:)), '-', 'color', Yeo_7net_colormap(FN,:), 'linewidth', 1)
    end
    plot(data_Figure2.tuningcurve.w, data_Figure2.global_FC{type_ind}, 'k-', 'linewidth', 2)
    hold off;
    if type_ind==2
        leg = legend([Yeo_names, 'whole brain'], 'fontsize', fontsize_axis-2, 'location', 'northeast', 'box', 'off');
        leg.ItemTokenSize = leg.ItemTokenSize/2;
    end
    set(ax5, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 0.3])
    xlabel('global recurrent strength, $w$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel('$FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    
    ax5_image= axes('Position', [ax5.Position(1)+ax5.Position(3)*0.02 ax5.Position(2)+ax5.Position(4)*0.7 image_width image_width]);
    imshow(bw);
    colormap(ax5_image, flipud(gray))
end

%%% panel letters
annotation(fig, 'textbox', [0.09, 0.99, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.3, 0.77, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.09, 0.58, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.09, 0.39, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.09, 0.20, 0.01, 0.01], 'string', 'E', 'edgecolor', 'none', ...
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
image_width = 0.065;

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
%     brighten_vec = linspace(0.8,0,length(data_Figure3.tuningcurve.stats{type_ind}.max));
    brighten_vec = zeros(length(data_Figure3.tuningcurve.stats{type_ind}.max),1);

    [~,xval_transition_ind] = sort(data_Figure3.tuningcurve.stats{type_ind}.xval_transition, 'ascend');
    
    subplot(2,length(types),type_ind)
    ax1 = gca;
    hold on;
    for ii=1:length(data_Figure3.tuningcurve.stats{type_ind}.max)
        plot(data_Figure3.tuningcurve.w, data_Figure3.tuningcurve.S{type_ind}(xval_transition_ind(ii),:), 'color', brighten(gain_cmap(ii,:),brighten_vec(ii)))
    end
    hold off;
    set(gca, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [0 1], ...
                'xtick', 0:0.5:2, 'xlim', [0 2]);
    xlabel({'global recurrent'; 'strength, $w$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel({'regional synaptic'; 'gating, $S$'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    
    ax1_image= axes('Position', [ax1.Position(1)-0.01 ax1.Position(2)+ax1.Position(4)*0.81 image_width image_width]);
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
    violins(type_ind).BoxWidth = 0.03;
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

w_interest = 0.45;
w_interest_ind = dsearchn(data_Figure4.tuningcurve.w', w_interest);

pvals_adj_ind = find(strcmpi(multiple_comparisons(:,1), 'Figure4'));

image_width = 0.04;

fig = figure('Position', [200 200 700 900]);
% =========================================================================
% A: correlation of timescales and dynamic range
% =========================================================================
ax0 = axes('Position', [0.1 0.71 0.4 0.23]);
ax0_position = ax0.Position;
set(ax0,'visible','off')

for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    bw = image_to_plot>0;

    ax1 = axes('Position', [ax0_position(1)+(ax0_position(3)+0.08)*(type_ind-1) ax0_position(2) ax0_position(3) ax0_position(4)]);
    data_to_plot_x = tiedrank(zscore(data_Figure4.tuningcurve.stats{type_ind}.dynamic_range));
    data_to_plot_y = tiedrank(data_Figure4.tuningcurve.stats{type_ind}.tau_neural(:,w_interest_ind));
    hold on;
    plot(data_to_plot_x, data_to_plot_y, '.', 'color', cmap(type_ind,:), 'markersize', 20)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'k-', 'linewidth', 2);
    hold off;
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'box', 'off', 'xlim', [1 120], 'ylim', [1 120])
    xlabel('rank of dynamic range', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel('rank of timescale', 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    text(max(get(ax1,'xlim'))+5, 34, sprintf('r = %.2g', rho), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
    if ~show_padj
        text(max(get(ax1,'xlim')), 25, utils.extract_pvalue_text(pval), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
    else
        text(max(get(ax1,'xlim'))+5, 25, utils.extract_pvalue_text(pvals_adj(pvals_adj_ind(0+type_ind)), 1), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
 
    ax1_image = axes('Position', [ax1.Position(1)+0.01 ax1.Position(2)+ax1.Position(4)*0.85 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
end

% =========================================================================
% B left: brain network
% =========================================================================
ax2 = axes('Position', [0.0 0.34 0.3 ax0_position(4)]);
 
W = connectome_human;
[edge_X, edge_Y, edge_Z] = adjacency_plot_und(threshold_proportional(W, frac), nodeLocations);  % get all the edges
[nodeLocations, edges] = utils.extract_scatterBrain_locs_edges(nodeLocations, edge_X, edge_Y, edge_Z, slice);
 
hold on;
if strcmpi(slice, 'axial')
    scatter3(ax2, nodeLocations(:,1), nodeLocations(:,2), max(nodeLocations(:,3))*ones(N,1), ...
             nodeSizes, node_cval, 'filled', 'markeredgecolor', 'k');
else
    scatter3(ax2, nodeLocations(:,1), nodeLocations(:,2), nodeLocations(:,3), ...
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
colormap(ax2, node_cmap)
view(ax2, 2)
axis(ax2, 'equal')
set(ax2, 'visible', 'off')

% =========================================================================
% B right: drift diffusion model schematic
% =========================================================================
ax3 = axes('Position', [ax2.Position(1)+ax2.Position(3)*1.0 ax2.Position(2) 0.4 ax2.Position(4)]);
circle_locs_x = [0.25, 0.57];
circle_locs_y = [0.5, 0.5];
hold on;
plot(circle_locs_x, circle_locs_y, 'k-', 'linewidth', 3)
plot(circle_locs_x(1), circle_locs_y(1), 'ko', 'markersize', 45, 'markerfacecolor', node_interests_colors(1,:))
plot(circle_locs_x(2), circle_locs_y(2), 'ko', 'markersize', 45, 'markerfacecolor', node_interests_colors(2,:))
hold off;
set(ax3, 'xtick', [], 'ytick', [], 'xlim', [0 1], 'ylim', [0 1])
 
% annotation lines
annotation('arrow', [ax3.Position(1)+ax3.Position(3)*circle_locs_x(1) ax3.Position(1)+ax3.Position(3)*circle_locs_x(1)], [ax3.Position(2)+ax3.Position(4)*0.84 ax3.Position(2)+ax3.Position(4)*0.65], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 18)
annotation('arrow', [ax3.Position(1)+ax3.Position(3)*circle_locs_x(2) ax3.Position(1)+ax3.Position(3)*circle_locs_x(2)], [ax3.Position(2)+ax3.Position(4)*0.84 ax3.Position(2)+ax3.Position(4)*0.65], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 18)
annotation('arrow', [ax3.Position(1)+ax3.Position(3)*circle_locs_x(1) ax3.Position(1)+ax3.Position(3)*circle_locs_x(1)], [ax3.Position(2)+ax3.Position(4)*0.16 ax3.Position(2)+ax3.Position(4)*0.35], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 18)
annotation('arrow', [ax3.Position(1)+ax3.Position(3)*circle_locs_x(2) ax3.Position(1)+ax3.Position(3)*circle_locs_x(2)], [ax3.Position(2)+ax3.Position(4)*0.16 ax3.Position(2)+ax3.Position(4)*0.35], ...
    'linewidth', 3, 'headlength', 10, 'headwidth', 18)
 
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
% C: time series
% =========================================================================
ax4 = axes('Position', [ax3.Position(1)+ax3.Position(3)*1.02 ax3.Position(2)+ax3.Position(4)*0.65 ax3.Position(3)*0.48 ax3.Position(4)*0.3]);
hold on;
plot(data_Figure4.timeseries.time, data_Figure4.timeseries.y(node_interests(1),:), 'color', node_interests_colors(1,:), 'linewidth', 1.5)
plot(data_Figure4.timeseries.time, data_Figure4.timeseries.y(node_interests(2),:), 'color', node_interests_colors(2,:), 'linewidth', 1.5)
yline(thres, 'k-', 'linewidth', 2);
yline(-thres, 'k-', 'linewidth', 2);
hold off;
text(4, thres, 'correct', ...
    'horizontalalignment', 'left', 'verticalalignment', 'middle','fontsize', fontsize_axis, 'fontweight', 'b', 'color', 'b')
text(4, -thres, 'incorrect', ...
    'horizontalalignment', 'left', 'verticalalignment', 'middle','fontsize', fontsize_axis, 'fontweight', 'b', 'color', 'r')
set(ax4, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 4], 'ylim', [-thres,thres], 'ytick', []);
xlabel('time (s)', 'fontsize', fontsize_axis, 'interpreter', 'latex')
ylabel({'regional'; 'decision evidence'}, 'fontsize', fontsize_axis, 'interpreter', 'latex')

% =========================================================================
% D: accuracy vs time
% =========================================================================
ax5 = axes('Position', [ax4.Position(1) ax4.Position(2)-ax4.Position(4)*1.9 ax4.Position(3) ax4.Position(4)]);
hold on;
plot(data_Figure4.timeseries.time, data_Figure4.timeseries.accuracy(node_interests(1),:), 'color', node_interests_colors(1,:), 'linewidth', 1.5)
plot(data_Figure4.timeseries.time, data_Figure4.timeseries.accuracy(node_interests(2),:), 'color', node_interests_colors(2,:), 'linewidth', 1.5)
hold off;
set(ax5, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 1.5]);
xlabel('time (s)', 'fontsize', fontsize_axis, 'interpreter', 'latex')
ylabel({'regional'; 'accuracy ($\%$)'}, 'fontsize', fontsize_axis, 'interpreter', 'latex')

% =========================================================================
% E: human - chimpanzee whole-brain accuracy
% =========================================================================
ax6 = axes('Position', [ax0_position(1) 0.05 ax0_position(3) ax0_position(4)]);
data_to_plot_x = data_Figure4.decision.time;
data_to_plot_y = mean(data_Figure4.decision.accuracy{1},1)-mean(data_Figure4.decision.accuracy{2},1);
[~, min_diff_ind] = min(data_to_plot_y);
hold on;
plot(data_to_plot_x, data_to_plot_y, 'k-', 'linewidth', 2)
plot(data_to_plot_x, 0*(ones(size(data_to_plot_x))), 'k:');
plot(data_to_plot_x(min_diff_ind)*ones(1,2), get(gca,'ylim'), 'k--');
hold off;
text(data_to_plot_x(min_diff_ind), min(get(gca,'ylim')), '$t_{\rm min}$', 'interpreter', 'latex', ...
     'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontsize', fontsize_label)
set(ax6, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'xlim', [0 1.5])
xlabel('time (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
ylabel({'human -- chimpanzee'; 'whole-brain accuracy ($\%$)'}, 'fontsize', fontsize_label, 'interpreter', 'latex')

% =========================================================================
% F: correlation of accuracy at tmin with dynamic range
% =========================================================================
for type_ind=1:length(types)
    if type_ind==1
        image_to_plot = human_female;
        yloc = ax6.Position(2)+ax6.Position(4)*0.6;
    elseif type_ind==2
        image_to_plot = chimpanzee;
        yloc = ax6.Position(2);
    end
    bw = image_to_plot>0;
    
    ax7 = axes('Position', [ax6.Position(1)+ax6.Position(3)+0.12 yloc ax6.Position(3)*0.9 ax6.Position(4)*0.4]);
    data_to_plot_x = zscore(data_Figure4.tuningcurve.stats{type_ind}.dynamic_range);
    data_to_plot_y = squeeze(data_Figure4.decision.accuracy{type_ind}(:,min_diff_ind));
    hold on;
    plot(data_to_plot_x, data_to_plot_y, '.', 'color', cmap(type_ind,:), 'markersize', 12)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'k-', 'linewidth', 2);
    hold off;
    set(ax7, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'ylim', [10 50], 'box', 'off')
    if type_ind==2
        xlabel('dynamic range (z score)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    if type_ind==2
        text(-2.55, min(data_to_plot_y)+7, {'regional accuracy at $t_{\rm min}$ ($\%$)'}, 'fontsize', fontsize_label, 'rotation', 90, 'interpreter', 'latex')
    end
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    text(max(get(ax7, 'xlim'))*0.99, max(get(ax7, 'ylim'))*0.9, sprintf('r = %.2f', rho), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
    if ~show_padj
        text(max(get(ax7, 'xlim'))*0.99, max(get(ax7, 'ylim'))*0.75, utils.extract_pvalue_text(pval), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
    else
        text(max(get(ax7, 'xlim'))*0.99, max(get(ax7, 'ylim'))*0.71, utils.extract_pvalue_text(pvals_adj(pvals_adj_ind(2+type_ind)), 1), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'bottom', 'horizontalalignment', 'right');
    end
 
    ax7_image = axes('Position', [ax7.Position(1)+0.01 ax7.Position(2)+ax7.Position(4)*0.75 image_width image_width]);
    imshow(bw);
    colormap(ax7_image, flipud(gray))
end

%%% titles
annotation(fig, 'textbox', [0.1, 0.985, 0.8, 0.01], 'string', 'neural timescales', 'edgecolor', 'none', ...
        'fontsize', fontsize_label, 'fontweight', 'b', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
annotation(fig, 'textbox', [0.1, 0.62, 0.8, 0.01], 'string', 'computational capacity', 'edgecolor', 'none', ...
        'fontsize', fontsize_label, 'fontweight', 'b', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    
%%% panel letters
annotation(fig, 'textbox', [0.03, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.03, 0.59, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.62, 0.59, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.62, 0.46, 0.01, 0.01], 'string', 'D', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.03, 0.31, 0.01, 0.01], 'string', 'E', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.54, 0.31, 0.01, 0.01], 'string', 'F', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

%% FIGURE 5

% load all data relevant to Figure 5
data_Figure5 = load(sprintf('%s/Figure5.mat', data_foldername));

types = {'human', 'chimp'};
titles = {'human', 'chimpanzee'};

cmap = [color_brown; color_green];

parc_lh = dlmread(parc_filename('lh'));
parcels_lh = unique(parc_lh(parc_lh>0));
num_parcels_lh = length(parcels_lh);
parc_rh = dlmread(parc_filename('rh'));
parcels_rh = unique(parc_rh(parc_rh>0));
num_parcels_rh = length(parcels_rh);

pvals_adj_ind = find(strcmpi(multiple_comparisons(:,1), 'Figure5'));

image_width = 0.06;

fig = figure('Position', [400 100 700 500]);
for type_ind = 1:length(types)
    type = types{type_ind};
    
    if type_ind==1
        image_to_plot = human_female;
    elseif type_ind==2
        image_to_plot = chimpanzee;
    end
    
    bw = image_to_plot>0;
    
    % =========================================================================
    % A: myelin map spatial arrangement
    % =========================================================================
    myelin_cmap = cbrewer('seq', 'YlGnBu', 130, 'pchip');
    myelin_cmap = [0.5 0.5 0.5; flipud(myelin_cmap(end-N+1:end,:))];
    
    surface_to_plot_lh = gifti(surface_file{type_ind}.lh);
    surface_to_plot_rh = gifti(surface_file{type_ind}.rh);
    
    boundary_method = 'midpoint';
    BOUNDARY_lh = findROIboundaries(surface_to_plot_lh.vertices,surface_to_plot_lh.faces,parc_lh,boundary_method);
    BOUNDARY_rh = findROIboundaries(surface_to_plot_rh.vertices,surface_to_plot_rh.faces,parc_rh,boundary_method);
    
    data_parcel = data_Figure5.tuningcurve.stats{type_ind}.myelin;
    
    data_to_plot_lh = zeros(size(surface_to_plot_lh.vertices,1),1);
    data_to_plot_rh = zeros(size(surface_to_plot_rh.vertices,1),1);
    for ii=1:num_parcels_lh
        parcel = parcels_lh(ii);
        
        data_to_plot_lh(find(parc_lh==parcel)) = data_parcel(ii);
    end
    for ii=1:num_parcels_rh
        parcel = parcels_rh(ii);
        
        data_to_plot_rh(find(parc_rh==parcel)) = data_parcel(num_parcels_lh+ii);
    end
    
%     clims = [min(data_parcel(data_parcel>0)), max(data_parcel(data_parcel>0))];
    clims = prctile(data_parcel(data_parcel>0), [4 96]);
    
    ind_zeros = find(parc_lh==0);
    data_to_plot_lh(data_to_plot_lh<clims(1)) = clims(1);
    data_to_plot_lh(data_to_plot_lh>clims(2)) = clims(2);
    data_to_plot_lh(ind_zeros) = clims(1)-0.01;
    ind_zeros = find(parc_rh==0);
    data_to_plot_rh(data_to_plot_rh<clims(1)) = clims(1);
    data_to_plot_rh(data_to_plot_rh>clims(2)) = clims(2);
    data_to_plot_rh(ind_zeros) = clims(1)-0.01;
    clims(1) = clims(1)-0.01;
    
    ax1 = axes('Position', [0.15+0.43*(type_ind-1) 0.7 0.35 0.3]);
    ax1_positions = ax1.Position;
    set(ax1, 'visible', 'off')
    
    ax1_1 = axes('Position', [ax1_positions(1)+0.01 ax1_positions(2)-0.03 0.16 0.26]);
    patch(ax1_1, 'Vertices', surface_to_plot_lh.vertices, 'Faces', surface_to_plot_lh.faces, 'FaceVertexCData', data_to_plot_lh, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for ii = 1:length(BOUNDARY_lh)
            plot3(BOUNDARY_lh{ii}(:,1), BOUNDARY_lh{ii}(:,2), BOUNDARY_lh{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end
    hold off;
    view([-90 0]);
    caxis(clims)
    material dull
    camlight('headlight');
    colormap(ax1_1, myelin_cmap)
    axis off
    axis image
    
    ax1_1_image = axes('Position', [ax1_positions(1)-0.03 ax1_1.Position(2)+ax1_1.Position(4)*0.66 image_width image_width]);
    imshow(bw);
    colormap(ax1_1_image, flipud(gray))
    
    ax1_2 = axes('Position', [ax1_1.Position(1) ax1_1.Position(2)-ax1_1.Position(4)/1.55 ax1_1.Position(3) ax1_1.Position(4)]);
    patch(ax1_2, 'Vertices', surface_to_plot_lh.vertices, 'Faces', surface_to_plot_lh.faces, 'FaceVertexCData', data_to_plot_lh, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for ii = 1:length(BOUNDARY_lh)
            plot3(BOUNDARY_lh{ii}(:,1), BOUNDARY_lh{ii}(:,2), BOUNDARY_lh{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end
    hold off;
    view([90 0]);
    caxis(clims)
    material dull
    camlight('headlight');
    colormap(ax1_2, myelin_cmap)
    axis off
    axis image
    
    ax1_3 = axes('Position', [ax1_1.Position(1)+ax1_1.Position(3)*1.1 ax1_1.Position(2) ax1_1.Position(3) ax1_1.Position(4)]);
    patch(ax1_3, 'Vertices', surface_to_plot_rh.vertices, 'Faces', surface_to_plot_rh.faces, 'FaceVertexCData', data_to_plot_rh, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for ii = 1:length(BOUNDARY_rh)
            plot3(BOUNDARY_rh{ii}(:,1), BOUNDARY_rh{ii}(:,2), BOUNDARY_rh{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end
    hold off;
    view([90 0]);
    caxis(clims)
    material dull
    camlight('headlight');
    colormap(ax1_3, myelin_cmap)
    axis off
    axis image
    
    ax1_4 = axes('Position', [ax1_3.Position(1) ax1_2.Position(2) ax1_1.Position(3) ax1_1.Position(4)]);
    patch(ax1_4, 'Vertices', surface_to_plot_rh.vertices, 'Faces', surface_to_plot_rh.faces, 'FaceVertexCData', data_to_plot_rh, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
    hold on;
    for ii = 1:length(BOUNDARY_rh)
            plot3(BOUNDARY_rh{ii}(:,1), BOUNDARY_rh{ii}(:,2), BOUNDARY_rh{ii}(:,3), 'Color', 'k', 'LineWidth',1, 'Clipping','off');
    end
    hold off;
    view([-90 0]);
    caxis(clims)
    material dull
    camlight('headlight');
    colormap(ax1_4, myelin_cmap)
    axis off
    axis image
    
    cbar = colorbar(ax1_4,'southoutside');
    ylabel(cbar, 'myelination', 'fontsize', fontsize_label, 'interpreter', 'latex')
    set(cbar, 'fontsize', fontsize_axis, 'ticklength', 0.01, ...
        'position', [ax1_positions(1)+ax1_positions(3)*0.32, ax1_3.Position(2)-0.15, ax1_3.Position(3)*0.8, 0.02], ...
        'ytick', [])
    annotation(fig, 'textbox', [cbar.Position(1)-0.039, cbar.Position(2)*1, 0.04, 0.01], 'string', 'low', 'edgecolor', 'none', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    annotation(fig, 'textbox', [cbar.Position(1)+cbar.Position(3)+0.005, cbar.Position(2)*1, 0.04, 0.01], 'string', 'high', 'edgecolor', 'none', ...
               'fontsize', fontsize_axis, 'horizontalalignment', 'center', 'verticalalignment', 'middle')    
   
    annotation('textbox', [ax1_1.Position(1)-0.04, ax1_1.Position(2), 0.1, 0.1], 'string', 'LH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annotation('textbox', [ax1_1.Position(1)+ax1_1.Position(3)*1.73, ax1_1.Position(2), 0.1, 0.1], 'string', 'RH', ...
       'fontsize', fontsize_axis, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    annotation('textbox', [ax1_1.Position(1)+ax1_1.Position(3)*0.7, ax1_1.Position(2)+ax1_1.Position(4)*0.7, 0.1, 0.1], 'string', titles{type_ind}, ...
       'fontsize', fontsize_label, 'LineStyle', 'none', 'horizontalalignment', 'center', ...
       'verticalalignment', 'middle', 'fontweight', 'b')
    
    % =========================================================================
    % B: correlation of myelin and dynamic range
    % =========================================================================
    ax2 = axes('Position', [ax1.Position(1) ax1.Position(2)-ax1.Position(4)*1.95 ax1.Position(3) ax1.Position(4)]);
    data_to_plot_x = zscore(data_Figure5.tuningcurve.stats{type_ind}.dynamic_range);
    data_to_plot_y = data_Figure5.tuningcurve.stats{type_ind}.myelin;
    hold on;
    plot(data_to_plot_x, data_to_plot_y, '.', 'color', cmap(type_ind,:), 'markersize', 20)
    plot(data_to_plot_x, polyval(polyfit(data_to_plot_x,data_to_plot_y,1), data_to_plot_x), ...
        'k-', 'linewidth', 2);
    hold off;
    set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], 'box', 'off')
    xlabel('dynamic range (z score)', 'fontsize', fontsize_label, 'interpreter', 'latex')
    if type_ind==1
        ylabel('myelination', 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    [rho,pval] = corr(data_to_plot_x, data_to_plot_y, 'type', 'pearson');
    text(max(get(ax2,'xlim')), max(get(ax2,'ylim'))*1, sprintf('r = %.2f', rho), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
    if ~show_padj
        text(max(get(ax2,'xlim')), max(get(ax2,'ylim'))*0.965, utils.extract_pvalue_text(pval), ...
            'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
    else
        text(max(get(ax2,'xlim')), max(get(ax2,'ylim'))*0.965, utils.extract_pvalue_text(pvals_adj(pvals_adj_ind(type_ind)), 1), ...
        'fontsize', fontsize_label-2, 'fontweight', 'bold', 'verticalalignment', 'top', 'horizontalalignment', 'right');
    end
    
    ax2_image = axes('Position', [ax2.Position(1) ax2.Position(2)+ax2.Position(4)*0.85 image_width image_width]);
    imshow(bw);
    colormap(ax2_image, flipud(gray))
end

%%% panel letters
annotation(fig, 'textbox', [0.09, 0.93, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.09, 0.47, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
    
%% FIGURE 6

% load all data relevant to Figure 6
data_Figure6 = load(sprintf('%s/Figure6.mat', data_foldername));
%%
types = {'human', 'macaque'};
titles = {'human', 'macaque'};

cmap = [color_brown; color_blue];

image_width = 0.06;
 
fig = figure('Position', [200 200 600 500]);
ax0 = axes('Position', [0.12 0.62 0.35 0.3]);
ax0_position = ax0.Position;
set(ax0,'visible','off')
 
% =========================================================================
% A: mean FC of large-scale networks
% =========================================================================
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
        network_names = Yeo_names;
    elseif type_ind==2
        image_to_plot = macaque;
        network_names = macaque_network_names;
    end
    bw = image_to_plot>0;
    
    ax1 = axes('Position', [ax0_position(1)+(ax0_position(3)+0.07)*(type_ind-1)+0.03 ax0_position(2) ax0_position(3) ax0_position(4)]);
    bar([data_Figure6.empirical_stats.meanFC{type_ind}.RSN, 0, data_Figure6.empirical_stats.meanFC{type_ind}.global], 'FaceColor', cmap(type_ind,:))
    hold on;
    errorbar(1:(length(data_Figure6.empirical_stats.meanFC{type_ind}.RSN)+1+1), [data_Figure6.empirical_stats.meanFC{type_ind}.RSN, 0, data_Figure6.empirical_stats.meanFC{type_ind}.global], [data_Figure6.empirical_stats.stdFC{type_ind}.RSN, NaN, data_Figure6.empirical_stats.stdFC{type_ind}.global], 'k', 'linestyle', 'none')
    hold off;
    set(ax1, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
        'xtick', [1:length(data_Figure6.empirical_stats.meanFC{type_ind}.RSN), length(data_Figure6.empirical_stats.meanFC{type_ind}.RSN)+2], ...
        'xticklabel', {network_names{:}, '{\bfwhole brain}'}, ...
        'xticklabelrotation', 45, 'ylim', [0 1])
    if type_ind==1
        ylabel('$FC$', 'fontsize', fontsize_label, 'interpreter', 'latex')
    end
    title(titles{type_ind}, 'fontsize', fontsize_label)
    box off
    
    ax1_image = axes('Position', [ax1.Position(1)+0.01 ax1.Position(2)+ax1.Position(4)*0.75 image_width image_width]);
    imshow(bw);
    colormap(ax1_image, flipud(gray))
end

% =========================================================================
% B: functional path length
% =========================================================================
ax2 = axes('Position', [ax0_position(1)+(ax0_position(3)+0.1)*(1-1) 0.05 ax0_position(3) ax0_position(4)*1.2]);
data_violin = utils.padconcatenation(data_Figure6.empirical_stats.pl{1}, data_Figure6.empirical_stats.pl{2}, 2);
violins = utils.violinplot(data_violin, {});
for type_ind=1:length(types)
    violins(type_ind).MedianPlot.SizeData = 30;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 10;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.01;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
end
set(ax2, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
        'xticklabel', titles)
ylabel({'regional functional'; 'path length'}, 'fontsize', fontsize_label, 'interpreter', 'latex')
 
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
        yloc = 0.38;
    elseif type_ind==2
        image_to_plot = macaque;
        yloc = 0.81;
    end
    bw = image_to_plot>0;
 
    ax2_image = axes('Position', [ax2.Position(1)+0.01+(type_ind-1)*0.17 ax2.Position(2)+ax2.Position(4)*yloc image_width image_width]);
    imshow(bw);
    colormap(ax2_image, flipud(gray))
end
 
% =========================================================================
% C: fMRI timescale
% =========================================================================
ax3 = axes('Position', [ax0_position(1)+(ax0_position(3)+0.15)*(2-1) 0.05 ax0_position(3) ax0_position(4)*1.2]);
data_violin = utils.padconcatenation(data_Figure6.empirical_stats.tau_BOLD{1}, data_Figure6.empirical_stats.tau_BOLD{2}, 2);
violins = utils.violinplot(data_violin, {});
for type_ind=1:length(types)
    violins(type_ind).MedianPlot.SizeData = 30;
    violins(type_ind).ViolinColor = cmap(type_ind,:);
    violins(type_ind).ViolinAlpha = 0.2;
    violins(type_ind).ScatterPlot.SizeData = 10;
    violins(type_ind).ScatterPlot.MarkerFaceAlpha = 1;
    violins(type_ind).BoxColor = [0 0 0];
    violins(type_ind).BoxWidth = 0.01;
    violins(type_ind).WhiskerPlot.LineStyle = 'none';
end
set(ax3, 'fontsize', fontsize_axis, 'ticklength', [0.02, 0.02], ...
        'xticklabel', titles)
ylabel('fMRI signal timescale (s)', 'fontsize', fontsize_label, 'interpreter', 'latex')
 
for type_ind = 1:length(types)
    if type_ind==1
        image_to_plot = human_female;
        yloc = 0.44;
    elseif type_ind==2
        image_to_plot = macaque;
        yloc = 0.81;
    end
    bw = image_to_plot>0;
 
    ax3_image = axes('Position', [ax3.Position(1)+0.01+(type_ind-1)*0.17 ax3.Position(2)+ax3.Position(4)*yloc image_width image_width]);
    imshow(bw);
    colormap(ax3_image, flipud(gray))
end
 
%%% panel letters
annotation(fig, 'textbox', [0.08, 0.98, 0.01, 0.01], 'string', 'A', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.04, 0.48, 0.01, 0.01], 'string', 'B', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')
annotation(fig, 'textbox', [0.54, 0.48, 0.01, 0.01], 'string', 'C', 'edgecolor', 'none', ...
        'fontsize', 20, 'fontweight', 'b', 'horizontalalignment', 'center')

