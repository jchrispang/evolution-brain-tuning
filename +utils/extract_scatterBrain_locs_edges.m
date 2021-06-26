function [locs, edges] = extract_scatterBrain_locs_edges(nodeLocations, edge_X, edge_Y, edge_Z, slice)
% extract_scatterBrain_locs_edges.m
%
% Extract locations and edges specific to brain view
%
% Inputs: nodeLocations : 3D locations of the nodes [Nx3]
%         edge_X        : edges along the x-axis [Mx1]
%         edge_Y        : edges along the y-axis [Mx1]
%         edge_Z        : edges along the z-axis [Mx1]
%         slice         : view slice (string)
%                         'axial', 'sagittal_left', 'sagittal_right',
%                         'coronal'
%
% Outputs: locs         : locations of the nodes according to chosen brain view [nx3]
%          edges        : edges according to chosen brain view [Mx3]
%
% Original: James Pang, QIMR Berghofer, 2021
% Revised:  James Pang, Monash University, 2021

%%
if nargin<5
    slice = 'axial';
end

if strcmpi(slice, 'axial')
    locs = cat(2, nodeLocations(:,1), nodeLocations(:,2), nodeLocations(:,3));
    edges = cat(2, edge_X, edge_Y, edge_Z);
elseif strcmpi(slice, 'sagittal_left')
    locs = cat(2, -nodeLocations(:,2), nodeLocations(:,3), -nodeLocations(:,1));
    edges = cat(2, -edge_Y, edge_Z, -edge_X);
elseif strcmpi(slice, 'sagittal_right')
    locs = cat(2, nodeLocations(:,2), nodeLocations(:,3), nodeLocations(:,1));
    edges = cat(2, edge_Y, edge_Z, edge_X);
elseif strcmpi(slice, 'coronal')
    locs = cat(2, nodeLocations(:,1), nodeLocations(:,3), -nodeLocations(:,2));
    edges = cat(2, edge_X, edge_Z, -edge_Y);
end
