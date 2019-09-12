function [seg, Basis] = mesh_seg(mesh, Para)
% Input a trianglar mesh:
% mesh.vertexPoss: a 3xnumV matrix that stores vertex positions
% mesh.faceVIds: a 3xnumF matrix that stores 1-based face vertex indices
% Para.kring : specifies the kring-neighbor that is used to fit the local
%              rigid velocity field the default value is 2
% Para.maxNumEigens: the maximum number of deformation field (each of which correspond
%                 to a eigenvector of the model matrix, the default value is 20
%                 Note that we perform eigen-gap analysis to determine the
%                 optimal one

% Para.maxNumSegs: the maximum number of segments that we used to determine
%                 the optimal number of segments (we do this based on the
%                 fitting error), the default value is 20

% The algorithm proceeds in three steps
% The first phase computes the eigen-basis
Basis = modal_analysis_vectors(mesh, Para.maxNumEigens);
% The second phase fits a velocity field at each vertex for each
% deformation field, we use the fitting error to pick seed points for
% patches
numV = size(mesh.vertexPoss, 2);
for eigenId = 1 : length(Basis.eigenVals)
    vectorfield = reshape(Basis.eigenVecs(:, eigenId), [3, numV]);
    Velocities{eigenId} = vector_field_analysis(mesh, vectorfield, Para.kring);
end
% The last phase applies a greedy approach for segmentation
seg = perform_segmentation(mesh, Basis, Velocities, Para.maxNumSegs);
% Visualization
trisurf(mesh.faceVIds', mesh.vertexPoss(1,:), mesh.vertexPoss(2,:),...
    mesh.vertexPoss(3,:), median(seg(mesh.faceVIds)))
daspect([1,1,1])