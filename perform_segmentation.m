function [seg] = perform_segmentation(mesh, Basis, Velocities, maxNumSegs)
% The weights used to combine different deformation fields
weights = Basis.eigenVals/Basis.eigenVals(1);

% Compute the fitting errors as the saliency score at each vertex
% Low fitting error means good candidates for the seed vertices
numV = size(mesh.vertexPoss, 2);
saliencies = zeros(1, numV);
for id = 1 : length(weights)
    saliencies = saliencies + weights(id)*Velocities{id}(1,:);
end
[s, next_seed_id] = min(saliencies);
seedIds = [];
closestSeedIds = -ones(1, numV);
seedDistances = 1e5*ones(1, numV);
curves = zeros(1, maxNumSegs);
for id = 1 : maxNumSegs
    % Use the velocity at each vertex to express vectors at other vertices,
    % small expression error means high likely-hood of being in the same
    % rigid clsuter
    errorMat = zeros(length(weights), numV);
    for eigenId = 1 : length(weights)
        errorMat(eigenId, :) = description_error(...
            mesh.vertexPoss,...
            Velocities{eigenId}(2:7, next_seed_id),...
            Basis.eigenVecs(:, eigenId));
    end
    % Add the results across all eigenvectors
    culErrorWeight = weights*errorMat;
    seedIds = [seedIds, next_seed_id];
    ids = find(culErrorWeight < seedDistances);
    seedDistances(ids) = culErrorWeight(ids);
    closestSeedIds(ids) = id;
    % Compute the next one seed vertex, which has low saliency and large
    % expression error with respect to existing seed vertices
    curves(id) = max(seedDistances);
    activeSeedIds = find(seedDistances > curves(id)/2);
    [s,offset] = min(saliencies(activeSeedIds));
    next_seed_id = activeSeedIds(offset);
    history{id} = closestSeedIds;
end
numSegs = curve_analysis(curves);
seg = history{numSegs};
%

function [optimalNumSegs] = curve_analysis(curves)
% Extract the 'turning' points
numS = length(curves);
ratios = zeros(1, numS);
ratios(2:numS) = curves(1:(numS-1))./curves(2:numS);
validIds = find(curves<min(curves)*8);
[s,off] = max(ratios(validIds));
optimalNumSegs = validIds(off);


function [errorVec] = description_error(vertexPoss, vel, vectorfield)
numV = size(vertexPoss, 2);
vectorfield = reshape(vectorfield, [3, numV]);
barc = vel(1:3);
c = vel(4:6);
dif = vectorfield - kron(barc, ones(1,numV))...
    - cross(kron(c, ones(1,numV)), vertexPoss);
errorVec = sum(dif.*dif);