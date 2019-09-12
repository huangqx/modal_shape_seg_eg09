function [Basis] = modal_analysis_vectors(mesh, maxNumEigens)
% Compute eigenmodes
edges = extract_edges(mesh);
edgeVec = mesh.vertexPoss(:, edges(1,:)) - mesh.vertexPoss(:, edges(2,:));
%
L = modal_matrix(edges, edgeVec);
% The first six eigenvalues are zero
[eigenVecs, eigenVals] = eigs(L, maxNumEigens + 6, -1e-10);
eigenVals = diag(eigenVals)';

% This is usually application specific
[numEigens] = analysis(eigenVals(7:(maxNumEigens+6)));
for id = 1 : numEigens
    Basis.eigenVals(id) = eigenVals(id+6);
    Basis.eigenVecs(:,id) = eigenVecs(:, id+6);
end
% 
function [num] = analysis(vals)
ratios = vals(2:length(vals))./vals(1:(length(vals)-1));
[s,off] = max(ratios);
num = off;
%
function [L] = modal_matrix(edges, edgeVec)
%
numV = max(max(edges'));
numE = size(edges, 2);
%
rows_J1 = ones(2,1)*(1:numE);
cols_J1 = edges;
vals_J1 = [1,-1]'*ones(1,numE);
J1 = sparse(rows_J1, cols_J1, vals_J1, numE, numV);
J1 = kron(J1, sparse(eye(3)));
%
rows_J2 = ones(3,1)*(1:(3*numE));
sIds = edges(1,:);
cols_J2 = kron([3*sIds-2;3*sIds-1;3*sIds], ones(1,3));
vals_J2 = zeros(3, 3*numE);
vals_J2(1, 2:3:(3*numE)) = edgeVec(3,:);
vals_J2(1, 3:3:(3*numE)) = -edgeVec(2,:);
vals_J2(2, 1:3:(3*numE)) = -edgeVec(3,:);
vals_J2(2, 3:3:(3*numE)) = edgeVec(1,:);
vals_J2(3, 1:3:(3*numE)) = edgeVec(2,:);
vals_J2(3, 2:3:(3*numE)) = -edgeVec(1,:);
J2 = sparse(rows_J2, cols_J2, vals_J2, 3*numE, 3*numV);
%
vals_invD = zeros(3, 3*numV);
for eId = 1 : numE
    ids = (3*sIds(eId)-2):(3*sIds(eId));
    e = edgeVec(:, eId);
    vals_invD(:,ids) = vals_invD(:,ids) + (e'*e)*eye(3) - e*e';
end
for vId = 1 : numV
    ids = (3*vId-2):(3*vId);
    vals_invD(:, ids) = inv(vals_invD(:, ids));
end
rows_invD = ones(3,1)*(1:(3*numV));
cols_invD = kron([3*(1:numV)-2;3*(1:numV)-1;3*(1:numV)], ones(1,3));
invD = sparse(rows_invD, cols_invD, vals_invD, 3*numV, 3*numV);
A = J1'*J1;
B = J1'*J2;
L = A - B*invD*B';
L = (L+L')/2;

function [edges] = extract_edges(mesh)
%
v1Ids = mesh.faceVIds(1, :);
v2Ids = mesh.faceVIds(2, :);
v3Ids = mesh.faceVIds(3, :);
rows = [v1Ids, v2Ids, v3Ids];
cols = [v2Ids, v3Ids, v1Ids];
numV = size(mesh.vertexPoss, 2);
A = sparse(rows, cols, ones(1, length(rows)), numV, numV);
A = A + A';
[rows, cols,vals] = find(A);
edges = [rows,cols]';