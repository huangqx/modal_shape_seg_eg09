function [Result] = vector_field_analysis(mesh, vectorfield, kring)
% Perform per-vertex analysis
edges = extract_edges(mesh);
numV = max(max(edges));
vectorfield = reshape(vectorfield, [3, numV]);
numE = size(edges, 2);
A = sparse(edges(1, :), edges(2,:), ones(1,numE), numV, numV);
A = max(A, A');
B = A;
for i = 1 : (kring-1)
    B = A + A*B;
end
%
Result = zeros(7, numV);
for vId = 1 : numV
    nIds = find(B(vId,:));
    [barc,c, e] = vf_fitting(vectorfield(:, nIds), mesh.vertexPoss(:, nIds));
    Result(1, vId) = e;
    Result(2:4, vId) = barc;
    Result(5:7, vId) = c;
end

%
function [barc, c, e] = vf_fitting(vfs, vertexPoss)
%
numP = size(vfs, 2);
J = zeros(3*numP, 6);
g = reshape(vfs, [3*numP,1]);
J(:,1:3) = kron(ones(numP,1), eye(3));
J(1:3:(3*numP), 5) = vertexPoss(3,:)';
J(1:3:(3*numP), 6) = -vertexPoss(2,:)';
J(2:3:(3*numP), 4) = -vertexPoss(3,:)';
J(2:3:(3*numP), 6) = vertexPoss(1,:)';
J(3:3:(3*numP), 4) = vertexPoss(2,:)';
J(3:3:(3*numP), 5) = -vertexPoss(1,:)';
A = J'*J;
b = J'*g;
x = A\b;
e = J*x - g;
barc = x(1:3);
c = x(4:6);
e = sum(e'*e)/sum(sum(vfs.*vfs));

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