function [edgeStruct]=G7_CreateGridUGMModel(NumFils, NumCols, K )
%
% NumFils, NumCols: image dimension
% K: number of states
nNodes = NumFils * NumCols;


adj = sparse(nNodes,nNodes);
 
% %BEGIN TODO
% Add Down Edges
% TODO
ind = 1:nNodes;
%We dont have down edge at the last row
exclude = sub2ind([NumFils NumCols],repmat(NumFils,[1 NumCols]),1:NumCols); 
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;
 
% Add Right Edges
% TODO
ind = 1:nNodes;
exclude = sub2ind([NumFils NumCols],1:NumFils,repmat(NumCols,[1 NumFils])); 
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+NumFils)) = 1;
 
% Add Up/Left Edges
% TODO
adj = adj + adj';

%End of TODO


edgeStruct = UGM_makeEdgeStruct(adj, K);
