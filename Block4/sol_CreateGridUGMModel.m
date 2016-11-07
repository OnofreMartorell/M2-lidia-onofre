function [edgeStruct]=CreateGridUGMModel(NumFils, NumCols, K )
%
% NumFils, NumCols: image dimension
% K: number of states




nNodes = ???;


adj = sparse(??);

%BEGIN TODO
% Add Down Edges
TODO
.
.
.

% Add Right Edges
TODO
.
.
.

% Add Up/Left Edges
TODO
.
.
.
%End of TODO

edgeStruct = UGM_makeEdgeStruct(adj,K);
