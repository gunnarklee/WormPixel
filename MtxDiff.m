function [DiffPts]= MtxDiff (mtx1, mtx2)
%find the points that differ between two BINARY matricies.
compMtx=mtx1-mtx2
[DiffPtsR,DiffPtsC]=find(compMtx)
DiffPts=[DiffPtsR,DiffPtsC]