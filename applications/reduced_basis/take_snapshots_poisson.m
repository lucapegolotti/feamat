function [S] = take_snapshots_poisson(fespace, rangemu, nmus, rangeomega, nomegas)
% Take snapshots for the poisson problem where the parameter is the
% diffusion coefficient (constant over the domain). The additional
% parameter sets the integral on the non constrained boundary
% input= 
%           fespace: finite element space
%           rangemu: range of sampling of the diffusion parameter
%           nmus: number of samples taken from rangemu
%           rangeomega: range of sampling for omega (integral over non
%               constrained boundary)
%           nomegas: number of samples taken from rangeomega
% output=
%           S: matrix of snapshots




end

