function [M] = assemble_mass(fespace, varargin)
% Assemble mass matrix
% input=
%           fespace: finite element space
%           varargin: if empty, the reation oeffiient is set to 1,
%           otherwise a funtion handle expressing the reaction coefficient
%           is present
% output=
%           M: matrix (sparse)


if nargin==1
    M = assemble_cu(1, fespace);
else
    M = assemble_cu(varargin{1}, fespace);
end