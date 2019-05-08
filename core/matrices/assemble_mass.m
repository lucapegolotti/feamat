function [M] = assemble_mass(fespace, varargin)
% Assemble mass matrix
% input=
%           fespace: finite element space
%           varargin: if empty, no elementlist is considered and the 
%                     reaction coefficient is set to 1,
%                     otherwise the second argument represents the list of 
%                     considered elememts and the first one is a function
%                     handle representing the reaction coefficient
% output=
%           M: matrix (sparse)


if nargin==1
    M = assemble_cu(1, fespace);
elseif nargin == 2
    M = assemble_cu(varargin{1}, fespace);
elseif nargin==3
    M = assemble_cu_elementlist(varargin{1},fespace,varargin{2});
end