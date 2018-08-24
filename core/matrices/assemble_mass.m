function [M] = assemble_mass(fespace)
% Assemble mass matrix
% input=
%           fespace: finite element space
% output=
%           M: matrix (sparse)

M = assemble_cu(1,fespace);