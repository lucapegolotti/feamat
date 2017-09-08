function [M] = assemble_mass(fespace)

M = assemble_advection(@(x) 1,fespace);