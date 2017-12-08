function [M] = assemble_mass(fespace)

M = assemble_cu(1,fespace);