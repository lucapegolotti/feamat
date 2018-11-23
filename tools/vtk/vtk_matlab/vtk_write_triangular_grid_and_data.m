function stat = vtk_write_triangular_grid_and_data(filename,vtk_title,grid_X,grid_TRI,data_struct,flipped)
% modified version of vtk_write_tetrahedral_grid and data by 
% Shawn W. Walker (2016) to write 2D data on unstructured grids. Refer
% to the documentation of vtk_write_tetrahedral_grid_and_data for
% additional information

% define write commands
write_ascii = @(fid,str) [fprintf(fid,str); fprintf(fid,'\n');];
write_data  = @(fid,dat,prec) [fwrite(fid,dat,prec); fprintf(fid,'\n');];

% check the filename
[P1, N1, E1] = fileparts(filename);
if isempty(E1)
    E1 = '.vtk'; % add the extension
elseif ~strcmp(E1,'.vtk')
    disp('warning: file extension is *not* ".vtk"!');
end
filename = fullfile(P1, [N1, E1]);

% open the file in binary mode
fopen_opts = {'wb','ieee-be'};
fid = fopen(filename,fopen_opts{:});
if fid == -1
    error('Unable to write file %s: permission denied.',filename);
end

% write the initial header
write_ascii(fid,'# vtk DataFile Version 2.0');
vtk_title_reduce = vtk_title(1:min(length(vtk_title),256));
write_ascii(fid,vtk_title_reduce);
write_ascii(fid,'BINARY\n'); % give extra line-feed
% write the vertex coordinates of the mesh
write_ascii(fid,['DATASET ', 'UNSTRUCTURED_GRID']);
if ~flipped
    grid_X = grid_X'; % need to transpose for fwrite
end
GD = size(grid_X,1);
if (GD~=3)
    stat = fclose(fid);
    error('Grid vertex coordinate data is not 3-D!');
end
Num_Vtx = size(grid_X,2);
write_ascii(fid, ['POINTS ', num2str(Num_Vtx), ' float']);
write_data(fid, grid_X, 'float32');

% write the triangular grid connectivity
if ~flipped
    Tri_Order = size(grid_TRI,2);
    Num_Tri   = size(grid_TRI,1);
else
    Tri_Order = size(grid_TRI,1);
    Num_Tri   = size(grid_TRI,2);
end
if (Tri_Order~=3)
    stat = fclose(fid);
    error('Grid triangles connectivity does not have 3 nodes per tet!');
end
Cell_Size = Num_Tri * ( Tri_Order + 1 ); % vtk needs this
write_ascii(fid,['CELLS  ', num2str(Num_Tri), '  ', num2str(Cell_Size)]);

% 1-based to 0-based indexing
min_index = min(grid_TRI(:));
if min_index==0
    disp('Triangle connectivity is already 0-based!');
elseif min_index > 0
    grid_TRI = grid_TRI - 1;
else
    fclose(fid);
    error('Triangle connectivity has negative indices!');
end

% modify and append extra data
if ~flipped
    DATA = uint32([(Tri_Order + 0*grid_TRI(:,1))';
                    grid_TRI']);
else
    DATA = uint32([(Tri_Order + 0*grid_TRI(1,:));
                    grid_TRI]);
end
write_data(fid,DATA, 'uint32');

% must write the cell types
% VTK has a cell type 5 for linear triangles.  
write_ascii(fid,['CELL_TYPES ', num2str(Num_Tri)]);
if ( Tri_Order == 3 )
    TRI_LABEL = 5; 
else
    error('Invalid triangle order!');
end
DATA = uint32(TRI_LABEL*ones(1,Num_Tri));
write_data(fid,DATA, 'uint32');

% write the POINT_DATA
if ~isempty(data_struct)
    write_ascii(fid,['POINT_DATA ', num2str(Num_Vtx)]);
    
    Num_DS = length(data_struct);
    for ii = 1:Num_DS
        if strcmpi(data_struct(ii).type,'scalar')
            write_ascii(fid,['SCALARS ', data_struct(ii).name, ' float']);
            write_ascii(fid,'LOOKUP_TABLE default');
            write_data(fid,data_struct(ii).data(:)','float32'); % scalar data easy to write
        elseif strcmpi(data_struct(ii).type,'vector')
            write_ascii(fid,['VECTORS ', data_struct(ii).name, ' float']);
            if ~flipped
                D1 = data_struct(ii).data';
            else
                D1 = data_struct(ii).data;
            end
            write_data(fid,D1,'float32');
        else
            stat = fclose(fid);
            error('Invalid data type!');
        end
    end
end

% Note: CELL_DATA is not implemented!

% end of file!
stat = fclose(fid);

end