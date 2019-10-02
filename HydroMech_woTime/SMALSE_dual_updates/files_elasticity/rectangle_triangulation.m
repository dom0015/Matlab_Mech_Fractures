function [ ELEMENTS ] = rectangle_triangulation( nx,ny )
%TRIANGULACE Triangulation of a rectangular domain.
%   nx,ny   numbers of nodes on both sides
idx=reshape(1:nx*ny,ny,nx);
idx1=idx(1:end-1,1:end-1); idx1=idx1(:)';
idx2=idx(1:end-1,2:end); idx2=idx2(:)';
idx3=idx(2:end,2:end); idx3=idx3(:)';
idx4=idx(2:end,1:end-1); idx4=idx4(:)';
EL1=[idx1; idx1];
EL2=[idx2; idx3];
EL3=[idx3; idx4];
ELEMENTS=[EL1(:) EL2(:) EL3(:)];
end

