function [POINTS,ELEMENTS,Dirichlet_boundaries,Neumann_boundaries,Neumann_normalx,Neumann_normaly,...
    fractures,fractures_positions,no_intersections,fractures_cell,fracture_matrice,fracture_elem_map,intersections] = create_geometry(Nxy,L1,L2,frac_start_end)
%CREATE_GEOMETRY Summary of this function goes here
%   Detailed explanation goes here
nx=Nxy; ny=Nxy;
[coords1,coords2]=meshgrid(linspace(0,L1,nx),linspace(0,L2,ny));
coords1=coords1(:); coords2=coords2(:);
POINTS=[coords1,coords2];
ELEMENTS=mesh.rectangle_triangulation(nx,ny);
% NEUMANN - boundary edges
Neumann_boundaries=ELEMENTS*0;
Neumann_normalx=ELEMENTS*0;
Neumann_normaly=ELEMENTS*0;
Neumann_boundaries((1:((Nxy-1)*2):end),1)=1;
Neumann_normaly((1:((Nxy-1)*2):end),1)=1;
Neumann_boundaries((end-(Nxy-1)*2+1):2:end,2)=2;
Neumann_normalx((end-(Nxy-1)*2+1):2:end,2)=-1;
Neumann_boundaries(((Nxy-1)*2):((Nxy-1)*2):end,2)=3;
Neumann_normaly(((Nxy-1)*2):((Nxy-1)*2):end,2)=-1;
Neumann_boundaries((2:2:(Nxy-1)*2),3)=4;
Neumann_normalx((2:2:(Nxy-1)*2),3)=1;
% DIRICHLET - boundary nodes
Dirichlet_boundaries=false(4,length(POINTS));
Dirichlet_boundaries(1,1:Nxy:end)=true;
Dirichlet_boundaries(2,(end-Nxy+1):end)=true;
Dirichlet_boundaries(3,Nxy:Nxy:end)=true;
Dirichlet_boundaries(4,1:Nxy)=true;
%% Add fractures geometry -------------------------------------------------
[fractures, fractures_positions, no_fractures] = mesh.create_fractures( frac_start_end, POINTS, nx-1 );
[fractures_cell,fracture_matrice,intersections,lengths] = mesh.fracture2cells_geometry( fractures );
no_intersections = size(intersections,1);
[ POINTS,ELEMENTS,Neumann_boundaries,fractures_cell,fracture_matrice,Dirichlet_boundaries]...
    = mesh.multi_fracture_tear_FETI( POINTS,ELEMENTS,fractures_cell ,Neumann_boundaries,fracture_matrice,Dirichlet_boundaries);
% for i=1:1
%     [ POINTS,ELEMENTS,Neumann_boundaries,fractures_cell,fracture_matrice,...
%         Dirichlet_boundaries,Neumann_normalx,Neumann_normaly,fractures]...
%         = smooth_fracture( POINTS,ELEMENTS,fractures_cell ,Neumann_boundaries,fracture_matrice,...
%         Dirichlet_boundaries,Neumann_normalx,Neumann_normaly,fractures);
% end
[frac_bdflag,frac_normx,frac_normy,fracture_elem_map] = mesh.map_fract2elem(POINTS,ELEMENTS,fracture_matrice);
Neumann_boundaries=Neumann_boundaries+frac_bdflag;
Neumann_normalx=Neumann_normalx+frac_normx;
Neumann_normaly=Neumann_normaly+frac_normy;

end

