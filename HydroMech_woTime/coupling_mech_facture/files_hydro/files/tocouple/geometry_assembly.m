%% GEOMETRY PARAMETERS - defined in elasticity
L1=10; L2=10;
Nxy=101;
frac_start_end={[0.2 0.4], [0.8 0.4]
                [0.2 0.2], [0.8 0.8]
                [0.2 0.3], [0.8 0.3]};
            
%% GEOMETRY ASSEMBLING - defined in elasticity
[POINTS,ELEMENTS,Dirichlet_boundaries,Neumann_boundaries,...
    Neumann_normalx,Neumann_normaly,fractures,fractures_positions,...
    no_intersections,fractures_cell,fracture_matrice,fracture_elem_map]...
    = create_geometry(Nxy,L1,L2,frac_start_end);