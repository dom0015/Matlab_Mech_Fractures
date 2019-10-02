function [hydro_problem] = hydro_preparation( hydro_problem,shared_data )
frac_start_end = shared_data.frac_start_end;
Nxy = shared_data.Nxy;
L1 = shared_data.L1;
Neumann_boundaries = shared_data.Neumann_boundaries;
fractures = shared_data.fractures;
POINTS = shared_data.POINTS;
ELEMENTS = shared_data.ELEMENTS;
no_intersections = shared_data.no_intersections;
fracture_matrice = shared_data.fracture_matrice;
mat_omega_const = hydro_problem.mat_omega_const;
mat_frac_const = hydro_problem.mat_frac_const;
alfa_inter_const = hydro_problem.alfa_inter_const;
cislo_ulohy = hydro_problem.cislo_ulohy;
DIRICHLET_PRESSURE = hydro_problem.DIRICHLET_PRESSURE;
n_windows = hydro_problem.n_windows;
hydro_model = hydro_problem.hydro_model;

%% FORCE AND BOUNDARY
f = @(x)(0+0*x(:,1)+0*x(:,2)); % zatizeni
g_N=@(x)(0+0*x(:,1)+0*x(:,2)); % Neumann

%% ONLY FOR HYDRO
no_fractures=size(frac_start_end,1);
bdFlag=Neumann_boundaries;
bdFlag(bdFlag<0)=0;
bdFlag=bdFlag(:,[2,3,1]);
intersections=mesh.find_intersections(fractures);
lengths=zeros(1,no_fractures);
for i=1:no_fractures
    lengths(i)=length(fractures{i});
end
k = mat_omega_const*ones(length(ELEMENTS),1);
mat_frac = mat_frac_const*ones(no_fractures,1); % material - fractures
alfa_inter = alfa_inter_const*ones(no_intersections,1);

%% Dirichletova okna ------------------------------------------------
% V kazdem radku:
% prvni hodnota - 1 dole; 2 vpravo; 3 nahore; 4 vlevo
% druha hodnota = a ... zacatek okna 
% treti hodnota = b ... konec okna ... 0 < a < b <= 1
% ctvrta hodnota - hodnota Dirichletovy podminky

Dirichlet_windows = a_hyd.D_windows( cislo_ulohy,DIRICHLET_PRESSURE,n_windows,0);
Dirichlet_windows(3:end,2)=Dirichlet_windows(3:end,2)+1/(Nxy+1);
Dirichlet_windows(:,[2,3])=Dirichlet_windows(:,[2,3])*L1;

%% MATRICES ASSEMBLING
[u0, A, b, freeNode, downEdge, rightEdge, upEdge, leftEdge, M ] = a_hyd.FEM_windows( POINTS, ELEMENTS, Nxy-1, bdFlag, k, f, Dirichlet_windows, g_N );

hydro_problem.no_fractures=no_fractures;
hydro_problem.mat_frac=mat_frac;
hydro_problem.POINTS=POINTS;
hydro_problem.intersections=intersections;
hydro_problem.alfa_inter=alfa_inter;
hydro_problem.lengths=lengths;
hydro_problem.A=A;
hydro_problem.M=M;
hydro_problem.freeNode=freeNode;
hydro_problem.b=b;
hydro_problem.u0=u0;
hydro_problem.ELEMENTS=ELEMENTS;

hydro_problem.Dirichlet_windows=Dirichlet_windows;
hydro_problem.downEdge=downEdge;
hydro_problem.rightEdge=rightEdge;
hydro_problem.upEdge=upEdge;
hydro_problem.leftEdge=leftEdge;
hydro_problem.helem=Nxy-1;
hydro_problem.fracture_matrice=fracture_matrice;
hydro_problem.model=hydro_model;