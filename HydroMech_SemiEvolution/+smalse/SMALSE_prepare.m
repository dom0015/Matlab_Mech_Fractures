function [elast_problem] = SMALSE_prepare(elast_problem)
%SMALSE_SOLVER Summary of this function goes here
%   Detailed explanation goes here
K_plus=elast_problem.A_plus;%*problem_setting.mat_scale;
B_e=elast_problem.B_e;
b=(elast_problem.b);%/problem_setting.mat_scale;
R=elast_problem.R;
c_e=elast_problem.c_e;
B_iupdate=elast_problem.B_iupdate;

if isempty(elast_problem.B_i)
    B_i=B_iupdate(0*b);
else
    B_i=elast_problem.B_i;
end
NODES = cell2mat(elast_problem.sub_nodes)';
NODES = NODES(:);
c_i=B_i*NODES;
if ~isempty(elast_problem.c_i)
    c_i=c_i-elast_problem.c_i;
end
elast_problem.c=[c_e;c_i];
elast_problem.B=[B_e;-B_i];

elast_problem.F=elast_problem.B*K_plus*elast_problem.B';
elast_problem.G=R'*elast_problem.B';
elast_problem.GGt=elast_problem.G*elast_problem.G';
elast_problem.Q=elast_problem.G'*elast_problem.G;
elast_problem.B_K_plus=elast_problem.B*K_plus;

elast_problem.idx_no_bounds=1:length(c_e);
elast_problem.idx_bounds=length(c_e)+1:length(elast_problem.c);

elast_problem.update_G=@(x,y,data)smalse.Update_all_dual_new(x,elast_problem.idx_no_bounds,elast_problem.idx_bounds,elast_problem.c,B_e,B_i,K_plus,y,R,data);
elast_problem.NODES=NODES;
elast_problem.B_i=B_i;

elast_problem.norms_mat=eigs(elast_problem.F,1);
end

