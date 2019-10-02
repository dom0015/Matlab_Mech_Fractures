function [divergence] = elasticity_divergence(problem_setting,displacement,print)

NODE = cell2mat(problem_setting.sub_nodes);
coord_x = NODE(:,1);
coord_y = NODE(:,2);
ELEM = problem_setting.sub_elem;
n=0;
for i=1:length(ELEM)
    ELEM{i} = ELEM{i} + n;
    n = n + size(problem_setting.sub_nodes{i},1);
end
ELEM=cell2mat(ELEM);
disp_x = displacement(1:2:end);
disp_y = displacement(2:2:end);

a_x=coord_x(ELEM); a_y=coord_y(ELEM);
x1=a_x(:,1); x2=a_x(:,2); x3=a_x(:,3); 
y1=a_y(:,1); y2=a_y(:,2); y3=a_y(:,3); 
u_x=disp_x(ELEM);
u_y=disp_y(ELEM);
ux1=u_x(:,1); ux2=u_x(:,2); ux3=u_x(:,3);
uy1=u_y(:,1); uy2=u_y(:,2); uy3=u_y(:,3);
divergence_x=-(y1.*ux2 - y2.*ux1 - y1.*ux3 + y3.*ux1 + y2.*ux3 - y3.*ux2)./(x1.*y2 - x2.*y1 - x1.*y3 + x3.*y1 + x2.*y3 - x3.*y2);
divergence_y= (x1.*uy2 - x2.*uy1 - x1.*uy3 + x3.*uy1 + x2.*uy3 - x3.*uy2)./(x1.*y2 - x2.*y1 - x1.*y3 + x3.*y1 + x2.*y3 - x3.*y2);
divergence = divergence_x + divergence_y;

if print
    figure(104); plot_divergence(divergence,NODE,ELEM)
    title("Divergence")
end

% remap to hydro elem
map = problem_setting.map;
[~,ord] = sort(map);
divergence(ord) = divergence;
