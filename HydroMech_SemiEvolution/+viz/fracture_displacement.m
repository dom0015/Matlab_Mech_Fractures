function [] = fracture_displacement(elast_problem,x_elast,fid)
%FRACTURE_DISPLACEMENT Summary of this function goes here
%   Detailed explanation goes here
% fid ... indices of matrices to be visualized
if nargin<3
    fid=1:length(elast_problem.fracture_matrice);
end
for i=fid
    x_elast_demap=x_elast(elast_problem.node_map_on);
    elast = [x_elast_demap(1:2:end), x_elast_demap(2:2:end)];
    above = elast_problem.fracture_matrice{i}.above_nodes;
    above = [above(:,1); above(end,2)];
    under = elast_problem.fracture_matrice{i}.under_nodes;
    under = [under(:,1); under(end,2)];
    LEN = length(elast(above,1));
    
    figure; plot(linspace(3,7,LEN),elast(above,1)); hold on; plot(linspace(3,7,LEN),elast(under,1))
    grid on
    %xlim([2.9, 7.1]); ylim([7.2, 8.7]*1e-5);
    legend('up','down')
    title(['Fracture ' num2str(i) ' displacement (1)']);
    
    figure; plot(linspace(3,7,LEN),elast(above,2)); hold on; plot(linspace(3,7,LEN),elast(under,2))
    grid on
    %xlim([2.9, 7.1]); ylim([-18, -3]*1e-6);
    legend('up','down')
    title(['Fracture ' num2str(i) ' displacement (2)']);
    
    NORM = elast_problem.fracture_normal{i};
    figure; plot(linspace(3,7,LEN),elast(above,:)*NORM); hold on;
    plot(linspace(3,7,LEN),elast(under,:)*NORM);
    legend('up','down')
    grid on
    %xlim([2.9, 7.1]); ylim([-7.1, -5.6]*1e-5)
    title(['Fracture ' num2str(i) ' normal displacement']);
end
end

