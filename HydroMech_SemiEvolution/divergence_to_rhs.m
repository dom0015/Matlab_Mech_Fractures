function [b] = divergence_to_rhs( hydro_problem, divergence_diff )

node = hydro_problem.POINTS;
elem = hydro_problem.ELEMENTS;

N=size(node,1);     % number of nodes
NT=size(elem,1);    % number of elements

%% VECTORIZATION (EDGES + AREA)
ve=zeros(NT,2,3);
ve(:,:,1)=node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2)=node(elem(:,1),:)-node(elem(:,3),:);
ve(:,:,3)=node(elem(:,2),:)-node(elem(:,1),:);

area=0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));

%% RIGHT HAND SIDE
bt1=area.*divergence_diff/3;
bt2=area.*divergence_diff/3;
bt3=area.*divergence_diff/3;

b=accumarray(elem(:),[bt1;bt2;bt3],[N 1]);