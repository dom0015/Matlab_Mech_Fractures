close all
clear all
addpath(genpath('files'))
%% %% independent of fracture aperture
geometry_assembly
tocouple_aperture_independent

%% %% depends on fracture aperture
d = exp([-6 -8.5 -7]);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end

[PRESSURE,u0_,GRAD]=tocouple_handle(D,no_fractures,mat_frac,fracture_matrice,...
    POINTS,intersections,alfa_inter,lengths,A,freeNode,b,u0,ELEMENTS);
PRESSURE{1}'
[PRESSURE,u0_,GRAD]=tocouple_handle(D,no_fractures,mat_frac,fracture_matrice,...
    POINTS,intersections,alfa_inter,lengths,A,freeNode,b,u0,ELEMENTS);
PRESSURE{1}'
%% VISUALIZATION ETC.
Q = extract_flow( A, u0_(1:length(A)), POINTS, ELEMENTS, Nxy-1, Dirichlet_windows, downEdge, rightEdge, upEdge, leftEdge );
startx=ones(10,1);
starty=linspace(5/10,10/10,10);
[ tmp1,tmp2,xx,yy,ftmp1,ftmp2,fxx,fyy,node_fluxx,node_fluxy] = streamlines_calc( u0_(1:length(A)),POINTS,double(ELEMENTS),linspace(0.1/10,9.9/10,23),linspace(0.1/10,9.9/10,23),...
    linspace(0,1,1000),linspace(0,1,1000),fractures_positions);

% Q = Q*mat_omega/h_elem;
disp(Q); 
%disp(sum(Q))

%% FIGURE
figure;
N=length(A);
h = trisurf(ELEMENTS,POINTS(:,1),POINTS(:,2),u0_(1:N));
alpha 0.8
h.EdgeColor = 'none';
colormap jet(1000)
colorbar
axis equal
hold on
view(0,90)
%set(gca,'ZDir','reverse')
quiver3(xx,yy,xx*0+p,-tmp1,-tmp2,tmp2*0,1.1,'Color','black','LineWidth',1.4)
%[x_g,y_g]=meshgrid(linspace(1,9,5),linspace(1,9,5));
n_sl=30;
x_g=[ones(n_sl,1)/10;0.9*ones(n_sl,1)];
y_g=[linspace(0.1,9.9,n_sl)';linspace(0.1,9.9,n_sl)']/10;
frac_points=[x_g(:),y_g(:)];

hlines=streamline(fxx,fyy,ftmp1,ftmp2,frac_points(:,1),frac_points(:,2),[1e-2 1e7]);
set(hlines,'LineWidth',1.5,'Color','k');%[183 0 255]/256)

hlines=streamline(fxx,fyy,-ftmp1,-ftmp2,frac_points(:,1),frac_points(:,2),[1e-2 1e7]);
set(hlines,'LineWidth',1.5,'Color','k');%[183 0 255]/256)

%quiver(xx,yy,tmp1,tmp2,1.5,'Color','black','LineWidth',1.5)
%streamline(centers(:,1),centers(:,2),grad_elem(:,1),grad_elem(:,2),startx,starty)
%contour(node(:,1),node(:,2),u0(1:length(A)),10)
M = N;
for i=1:length(lengths)
    coord_x = fractures_positions{i}(:,1);
    coord_y = fractures_positions{i}(:,2);
    vals = u0_(M+1:M+lengths(i));
    plot3(coord_x,coord_y,vals+1e6,'k','LineWidth',2); hold on
    diff_x = coord_x(2:end) - coord_x(1:end-1);
    diff_y = coord_y(2:end) - coord_y(1:end-1);
    steps_arg = sqrt(diff_x.^2 + diff_y.^2);
    steps_val = vals(2:end) - vals(1:end-1);
    %figure; plot(cumsum(steps_arg),steps_val)
    %figure; plot(cumsum(steps_arg),steps_arg)
    %figure; plot(cumsum(steps_arg),mat_fract(i)*steps_val./steps_arg)
    %figure; plot([0; cumsum(steps_arg)],vals)
    M=M+lengths(i);
end