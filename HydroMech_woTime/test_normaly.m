x=[0; 0];
y=rand(2,1);
z=[1; 1];

M=[0 1; -1 0];

n_xy=M*(y-x);
n_yz=M*(z-y);

n_average=(n_xy+n_yz)/2;
n_xz=M*(z-x)/2;

figure
plot([x(1) y(1) z(1)], [x(2) y(2) z(2)], '*'); hold on
plot([y(1) y(1)+n_xy(1)],[y(2) y(2)+n_xy(2)], 'r'); 
plot([y(1) y(1)+n_yz(1)],[y(2) y(2)+n_yz(2)], 'm'); 
plot([y(1) y(1)+n_xz(1)],[y(2) y(2)+n_xz(2)], 'k'); 
axis equal

figure
plot([x(1) y(1) z(1)], [x(2) y(2) z(2)], '*'); hold on
plot([y(1) y(1)+n_xy(1)],[y(2) y(2)+n_xy(2)], 'r'); 
plot([y(1) y(1)+n_yz(1)],[y(2) y(2)+n_yz(2)], 'm'); 
plot([y(1) y(1)+n_average(1)],[y(2) y(2)+n_average(2)], 'k'); 
axis equal
