d=load('uloha.mat');
freeNode=d.freeNode_u;
Km=d.A(freeNode,freeNode);
Mm=d.M(freeNode,freeNode);
Kf=d.F_stif;
Mf=d.F_mass;
M=d.Au;
Bm=-d.B(:,freeNode);
Bf=d.G';
bm=d.b(freeNode);
u_old=d.u_old;

cm=1e-10;
T=1000;
cf=1e-10;

Am=Km+cm/T*Mm;
Af=Kf+cf/T*Mf;

n1=length(Am);
n2=length(Af);
n3=length(M);
O=zeros(n1,n2);

MAT=[Am  O   Bm'
    O'  Af  Bf'
    Bm  Bf -M];
MAT_time=[cm/T*Mm   O       0*Bm'
    O'        cf/T*Mf 0*Bf'
    0*Bm      0*Bf    0*M];

b_orig=[bm;zeros(n2+n3,1)];
b=b_orig+MAT_time*u_old;


[ y ] = my_shur_solver1( b,Am,Af,-M,Bm,Bf );
norm(MAT*y-b)

