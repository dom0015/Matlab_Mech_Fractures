addpath("/opt/petsc/linux-c-opt/share/petsc/matlab");
A = PetscBinaryRead('~/devel/permon/examples/Membranes_F');
B = PetscBinaryRead("~/devel/permon/examples/Membranes_G");
b = PetscBinaryRead("~/devel/permon/examples/Membranes_b0");

%[x,flg,k] = smalxe(A,b,B,zeros(size(B,1),1),zeros(size(A,1),1),10,1,10,10,1e-1,1e-12,100,@cgCommon,500);k


mprgpctx.lb = PetscBinaryRead("~/devel/permon/examples/Membranes_c");
mprgpctx.abarmult = 2;
mprgpctx.propConst = 1;
mprgpctx.settol = 10*eps;
mprgpctx.infeastol = 0;
mprgpctx.maxeig=6.928472126187e+01;
[x,flg,k] = smalxe(A,b,B,zeros(size(B,1),1),zeros(size(A,1),1),10,1,10,10,1e-1,1e-12,100,@mprgp,100,mprgpctx);k

