function [flux_obtain_matrix,pressure_obtain_matrix, u0, A, b, freeNode] = FEM_wells( node, elem, bdFlag, k,boundary_table,neumann_wells)
%FEM solution of -div(k*grad(p))=f
%   node ... coordinates of te nodes
%   elem ... vertices of the elements
%   dbFlag ... types of the sides
%               0 ... non-boundary side
%               1 ... Dirichlet boundary condition
%               2 ... Neumann boundary condition

N=size(node,1);     % number of nodes
NT=size(elem,1);    % number of elements

%% VECTORIZATION (EDGES + AREA)
ve=zeros(NT,2,3);
ve(:,:,1)=node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2)=node(elem(:,1),:)-node(elem(:,3),:);
ve(:,:,3)=node(elem(:,2),:)-node(elem(:,1),:);

area=0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));

%% COORDINATE FORMAT - ASSEMBLING 
% ii=zeros(9*NT,1,'double'); 
% jj=zeros(9*NT,1,'double');
% sA=zeros(9*NT,1);
% index=0;
% for i=1:3
%     for j=1:3
%         ii(index+1:index+NT)=elem(:,i);
%         jj(index+1:index+NT)=elem(:,j);
%         sA(index+1:index+NT)=k.*dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
%         index=index+NT;
%     end
% end
% clear ve
% A=sparse(ii,jj,sA,N,N);
% clear ii jj sA

%% COORDINATE FORMAT - ASSEMBLING PER PARTES
% ii=zeros(NT,1,'double'); 
% jj=zeros(NT,1,'double');
% A=sparse(N,N);
% for i=1:3
%     for j=1:3
%         ii(:)=elem(:,i);
%         jj(:)=elem(:,j);
%         sA=k.*dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
%         A=A+sparse(ii,jj,sA,N,N);
%     end
% end
% clear ve ii jj sA

%% USING THE SYMMETRY OF THE MATRIX I
% ii=zeros(9*NT,1,'double'); 
% jj=zeros(9*NT,1,'double');
% sA=zeros(9*NT,1);
% index=0;
% for i=1:3
%     ii(index+1:index+NT)=elem(:,i);
%     jj(index+1:index+NT)=elem(:,i);
%     sA(index+1:index+NT)=k.*dot(ve(:,:,i),ve(:,:,i),2)./(4*area);
%     index=index+NT;
%     for j=i+1:1:3
%         tmp=k.*dot(ve(:,:,i),ve(:,:,j),2)./(4*area);      
%         ii(index+1:index+2*NT)=[elem(:,i); elem(:,j)];
%         jj(index+1:index+2*NT)=[elem(:,j); elem(:,i)];
%         sA(index+1:index+2*NT)=[tmp; tmp];
%         index=index+2*NT;
%     end
% end
% clear ve
% A=sparse(ii,jj,sA,N,N);
% clear ii jj sA

%% USING THE SYMMETRY OF THE MATRIX II
% ii=zeros(6*NT,1,'double'); 
% jj=zeros(6*NT,1,'double');
% sA=zeros(6*NT,1);
% index=0;
% for i=1:3
%     ii(index+1:index+NT)=elem(:,i);
%     jj(index+1:index+NT)=elem(:,i);
%     sA(index+1:index+NT)=k.*dot(ve(:,:,i),ve(:,:,i),2)./(4*area);
%     index=index+NT;
% end
% for i=1:2
%     for j=i+1:3
%         ii(index+1:index+NT)=elem(:,i);
%         jj(index+1:index+NT)=elem(:,j);
%         sA(index+1:index+NT)=k.*dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
%         index=index+NT;
%     end
% end
% clear ve
% A=sparse([ii; jj(3*NT+1:end)],[jj; ii(3*NT+1:end)],[sA; sA(3*NT+1:end)],N,N);
% clear ii jj sA

%% USING THE SYMMETRY + PER PARTES
ii=zeros(NT,1,'double'); 
jj=zeros(NT,1,'double');
A=sparse(N,N);
for i=1:3
    ii(:)=elem(:,i);
    sA=k.*dot(ve(:,:,i),ve(:,:,i),2)./(4*area);
    A=A+sparse(ii,ii,sA,N,N);
end
for i=1:2
    for j=i+1:3
        ii(:)=elem(:,i);
        jj(:)=elem(:,j);
        sA=k.*dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
        A=A+sparse(ii,jj,sA,N,N);
        A=A+sparse(jj,ii,sA,N,N);
    end
end

clear ve clear ii jj sA

%% RIGHT HAND SIDE - 3-POINTS QUADRATURE
% tic;
% mid1=(node(elem(:,2),:)+node(elem(:,3),:))/2; % midpoints
% mid2=(node(elem(:,1),:)+node(elem(:,3),:))/2; % midpoints
% mid3=(node(elem(:,1),:)+node(elem(:,2),:))/2; % midpoints
% 
% bt1=area.*(f(mid2)+f(mid3))/6;
% bt2=area.*(f(mid1)+f(mid3))/6;
% bt3=area.*(f(mid1)+f(mid2))/6;
% 
% b=accumarray(elem(:),[bt1;bt2;bt3],[N 1]);
% toc
% whos A
%% BOUNDARY EDGES EXTRACTION
totalEdge=[elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
Dirichlet=totalEdge(bdFlag(:)>0,:);
Neumann=totalEdge(bdFlag(:)==-1,:);


%% DIRICHLET BOUNDARY CONDITIONS
isBdNode=false(N,1);
isBdNode(Dirichlet)=true;
bdNode=find(isBdNode);
freeNode=find(~isBdNode);
u=zeros(N,1);
flux_obtain_matrix=zeros(size(boundary_table,1),N);
pressure_obtain_matrix=zeros(size(neumann_wells,1),N);

for i=1:size(boundary_table,1)
    Dirichlet_=totalEdge(bdFlag(:)==boundary_table(i,1),:);
    isBdNode_=false(N,1);
    isBdNode_(Dirichlet_)=true;
    flux_obtain_matrix(i,:)=double(isBdNode_)'*A;
    bdNode_=find(isBdNode_);
    u(bdNode_)=boundary_table(i,2);
end

for i=1:size(neumann_wells,1)
    Dirichlet_=totalEdge(bdFlag(:)==neumann_wells(i),:);
    isBdNode_=false(N,1);
    isBdNode_(Dirichlet_)=true;
    pressure_obtain_matrix(i,:)=double(isBdNode_)'/sum(double(isBdNode_));
end

b=-A*u;
%% NEUMANN BOUNDARY CONDITIONS
% if (~isempty(Neumann))
%     Nve=node(Neumann(:,1),:) - node(Neumann(:,2),:);
%     edgeLength=sqrt(sum(Nve.^2,2));
%     mid=(node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
%     b=b+accumarray([Neumann(:),ones(2*size(Neumann,1),1)], repmat(edgeLength.*g_N(mid)/2,2,1),[N,1]);
% end

%% SOLUTION u_h
clear bdFlag elem mid1 mid2 mid3 node area bt1 bt2 bt3 k isBdNode Neumann Dirichlet bdNode f g_D g_N N NT i index j
u0=u;
%u(freeNode)=A(freeNode,freeNode)\b(freeNode);

end

