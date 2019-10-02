function [ u0, A, b, freeNode, downEdge, rightEdge, upEdge, leftEdge, M ] = FEM_windows( node, elem, h_elem, bdFlag, k, f, Dirichlet_windows, g_N )
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
M=sparse(N,N);
for i=1:3
    ii(:)=elem(:,i);
    sA=k.*dot(ve(:,:,i),ve(:,:,i),2)./(4*area);
    A=A+sparse(ii,ii,sA,N,N);
    sM=k.*(2*area)/12;
    M=M+sparse(ii,ii,sM,N,N);
end
for i=1:2
    for j=i+1:3
        ii(:)=elem(:,i);
        jj(:)=elem(:,j);
        sA=k.*dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
        A=A+sparse(ii,jj,sA,N,N);
        A=A+sparse(jj,ii,sA,N,N);
        sM=k.*(2*area)/24;
        M=M+sparse(ii,jj,sM,N,N);
        M=M+sparse(jj,ii,sM,N,N);
    end
end

clear ve clear ii jj sA

%% RIGHT HAND SIDE - 3-POINTS QUADRATURE

mid1=(node(elem(:,2),:)+node(elem(:,3),:))/2; % midpoints
mid2=(node(elem(:,1),:)+node(elem(:,3),:))/2; % midpoints
mid3=(node(elem(:,1),:)+node(elem(:,2),:))/2; % midpoints

bt1=area.*(f(mid2)+f(mid3))/6;
bt2=area.*(f(mid1)+f(mid3))/6;
bt3=area.*(f(mid1)+f(mid2))/6;

b=accumarray(elem(:),[bt1;bt2;bt3],[N 1]);

%% BOUNDARY EDGES EXTRACTION
totalEdge=[elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
downEdge=totalEdge(bdFlag(:)==1,:);    down=(bdFlag==1);   t_down=2*ones(h_elem,1);  v_down=zeros(h_elem,1);  i_down=v_down;
rightEdge=totalEdge(bdFlag(:)==2,:);   right=(bdFlag==2); t_right=2*ones(h_elem,1); v_right=zeros(h_elem,1); i_right=v_right;
upEdge=totalEdge(bdFlag(:)==3,:);      up=(bdFlag==3); t_up=2*ones(h_elem,1);    v_up=zeros(h_elem,1);    i_up=v_up;
leftEdge=totalEdge(bdFlag(:)==4,:);    left=(bdFlag==4); t_left=2*ones(h_elem,1);  v_left=zeros(h_elem,1);  i_left=v_left;
values_D=zeros(N,1);

%% Dirichletova okna ------------------------------------------------
for i=1:size(Dirichlet_windows,1)
    a_=Dirichlet_windows(i,2);
    b_=Dirichlet_windows(i,3);
    okna_stred=(node(2:(h_elem+1),2)+node(1:(h_elem),2))/2;
    temp_=(1:h_elem)';
    
    okno=temp_((okna_stred>a_)&(okna_stred<b_));
    v=Dirichlet_windows(i,4);
    switch Dirichlet_windows(i,1)
        case 1	% okno dole
            o=downEdge(okno,:);
            values_D(o(:))=v;
            t_down(okno)=1; v_down(okno)=v; i_down(okno)=i;
        case 2  % okno vpravo
            o=rightEdge(okno,:);
            values_D(o(:))=v;
            t_right(okno)=1; v_right(okno)=v; i_right(okno)=i;
        case 3  % okno nahore
            o=upEdge(okno,:);
            values_D(o(:))=v;
            t_up(okno)=1; v_up(okno)=v; i_up(okno)=i;
        case 4  % okno vlevo
            o=leftEdge(okno,:);
            values_D(o(:))=v;
            t_left(okno)=1; v_left(okno)=v; i_left(okno)=i;
    end
end
bdFlag(down) =t_down;
bdFlag(right)=t_right;
bdFlag(up)   =t_up;
bdFlag(left) =t_left;
% values_D(down)=v_down;
% values_D(diag(nodes2edge(aa(right),bb(right))))=v_right;
% values_D(diag(nodes2edge(aa(up),bb(up))))=v_up;
% values_D(diag(nodes2edge(aa(left),bb(left))))=v_left;
% id_D(diag(nodes2edge(aa(down),bb(down))))=i_down;
% id_D(diag(nodes2edge(aa(right),bb(right))))=i_right;
% id_D(diag(nodes2edge(aa(up),bb(up))))=i_up;
% id_D(diag(nodes2edge(aa(left),bb(left))))=i_left;

Dirichlet=totalEdge(bdFlag(:)==1,:);
Neumann=totalEdge(bdFlag(:)==2,:);
clear totalEdge

%% DIRICHLET BOUNDARY CONDITIONS
isBdNode=false(N,1);
isBdNode(Dirichlet)=true;
%bdNode=find(isBdNode);
freeNode=find(~isBdNode);
u=values_D;%zeros(N,1);
%u(bdNode)=g_D(node(bdNode,:));
b=b-A*u;
%% NEUMANN BOUNDARY CONDITIONS
if (~isempty(Neumann))
    Nve=node(Neumann(:,1),:) - node(Neumann(:,2),:);
    edgeLength=sqrt(sum(Nve.^2,2));
    mid=(node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
    b=b+accumarray([Neumann(:),ones(2*size(Neumann,1),1)], repmat(edgeLength.*g_N(mid)/2,2,1),[N,1]);
end

%% SOLUTION u_h
%clear bdFlag elem mid1 mid2 mid3 area bt1 bt2 bt3 isBdNode Neumann Dirichlet bdNode f g_D g_N N NT i index j
u0 = u;
%u(freeNode)=A(freeNode,freeNode)\b(freeNode);


%% EXTRACT FLOW FROM THE SOLUTION (ze slabe formulace)
% t=A*u;
% t=reshape(t,h_elem+1,h_elem+1);
% for i=1:size(Dirichlet_windows,1)
%     a_=Dirichlet_windows(i,2);
%     b_=Dirichlet_windows(i,3);
%     okna_stred=(node(2:(h_elem+1),2)+node(1:(h_elem),2))/2;
%     temp_=(1:h_elem)';
%     okno=temp_((okna_stred>a_)&(okna_stred<b_));
%     switch Dirichlet_windows(i,1)
%         case 1	% okno dole
%             o=downEdge(okno,:);
%         case 2  % okno vpravo
%             o=rightEdge(okno,:);
%         case 3  % okno nahore
%             o=upEdge(okno,:);
%         case 4  % okno vlevo
%             o=leftEdge(okno,:);
%     end
%     Q(i)=sum((t(unique(o))));
% end
%sum(Q)
end

