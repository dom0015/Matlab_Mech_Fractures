function [frac_bdflag,frac_normx,frac_normy,fracture_elem_map] = map_fract2elem(points,elem,fracture_matrice)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n=max(elem(:));
m=size(elem,1);
elem_matrix1=sparse([elem(:,1);elem(:,2)],[1:m,1:m],ones(2*m,1),n,m);
elem_matrix2=sparse([elem(:,2);elem(:,3)],[1:m,1:m],ones(2*m,1),n,m);
elem_matrix3=sparse([elem(:,3);elem(:,1)],[1:m,1:m],ones(2*m,1),n,m);
frac_bdflag=elem*0;
frac_normx=elem*0;
frac_normy=elem*0;

fracture_elem_map=cell(length(fracture_matrice),1);
for i=1:length(fracture_matrice)
    tmp=fracture_matrice{i};
    
    tmp_norm=points(tmp.above_nodes(:,2),:)-points(tmp.above_nodes(:,1),:);
    tmp_norm=tmp_norm./repmat(sqrt(tmp_norm(:,1).^2+tmp_norm(:,2).^2),1,2);
    
    norm_up=tmp_norm*[0 1;-1 0];
    norm_down=tmp_norm*[0 -1;1 0];
    
    mm=length(tmp.above_nodes);
    frac_spmat=sparse([1:mm,1:mm],tmp.above_nodes(:),ones(2*mm,1),mm,n);
    
    [i1,j1,~]=find(frac_spmat*elem_matrix1>1);
    [i2,j2,~]=find(frac_spmat*elem_matrix2>1);
    [i3,j3,~]=find(frac_spmat*elem_matrix3>1);
    ii=[i1;i2;i3];
    jj=[j1;j2;j3];
    kk=[1+0*j1;2+0*j2;3+0*j3];
    [ii,idx]=sort(ii);
    jj=jj(idx);
    kk=kk(idx);
    
    frac_bdflag((kk-1)*m+jj)=-i;
    frac_normx((kk-1)*m+jj)=norm_up(:,1);
    frac_normy((kk-1)*m+jj)=norm_up(:,2);
    
    tmp_frac_elem_map.upi=jj;
    tmp_frac_elem_map.upj=kk;
    
    frac_spmat=sparse([1:mm,1:mm],tmp.under_nodes(:),ones(2*mm,1),mm,n);
    
    [i1,j1,~]=find(frac_spmat*elem_matrix1>1);
    [i2,j2,~]=find(frac_spmat*elem_matrix2>1);
    [i3,j3,~]=find(frac_spmat*elem_matrix3>1);
    ii=[i1;i2;i3];
    jj=[j1;j2;j3];
    kk=[1+0*j1;2+0*j2;3+0*j3];
    [ii,idx]=sort(ii);
    jj=jj(idx);
    kk=kk(idx);
    
    frac_bdflag((kk-1)*m+jj)=-i-0.5;
    frac_normx((kk-1)*m+jj)=norm_down(:,1);
    frac_normy((kk-1)*m+jj)=norm_down(:,2);
    tmp_frac_elem_map.downi=jj;
    tmp_frac_elem_map.downj=kk;
    fracture_elem_map{i}=tmp_frac_elem_map;
end

end

