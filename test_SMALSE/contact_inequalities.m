function [B] = contact_inequalities(u,points,fracture_matrice,maping,n_elast)

mapingnode=maping(2:2:end)/2;
coods=u(1:n_elast);
points=points';
coods(maping)=coods(maping)+points(:);
cood1=coods(1:2:end);
cood2=coods(2:2:end);

n=length(fracture_matrice);
B=cell(n,1);
for i=1:n
    upn=mapingnode(fracture_matrice{i}.above_nodes);
    downn=mapingnode(fracture_matrice{i}.under_nodes);
    
    
    dlx=cood1(downn(:,1));
    dly=cood2(downn(:,1));
    ulx=cood1(upn(:,1));
    uly=cood2(upn(:,1));
    drx=cood1(downn(:,2));
    dry=cood2(downn(:,2));
    urx=cood1(upn(:,2));
    ury=cood2(upn(:,2));
    
    middleleftx=(ulx+dlx)/2;
    middlelefty=(uly+dly)/2;
    middlerightx=(urx+drx)/2;
    middlerighty=(ury+dry)/2;
    
    normalsx=-(middlerighty-middlelefty);
    normalsy=middlerightx-middleleftx;
    normals_norm=sqrt(normalsx.^2+normalsy.^2);
    normalsx=normalsx./normals_norm;
    normalsy=normalsy./normals_norm;
    
    ljidx=[upn(2:end,1)'*2-1;upn(2:end,1)'*2;downn(2:end,1)'*2-1;downn(2:end,1)'*2];
    tmpleft=1:(length(normalsx)-1);
    liidx=repmat(tmpleft,4,1);
    lvals=[normalsx(2:end)';normalsy(2:end)';-normalsx(2:end)';-normalsy(2:end)'];
    
    rjidx=[upn(1:(end-1),2)'*2-1;upn(1:(end-1),2)'*2;downn(1:(end-1),2)'*2-1;downn(1:(end-1),2)'*2];
    tmpright=(1:(length(normalsx)-1))+tmpleft(end);
    riidx=repmat(tmpright,4,1);
    rvals=[normalsx(1:(end-1))';normalsy(1:(end-1))';-normalsx(1:(end-1))';-normalsy(1:(end-1))'];
    
    B{i}=sparse([liidx(:);riidx(:)],[ljidx(:);rjidx(:)],[lvals(:);rvals(:)],tmpright(end),n_elast); 
end
B=cat(1,B{:});
end