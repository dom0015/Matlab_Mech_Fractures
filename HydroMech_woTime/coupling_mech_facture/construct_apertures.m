function [D] = construct_apertures(Bu,fracture_matrice)

n=length(fracture_matrice);
D=cell(n,1);
startidx=1;
endidx=0;
for i=1:n
    upn=fracture_matrice{i}.above_nodes;
    idx_same=upn(1:end-1,2)==upn(2:end,1);
    endidx=endidx+length(idx_same)+sum(idx_same==0);
    
    tmpbu=Bu(startidx:endidx);
    tmp_apertures=zeros(length(idx_same)+1,1);
    tmp_apertures(1:end-1)=tmpbu(1:length(idx_same));
    idx_d=find(~idx_same);
    tmpbu(idx_d-1)=tmpbu((length(idx_same)+1):end);
    tmp_apertures(2:end)=tmp_apertures(2:end)+tmpbu(1:length(idx_same));
    tmp_apertures=tmp_apertures/2;
    D{i}=tmp_apertures;
    startidx=startidx+(endidx-startidx)+1;
end
end