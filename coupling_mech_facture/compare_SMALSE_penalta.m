
kk=zeros(length(D_penalta),1);
for i=1:length(D_penalta)
kk(i)=length(D_penalta{i});
end
ends=[0 ;cumsum(kk)];
figure
hold on

for i=1:length(D_penalta)
%plot((ends(i)+1):ends(i+1),D_penalta{i},'b-')
plot((ends(i)+1):ends(i+1),D_penalta{i}-D{i},'r-')
plot([ends(i+1) ends(i+1)],[-1 1],'k-')
end
ylim([-0.4e-9 0.4e-9])