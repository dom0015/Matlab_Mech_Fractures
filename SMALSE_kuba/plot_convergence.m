function  plot_convergence(all_iter,errors)
%PLAT_CONVERGENCE Summary of this function goes here
%   Detailed explanation goes here
figure;
plot(errors)
hold on
set(gca,'yScale','log')
mn=min(errors(:));
mx=max(errors(:));

tmp=cumsum(all_iter(:,1)-all_iter(:,3));
% for i=1:size(tmp,1)
%     plot([tmp(i) tmp(i)],[mn mx],'Color',ones(1,3)*0.5)
% end
ylim([mn mx])
xlabel('inner iter')
grid on
box on