

figure
subplot(2,1,1)
plot(v_grid/ybar,cdf_mat(:,1,ceil(nuG/2),1),'linewidth',2); hold on;
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),ceil(nuG/2),1),'r--','linewidth',2); hold on;
plot(v_grid/ybar,cdf_mat(:,end,ceil(nuG/2),1),'k-.','linewidth',2); hold on;
legend('low A','mean A','high A')
title('Estimated CDF (z^{rs}_ = 1)','fontsize',12)
xlabel('Debt-GDP','fontsize',12)
grid on
subplot(2,1,2)
plot(v_grid/ybar,cdf_mat(:,1,ceil(nuG/2),2),'linewidth',2); hold on;
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),ceil(nuG/2),2),'r--','linewidth',2); hold on;
plot(v_grid/ybar,cdf_mat(:,end,ceil(nuG/2),2),'k-.','linewidth',2); hold on;
legend('low A','mean A','high A')
title('Estimated CDF (z^{rs}_ = 2)','fontsize',12)
xlabel('Debt-GDP','fontsize',12)
grid on

figure
subplot(2,1,1)
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),1,1),'linewidth',2); hold on;
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),ceil(nuG/2),1),'r--','linewidth',2); hold on;
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),end,1),'k-.','linewidth',2); hold on;
grid on
legend('low G','mean G','high G')
title('Estimated CDF','fontsize',12)
xlabel('Debt-GDP','fontsize',12)
subplot(2,1,2)
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),1,2),'linewidth',2); hold on;
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),ceil(nuG/2),2),'r--','linewidth',2); hold on;
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),end,2),'k-.','linewidth',2); hold on;
grid on
legend('low G','mean G','high G')
title('Estimated CDF','fontsize',12)
xlabel('Debt-GDP','fontsize',12)

figure
subplot(2,1,1)
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),ceil(nuG/2),1),'linewidth',2); hold on;
plot(v_grid/ybar,cdf_mat(:,ceil(nuA/2),ceil(nuG/2),2),'k-.','linewidth',2); hold on;
grid on
legend('zrs = 1','zrs = 2')
title('Estimated CDF','fontsize',12)
xlabel('Debt-GDP','fontsize',12)


figure
subplot(2,1,1)
plot(v_grid/ybar,pr_mat(:,1,ceil(nuG/2),1),'linewidth',2); hold on;
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),ceil(nuG/2),1),'r--','linewidth',2); hold on;
plot(v_grid/ybar,pr_mat(:,end,ceil(nuG/2),1),'k-.','linewidth',2); hold on;
legend('low A','mean A','high A')
title('Estimated CDF (z^{rs}_ = 1)','fontsize',12)
xlabel('Debt-GDP','fontsize',12)
grid on
subplot(2,1,2)
plot(v_grid/ybar,pr_mat(:,1,ceil(nuG/2),2),'linewidth',2); hold on;
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),ceil(nuG/2),2),'r--','linewidth',2); hold on;
plot(v_grid/ybar,pr_mat(:,end,ceil(nuG/2),2),'k-.','linewidth',2); hold on;
legend('low A','mean A','high A')
title('Estimated CDF (z^{rs}_ = 2)','fontsize',12)
xlabel('Debt-GDP','fontsize',12)
grid on

figure
subplot(2,1,1)
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),1,1),'linewidth',2); hold on;
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),ceil(nuG/2),1),'r--','linewidth',2); hold on;
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),end,1),'k-.','linewidth',2); hold on;
grid on
legend('low G','mean G','high G')
title('Estimated CDF','fontsize',12)
xlabel('Debt-GDP','fontsize',12)
subplot(2,1,2)
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),1,2),'linewidth',2); hold on;
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),ceil(nuG/2),2),'r--','linewidth',2); hold on;
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),end,2),'k-.','linewidth',2); hold on;
grid on
legend('low G','mean G','high G')
title('Estimated CDF','fontsize',12)
xlabel('Debt-GDP','fontsize',12)

figure
subplot(2,1,1)
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),ceil(nuG/2),1),'linewidth',2); hold on;
plot(v_grid/ybar,pr_mat(:,ceil(nuA/2),ceil(nuG/2),2),'k-.','linewidth',2); hold on;
grid on
legend('zrs = 1','zrs = 2')
title('Estimated CDF','fontsize',12)
xlabel('Debt-GDP','fontsize',12)

% set(gcf,'PaperPositionMode','manual')
% set(gcf,'PaperUnits','inches','PaperPosition',[0 4 7 7])
% % left, bottom, horitonzal length, vertial length   
% print -depsc bstar_dist.eps
% 
% 