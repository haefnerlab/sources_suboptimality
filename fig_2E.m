cc=cbrewer('div','RdYlGn',11);
cc=cc([7,9,11],:);

cc(2,:)=[0,0,0];
cc(3,:)=[1,0,0];

cc1=cbrewer('seq','YlOrRd',9);
cc1=cc1(5:end,:);
cc(3,:)=cc1(1,:);


mu_lik=0.1;
sig_lik=0.3;
mu_prior=-0.3;
sig_prior=0.31;
gam=(sig_prior^2)./((sig_lik^2)+(sig_prior^2));
mu_post=mu_lik*gam+mu_prior*(1-gam);
sig_post=sig_lik*sqrt(gam);
bet=0.45;

scl1=bet/normcdf(mu_prior./sig_prior,0,1);
scl2=(1-bet)/normcdf(-mu_prior./sig_prior,0,1);

plot([-2,2],[0.5,0.5],'k:','linewidth',1.5);
hold on

% plot_hatched([0.5,0.5,1.5,1.5],[0,1-normcdf(0,mu_post,sig_post),1-normcdf(0,mu_post,sig_post),0],2,cc(2,:),'-',45,5,1);
t1=1-normcdf(0,mu_post,sig_post);
t1=scl1*t1./(scl1*t1+scl2*(1-t1));
patch([0.5,0.5,1.5,1.5],[0,t1,t1,0],[0.75,0.75,0.75]);

t2=1-normcdf(0,mu_prior,sig_prior);


plot_hatched([-1.5,-1.5,-0.5,-0.5],[0,t2,t2,0],2,cc(3,:),'-',315,20,1);
plot_hatched([-1.5,-1.5,-0.5,-0.5],[t2,bet,bet,t2],2,cc1(4,:),'-',45,20,1);

plot([-1.5,-1.5,-0.5,-0.5,-1.5],[t2,bet,bet,t2,t2],'color',cc1(4,:),'linewidth',2);



bar([1],[t1],'EdgeColor',cc(2,:),'FaceAlpha',0,'LineWidth',2,'BarWidth',1);
% bar([-1],[t2],'EdgeColor',cc(3,:),'FaceAlpha',0,'LineWidth',2,'BarWidth',1);

xticks([])
xlim([-2,2])
yticks([0,0.5,1]);
ylim([0,1]);
standardize_figure(1,[1.5,1.5])

saveas(gcf,['plots/fig_2E'],'pdf');







% cc=cbrewer('div','RdYlGn',11);
% cc=cc([7,9,11],:);
% 
% 
% mu_lik=0.2;
% sig_lik=0.15;
% mu_prior=-0.25;
% sig_prior=0.25;
% gam=(sig_prior^2)./((sig_lik^2)+(sig_prior^2));
% mu_post=mu_lik*gam+mu_prior*(1-gam);
% sig_post=sig_lik*sqrt(gam);
% 
% plot_hatched([0.5,0.5,1.5,1.5],[0,1-normcdf(0,mu_prior,sig_prior),1-normcdf(0,mu_prior,sig_prior),0],2,cc(3,:),'-',45,5,1);
% hold on
% plot([0.5,0.5,1.5,1.5,0.5],[1-normcdf(0,mu_prior,sig_prior),0.3,0.3,1-normcdf(0,mu_prior,sig_prior),1-normcdf(0,mu_prior,sig_prior)],'r','linewidth',2);
% plot([0,2],[0.5,0.5],'k:','linewidth',1.5);
% % bar([1],[1-normcdf(0,mu_prior,sig_prior)],'EdgeColor',cc(3,:),'FaceAlpha',0,'LineWidth',2,'BarWidth',0.5);
% 
% xticks([])
% xlim([0,2])
% yticks([0,0.5,1]);
% ylim([0,1]);
% standardize_figure(1,[1.25,1.5])
% 
% saveas(gcf,['../../plots/fig1c4'],'pdf');