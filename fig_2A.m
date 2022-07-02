cc=cbrewer('div','RdYlGn',11);
cc=cc([7,9,11],:);
cc(2,:)=[0,0,0];
cc(3,:)=[1,0,0];

mu_lik=0.1;
sig_lik=0.3;
mu_prior=0;
sig_prior=0.31;
gam=(sig_prior^2)./((sig_lik^2)+(sig_prior^2));
mu_post=mu_lik*gam+mu_prior*(1-gam);
sig_post=sig_lik*sqrt(gam);



x=linspace(-1,1,101);
y_prior=normpdf(x,mu_prior,sig_prior);
y_lik=normpdf(x,mu_lik,sig_lik);
y_post=normpdf(x,mu_post,sig_post);









% plot_hatched([x(x>=0),1,0],[y_post(x>=0),0,0],2,cc(2,:),'-',45,5,1);

patch([x(x>=0),1,0],[y_post(x>=0),0,0],[0.75,0.75,0.75]);


hold on
plot_hatched([x(x>=0),1,0],[y_prior(x>=0),0,0],2,cc(3,:),'-',315,20,1);

% plot_hatched(x(x>=0),y_prior(x>=0),2,cc(3,:),'-',315,25,1);


plot(x,y_post,'linewidth',2,'color',cc(2,:));

plot(x,y_prior,'linewidth',2,'color',cc(3,:))
plot(x,zeros(size(x)),'linewidth',2,'color',[0,0,0])


% plot(x,y_lik,'linewidth',2,'color',cc(1,:))
plot([0,0],[0,3],'k:','linewidth',1.5)


xlabel('')
ylabel('')
xticks([-1,0,mu_lik,1])
yticks([0,1,2])
xticklabels({'-1','0','','1'});
ylim([0,2]);
xlim([-1,1]);
set(gca,'TickDir','out')
set(gcf,'color','white')
set(gca,'linewidth',2)
box off
axis fill
set(gca,'fontsize',12,'fontweight','normal','fontname','Helvetica')



standardize_figure(1,[2,1.5])
% size_in=[3,3];
% set(gca,'labelfontsizemultiplier',1);
% 
% set(gcf,'PaperUnits','inches')
% set(gcf,'Units','normalized')
% set(gcf,'PaperPosition',[0,0,size_in]);
% set(gcf,'PaperSize',size_in);
% set(gcf,'resize','off')

saveas(gcf,['plots/fig_2A'],'pdf');




