cc=cbrewer('div','PiYG',11);
cc=cc(end-1,:);


siga_2=1;
mu=0.25;

xs=linspace(-3,3,101);
ys=normcdf(xs,mu,sqrt(siga_2));


plot([min(xs),max(xs)],[0.5,0.5],'k:','linewidth',1.5);
hold on
plot([0,0],[0,1],'k:','linewidth',1.5);

plot([min(xs),norminv(0.25)*sqrt(siga_2)+mu],[0.25,0.25],'b:','linewidth',2);
plot([norminv(0.25)*sqrt(siga_2)+mu,norminv(0.25)*sqrt(siga_2)+mu],[0,0.25],'b:','linewidth',2);
plot([min(xs),norminv(0.75)*sqrt(siga_2)+mu],[0.75,0.75],'b:','linewidth',2);
plot([norminv(0.75)*sqrt(siga_2)+mu,norminv(0.75)*sqrt(siga_2)+mu],[0,0.75],'b:','linewidth',2);


sl=0.5./((norminv(0.75)-norminv(0.25))*sqrt(siga_2));
c=0.5-mu*sl;
xmin=norminv(0.025)*sqrt(siga_2)+mu;
xmax=norminv(0.975)*sqrt(siga_2)+mu;




plot(xs,ys,'linewidth',2,'color','k');
plot([xmin,xmax],c+sl*[xmin,xmax],'b','linewidth',1.5);

xlabel('');
ylabel('');



xticks([-1:0.5:1]*max(xs));
yticks([0:0.25:1]);
xticklabels([]);
yticklabels({'0','0.25','0.5','0.75','1'});
xlim([min(xs),max(xs)]);
ylim([0,1]);

set(gcf,'color','white');hold on
set(gca,'linewidth',2)
box off
set(gca,'TickDir','out')

axis fill




set(gca,'fontsize',11,'fontweight','normal','fontname','Helvetica Neue')


standardize_figure(1,[2,2])


saveas(gcf,['plots/fig_1B'],'pdf');