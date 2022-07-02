cc=cbrewer('div','PiYG',11);
cc=cc(end-1,:);


siga_2=1;
mu=0.25;

xs=linspace(-3,3,101);
ys=normcdf(xs,mu,sqrt(siga_2));


plot([-1,1]*max(xs),[0.5,0.5],'k:','linewidth',1.5);
hold on
plot([0,0],[0,1],'k:','linewidth',1.5);

plot([-1,0]*max(xs),[normcdf(0,mu,sqrt(siga_2)),normcdf(0,mu,sqrt(siga_2))],'r:','linewidth',2);
plot([mu,mu],[0,0.5],'r:','linewidth',2);


plot(xs,ys,'linewidth',2,'color','k');

xlabel('');
ylabel('');



xticks([-1:0.5:1]*max(xs));
yticks([0:0.25:1]);
xticklabels([]);
yticklabels({'0','0.25','','0.75','1'});



set(gcf,'color','white');hold on
set(gca,'linewidth',2)
box off
set(gca,'TickDir','out')

axis fill




set(gca,'fontsize',11,'fontweight','normal','fontname','Helvetica Neue')


standardize_figure(1,[2,2])


saveas(gcf,['plots/fig_1A'],'pdf');