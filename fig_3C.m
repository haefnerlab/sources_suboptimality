cc=cbrewer('div','RdYlGn',11);
cc=cc(10,:);
cc=[0,0,1];
cc=cbrewer('div','PiYG',11);
cols=cc([10,2],:);
cc=cols(1,:);
% yyaxis left


set(gca,'YColor',cc)
x1=linspace(-4,0,1001);
x2=linspace(0,4,1001);
plot(x1,normpdf(x1,1,0.75),'linewidth',2,'color',cc)
hold on
plot_hatched(x2,normpdf(x2,1,0.75),2,cc,'-',45,10,1);
% plot(x2,normpdf(x2,1,1),'linewidth',2,'color',cc)
cc=cbrewer('div','RdYlGn',11);
cc=cc([2,3,10],:);
cc=cols(2,:);
% yyaxis right


% set(gca,'YColor',cc(1,:))


pd=normpdf(1,1,sqrt(10/9))./(normpdf(1,1,sqrt(10/9))+normpdf(1,-1,sqrt(10/9)));
plot(x1,(1-pd)*normpdf(x1,-1,0.1),'linewidth',2,'color',cc(1,:),'linestyle','-');
plot(x2,pd*normpdf(x2,1,0.1),'linewidth',2,'color',cc(1,:),'linestyle','-');
hold on
plot_hatched(x2,pd*normpdf(x2,1,0.1),2,cc,'-',315,25,1);
yticks([]);

ylim([0,4]);


% stem(1,normpdf(1,1,1),'linewidth',2,'color',cc(1,:),'marker','^','markersize',6,'markerfacecolor','auto');
% stem(-1,normpdf(-1,1,1),'linewidth',2,'color',cc(1,:),'marker','^','markersize',6,'markerfacecolor','auto');

% yyaxis left
% plot([-1,0],normpdf(-1,1,1)*[1,1],':','color',cc(1,:),'linewidth',2);
% plot([0,1],normpdf(1,1,1)*[1,1],':','color',cc(1,:),'linewidth',2);

set(gcf,'color','white');hold on
set(gca,'linewidth',2)
box off
set(gca,'TickDir','out')

axis square


xticks([-1,0,1])
xticklabels({'','0',''})
% text(1,0.05,'o_a','fontsize',12,'fontname','Helvetica Neue');
yticks([]);
ylabel('posterior probability');




% xlabel('inferred tone position');

xlim([-5,5]);
ylim([0,3.5]);
plot([0,0],[0,3.5],'k:','linewidth',1.5);

% set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')

set(gca,'fontsize',11,'fontweight','normal','fontname','Helvetica Neue')
axis fill
standardize_figure(1,[4,4]);
% size_in=[4,4];
% set(gca,'labelfontsizemultiplier',1);
% 
% set(gcf,'PaperUnits','inches')
% set(gcf,'Units','normalized')
% set(gcf,'PaperPosition',[0,0,size_in]);
% set(gcf,'PaperSize',size_in);
% set(gcf,'resize','off')

saveas(gcf,['plots/fig_3C'],'pdf');