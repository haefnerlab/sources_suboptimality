x=linspace(0,1,1001);
% cc=cbrewer('div','RdGy',11);
cc=cbrewer('seq','YlGnBu',9);
cc=cc(5:end,:);
cc=cc([1,2,3],:);
% cc=cc([1,2,3],:);
y1=betainc(x,1,1);
y2=betainc(x,5,5);
y3=x>0.5+0.5*(x==0.5);


lr=0;
lb=0.5;
plot(x,lr*lb+(1-lr)*y1,'linewidth',2,'color',cc(3,:))
hold on
plot(x,lr*lb+(1-lr)*y2,'linewidth',2,'color',cc(2,:));
hold on
plot(x,lr*lb+(1-lr)*y3,'linewidth',2,'color',cc(1,:));

% plot([0,0.5],[lr*lb+(1-lr),lr*lb+(1-lr)],'k:');


axis fill
set(gca,'linewidth',2)
box off
set(gca,'TickDir','out')
set(gcf,'color','white')
% xticklabels([]);
% yticklabels([]);
xticks([0,0.5,1]);
yticks([0,0.5,1]);
yticklabels({'0','0.5','1'});
xlim([0,1]);
ylim([0,1]);
xlabel('')
ylabel('')
set(gca,'fontsize',12,'fontweight','normal','fontname','Helvetica Neue')


standardize_figure(1,[1.5,1.5])

% size_in=[3.25,3.25];
% set(gca,'labelfontsizemultiplier',1);
% 
% set(gcf,'PaperUnits','inches')
% set(gcf,'Units','normalized')
% set(gcf,'PaperPosition',[0,0,size_in]);
% set(gcf,'PaperSize',size_in);
% set(gcf,'resize','off')

saveas(gcf,['plots/fig_2F'],'pdf');