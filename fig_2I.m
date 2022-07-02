bias = @(alp,mu_sigp,beta,siga_2) normcdf(sqrt(siga_2)*(norminv(beta)-(mu_sigp.*(1-sqrt(1-alp))))./sqrt(alp));
mus=linspace(-2,2,1001);
betas=linspace(-1,1,1001);
[q1,q2]=meshgrid(mus,betas);
alp=0.6;
siga_2=2/3;


y=bias(alp,q1,q2+normcdf(q1),siga_2);
xs=mus;


imagesc(mus,betas,y+0.5)
hold on
p=patch([xs,2,-2],[1-normcdf(xs),1,1],[1,1,1]);
% hatchfill2(p,'speckle','speckledensity',10000);
hatchfill2(p,'cross','hatchcolor',[0.5,0.5,0.5],'hatchdensity',10);
% plot_hatched([xs,2,-2],[1-normcdf(xs),1,1],2,[0,0,1],'-',45,25,1);
hold on
p=patch([xs,-2,-2],[-normcdf(xs),-1,0],[1,1,1]);
% hatchfill2(p,'speckle','speckledensity',10000);
hatchfill2(p,'cross','hatchcolor',[0.5,0.5,0.5],'hatchdensity',10);
% plot_hatched([xs,-2,-2],[-normcdf(xs),-1,0],2,[0,0,1],'-',45,25,1);
% plot([-2,2,2,-2,-2],[-1,-1,1,1,-1],'k','linewidth',4);


set(gcf,'color','white')
set(gca,'Ydir','normal')
axis square
set(gca,'TickDir','out')
set(gca,'linewidth',2)
box off
set(gca,'fontsize',11,'fontname','Helvetica Neue')
xticks([-2,-1,0,1,2])
yticks([-1,-0.5,0,0.5,1])
xticklabels({'-2','','0','','2'})
yticklabels({'-1','','0','','1'})

xlabel('perceptual bias');
ylabel('categorical bias');


c=colorbar;
c.Label.String = 'empirical bias';
c.Label.FontSize=11;
c.Label.FontName='Helvetica';
c.Ticks=[0.5,1,1.5];
c.TickLabels={'0','0.5','1'};
colormap gray

beta_bias=@(alp,mu,bias,siga_2) normcdf(sqrt(alp/siga_2).*norminv(bias)+mu.*(1-sqrt(1-alp)));
cc=cbrewer('seq','YlOrRd',9);
cc=cc(5:end,:);




x=mus;
y=beta_bias(alp,mus,0.1,siga_2);
hold on
plot(x,y-normcdf(x),'color',cc(1,:),'linewidth',2)
x=mus;
y=beta_bias(alp,mus,0.25,siga_2);
plot(x,y-normcdf(x),'color',cc(2,:),'linewidth',2)
x=mus;
y=beta_bias(alp,mus,0.5,siga_2);
plot(x,y-normcdf(x),'color',cc(3,:),'linewidth',2)
x=mus;
y=beta_bias(alp,mus,0.75,siga_2);
plot(x,y-normcdf(x),'color',cc(4,:),'linewidth',2)
x=mus;
y=beta_bias(alp,mus,0.9,siga_2);
plot(x,y-normcdf(x),'color',cc(5,:),'linewidth',2)


xlim([-2,2]);
ylim([-1,1]);








standardize_figure(1,[3,3])

axis square

% size_in=[3,3];
% set(gca,'labelfontsizemultiplier',1);
% 
% set(gcf,'PaperUnits','inches')
% set(gcf,'Units','normalized')
% set(gcf,'PaperPosition',[0,0,size_in]);
% set(gcf,'PaperSize',size_in);
% set(gcf,'resize','off')

saveas(gcf,['plots/fig_2I'],'pdf');
