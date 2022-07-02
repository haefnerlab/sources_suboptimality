nsamps=logspace(0,2,1001);
sigs=logspace(-1,1,1001);
alp=0.6;

thresh = @(alp,sig,nsamp) (norminv(0.75)-norminv(0.25)).*sig.*sqrt(1+(psi(1,0.5*(nsamp+1))./(alp.*psi(1,1))));

[q1,q2]=meshgrid(sigs,nsamps);

y=thresh(alp,q1,q2);
imagesc(log10(sigs),log10(nsamps),log10(y))

set(gcf,'color','white')
set(gca,'Ydir','normal')
axis square
set(gca,'TickDir','out')
set(gca,'linewidth',2)
box off
set(gca,'fontsize',12,'fontname','Helvetica Neue')
xticks(0.5*[-2,-1,0,1,2])
yticks([0,0.5,1,1.5,2])
xticklabels({'0.1','','1','','10'})
yticklabels({'1','','10','','100'})

xlabel('sensory uncertainty');
ylabel('number of samples');


c=colorbar;
c.Label.String = 'empirical threshold';
c.Label.FontSize=11;
c.Label.FontName='Helvetica';
c.Ticks=[-1,0,1];
c.TickLabels={'0.01','1','100'};
colormap gray
caxis([min(log10(sigs)),max(log10(sigs))]);



cc=cbrewer('seq','YlGnBu',9);
cc=cc(5:end,:);


sig_thres = @(alp,thresh,nsamp) thresh./((norminv(0.75)-norminv(0.25)).*sqrt(1+(psi(1,0.5*(nsamp+1))./(alp.*psi(1,1)))));
x=nsamps;
% tts=0.5*[-1.5,-0.75,0,0.75,1.5];
tts=log10([0.25,0.5,1,2,4]*(norminv(0.75)-norminv(0.25)));
for i=1:5
y=sig_thres(alp,10^(tts(i)),nsamps);
hold on
plot(log10(y),log10(x),'color',cc(i,:),'linewidth',2)
end


xlim([min(log10(sigs)),max(log10(sigs))]);
ylim([min(log10(nsamps)),max(log10(nsamps))]);

standardize_figure(1,[3.05,3])

axis square

% size_in=[3,3];
% set(gca,'labelfontsizemultiplier',1);
% 
% set(gcf,'PaperUnits','inches')
% set(gcf,'Units','normalized')
% set(gcf,'PaperPosition',[0,0,size_in]);
% set(gcf,'PaperSize',size_in);
% set(gcf,'resize','off')

saveas(gcf,['plots/fig_2J'],'pdf');