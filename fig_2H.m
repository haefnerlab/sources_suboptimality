cc=cbrewer('seq','YlGnBu',9);
cc=cc(5:end,:);

xs=linspace(-3,3,101);

% tts=[-1.5,-0.75,0,0.75,1.5];
tts=log10([0.25,0.5,1,2,4]*(norminv(0.75)-norminv(0.25)));
tss=(10.^(tts)./(norminv(0.75)-norminv(0.25)));




plot([-3,3],[0.5,0.5],'k:');
hold on
% plot([0,0],[0,1],'k:');


for i=1:5
ys=normcdf(xs,0,tss(i));
plot(xs,ys,'linewidth',2,'color',cc(i,:));
end



xticks([-3:1.5:3]);
yticks([0:0.25:1]);
xtickangle(0)
ytickangle(0)

standardize_figure(1,[2,2])

saveas(gcf,['plots/fig_2H'],'pdf');