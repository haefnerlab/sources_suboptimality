cc=cbrewer('seq','YlOrRd',9);
cc=cc(5:end,:);

xs=linspace(-3,3,101);

bss=[0.1,0.25,0.5,0.75,0.9];




plot([-3,3],[0.5,0.5],'k:');
hold on
plot([0,0],[0,1],'k:');


for i=1:5
ys=normcdf(xs,-norminv(bss(i)),1);
plot(xs,ys,'linewidth',2,'color',cc(i,:));
end



xticks([-3:1.5:3]);
yticks([0:0.25:1]);

standardize_figure(1,[2,2])

saveas(gcf,['plots/fig_2G'],'pdf');