load('fixed_agg_summary_v22.mat', 'mus','covs')




a=[-1-(0.02/2):0.02:1+(0.02/2)];
xs=0.5*(a(1:end-1)+a(2:end));
clear posts
for i=1:460
    
    
    
    
post=normpdf(a,mus(i,1)/4,sqrt(covs(2,2,i)/16));



posts(i,:)=post;
% plot(a,post,'color',[0,0,0,0.05]);
% hold on
i
end


%%

plot(a,mean(posts,1),'color',[0,0,0],'linewidth',2);
hold on

plot(xs,normpdf(xs,0,0.25),'k:','linewidth',2)






% plot([0,0],[0,4],'k:','linewidth',1.5);

ylim([0,4]);



xticks([-1:0.5:1])
yticks([0:2:4])

yticklabels([])




% standardize_figure(1,[2,2])
% xlabel('choice bias')
% ylabel('posterior probability')


standardize_figure(1,[2,1.5])
saveas(gcf,['plots/fig_6E'],'pdf');
