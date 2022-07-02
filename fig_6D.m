% load('post_samps_v22/fixed/subjdata_fixed_agg_v22_fit.mat')
% for i=1:460
% h=num_hess(fixed_model{i});
% covs(:,:,i)=0.5*(inv(h)+inv(h)');
% mus(i,:)=fixed_model{i}.params_phi_map;
% i
% end

load('fixed_agg_summary_v22.mat', 'mus','covs')




xs=[1:0.1:100];a=0.5*(xs(1:end-1)+xs(2:end));
clear posts
for i=1:460
post=get_post_nsamp_cov(mus(i,2),covs(2,2,i),xs);
posts(i,:)=post/mean(diff(xs));
% plot(a,post,'color',[0,0,0,0.05]);
% hold on
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% drawnow
i
end


%%

plot(a,mean(posts,1),'color',[0,0,0],'linewidth',2);
hold on

% plot([1,100],[0.01,0.01],'m-.','linewidth',2);


dd=max(xs)-min(xs);

xticks([1,25,50,75,100])
yticks([0.5,1:5]*(1/dd))
ylim([0.1,50]*(1/dd));
xlim([1,100]);
xticklabels({'1','25','','','100'});





% yticks([0.1,0.2,0.5,1:5,10,25,50]*(1/dd));
% yticklabels({'0.1','0.2','0.5','1','','','','5','10','25','50'});


yticks([1,10,25,50]*(1/dd));
ylim([0,50]*(1/dd));

yticklabels(round([1,10,25,50]*(1/dd),2));


set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'YScale','linear')
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')


standardize_figure(1,[2,1.5])
% xlabel('num. of samples')
% ylabel('likelihood fn.')


% standardize_figure(1,[2,2])
saveas(gcf,['plots/fig_6D'],'pdf');
