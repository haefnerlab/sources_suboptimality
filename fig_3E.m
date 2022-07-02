cc=cbrewer('div','RdBu',11);
cc=cc([2,10],:);


load('nsamp_post_100_3s_null_v22.mat','lbf_model1','lbf_model2')




xs=10*logspace(1,4,21);


plot(xs,smooth(mean(lbf_model1-lbf_model2>(0.5/log10(exp(1))),2),10),'color','k','linewidth',2);
hold on

% load('nsamp_post_100_3s_null3_lbf.mat','lbf_model','lbf_model1');

load('nsamp_post_100_3s_null_v22.mat','lbf_model3')

xs=10*logspace(1,4,21);
plot(xs,smooth(mean(lbf_model1-lbf_model3>(0.5/log10(exp(1))),2),10),'--','color','k','linewidth',2);

% plot(xs,smooth(mean(lbf_model1-lbf_model>0,2),3),'color',cc(1,:),'linewidth',2)

% plot(xs,smooth(mean(lbf_model1-lbf_model2>(0.5/log10(exp(1))),2),5),'color',cc(2,:),'linewidth',2);

set(gca,'XScale','log');
xticks([100,1e3,1e4,1e5]);
xticklabels({'10^2','10^3','10^4','10^5'});
set(gca,'XMinorTick','off')
ylim([0,1]);
yticks([0:0.25:1]);

xlabel('');
ylabel({'probability for substantial','evidence in dissociating'});


standardize_figure(1,[3,2])

saveas(gcf,['plots/fig_3E'],'pdf');

