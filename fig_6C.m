% load('subjdata_all_agg_fit_v2.mat', 'lbf_fixed','lbf_fixed_null','lbf_fixed_null1')
load(['post_samps_v22/fixed/subjdata_fixed_agg_v22_fit.mat'], 'lbf_fixed','lbf_fixed_null','lbf_fixed_null1');
t1=logsumexp(lbf_fixed,2)-log(numel(lbf_fixed));
t2=logsumexp(lbf_fixed_null1,2)-log(numel(lbf_fixed_null1));
t3=logsumexp(lbf_fixed_null,2)-log(numel(lbf_fixed_null));
lbfs=[t1,t2,t3];

lb=sort(lbf_fixed);
st=46;
lbf=-1*(logsumexp(-lb(st:end),2)-log(numel(lb(st:end))));
lb=sort(lbf_fixed_null);
lbf1=-1*(logsumexp(-lb(st:end),2)-log(numel(lb(st:end))));
lb=sort(lbf_fixed_null1);
lbf2=-1*(logsumexp(-lb(st:end),2)-log(numel(lb(st:end))));


lbfs=[lbf,lbf1,lbf2];



figure(1);

% for i=1:460
%     plot(log(1+abs(log10(exp(1))*(lbf_fixed(i)-[lbf_fixed_null1(i),lbf_fixed_null(i)]))),[1,1.5],'-','color',[0,0,0,0.05],'markersize',10,'linewidth',1);
%     hold on
%     drawnow
% end
plot(log(1+abs(log10(exp(1))*(lbfs(1)-lbfs(2:3)))),[1,1.5],'.-','color',[0,0,0],'markersize',10,'linewidth',0.75);
hold on

% plot([0,0],[0.9,1.6],'m-.','linewidth',1);
plot(log(1+abs([0.5,0.5])),[0.9,1.6],'k:','linewidth',1.5);
plot(-log(1+abs([0.5,0.5])),[0.9,1.6],'k:','linewidth',1.5);
% plot([1,1],[0.9,1.6],'k-.','linewidth',0.5);
% plot([2,2],[0.9,1.6],'k-.','linewidth',0.5);

tt=[-4:4];
xticks(sign(tt).*log(1+abs(tt)));
yticks([1,1.5])
set(gca,'TickDir','out')
xticklabels({'10^{-4}','','','','10^0','','','','10^4'});
yticklabels('');
xlim([-log(15),log(15)])
ylim([0.9,1.6]);


% xlabel({'\Delta log_{10} Bayes Factor','(approx. inference - .)'});



% 
standardize_figure(1,[3,1.5])
saveas(gcf,['plots/fig_6C'],'pdf');

