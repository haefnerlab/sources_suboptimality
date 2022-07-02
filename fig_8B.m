% load('subjdata_all_agg_fit_v2.mat', 'lbf_vardur','lbf_vardur_null','lbf_vardur_null1')
% 
% load('subjdata_all_agg_ve_vardur.mat')


load(['post_samps_v22/vardur/subjdata_vardur_agg_v22_fit.mat'], 'lbf_vardur','lbf_vardur_null','lbf_vardur_null1');


lb=sort(lbf_vardur);
st=46;
lbf=-1*(logsumexp(-lb(st:end),2)-log(numel(lb(st:end))));
lb=sort(lbf_vardur_null);
lbf1=-1*(logsumexp(-lb(st:end),2)-log(numel(lb(st:end))));
lb=sort(lbf_vardur_null1);
lbf2=-1*(logsumexp(-lb(st:end),2)-log(numel(lb(st:end))));


lbfs=[lbf,lbf1,lbf2];

% 
% % ves=[ve(:),ve_null1(:),ve_null(:)];
% 
% 
% %TODO: Correct this by refitting
% % id_del=((sum(ves<0.9,2)>0));
% 
% 
% t1=logsumexp(lbf_vardur,2)-log(numel(lbf_vardur));
% t2=logsumexp(lbf_vardur_null1,2)-log(numel(lbf_vardur_null1));
% 
% %TODO fix this by refitting
% 
% 
% lbf_vardur_null(lbf_vardur_null==0)=[];
% t3=logsumexp(lbf_vardur_null,2)-log(numel(lbf_vardur_null));
% lbfs=[t1,t2,t3];

figure(1);


plot(log(1+abs(log10(exp(1))*(lbfs(1)-lbfs(2:3)))),[1,1.5],'.-','color',[0,0,0],'markersize',10,'linewidth',0.75);
hold on

plot(log(1+abs([0.5,0.5])),[0.9,1.6],'k:','linewidth',1.5);
plot(-log(1+abs([0.5,0.5])),[0.9,1.6],'k:','linewidth',1.5);
% plot([1,1],[0.9,1.6],'k-.','linewidth',0.5);
% plot([2,2],[0.9,1.6],'k-.','linewidth',0.5);

tt=[-4:4];
xticks(sign(tt).*log(1+abs(tt)));
yticks([1,1.5])
set(gca,'TickDir','out')
xticklabels({'10^{-4}','','10^{-2}','','10^0','','10^2','','10^4'});
yticklabels('');
xlim([-log(1+4.5),log(1+4.5)])
ylim([0.9,1.6]);


% xlabel({'\Delta log_{10} Bayes Factor','(approx. inference - .)'});



% 
standardize_figure(1,[3,1.5])
saveas(gcf,['plots/fig_8B'],'pdf');

