% load subjdata_all_lbf.mat
load('lbf_vardur_v22','lbf','lbf1','lbf2');

lbf_vardur=lbf;
lbf_vardur_null=lbf1;
lbf_vardur_null1=lbf2;

lbfs=[lbf_vardur(:),lbf_vardur_null(:),lbf_vardur_null1(:)];

figure(1);
for i=1:size(lbfs,1)
    ttmp=log10(exp(1))*(lbfs(i,1)-lbfs(i,2:3));
plot(sign(ttmp).*log(1+abs(ttmp)),[1,1.5],'.-','color',[0.75,0.75,0.75],'markersize',10,'linewidth',0.75);
hold on
drawnow
end


ids=[2];
cc=cbrewer('div','BrBG',11);
cc=cc([3],:);

for ii=1:numel(ids)

% plot(log10(exp(1))*(lbfs(ids(ii),1)-lbfs(ids(ii),2:3)),[1,1.5],'.-','color',cc(ii,:),'markersize',10,'linewidth',1.5);
ttmp=log10(exp(1))*(lbfs(ids(ii),1)-lbfs(ids(ii),2:3));
plot(sign(ttmp).*log(1+abs(ttmp)),[1,1.5],'.-','color',[0,0,0],'markersize',10,'linewidth',0.75);
hold on
drawnow
end



% load('subjdata_all_agg_fit.mat', 'lbf_vardur','lbf_vardur_null','lbf_vardur_null1')
% t1=logsumexp(lbf_vardur,2)-log(numel(lbf_vardur));
% t2=logsumexp(lbf_vardur_null1,2)-log(numel(lbf_vardur_null1));
% t3=logsumexp(lbf_vardur_null,2)-log(numel(lbf_vardur_null));
% lbfs=[t1,t2,t3];
% 
% 
% plot(log10(exp(1))*(lbfs(1)-lbfs(2:3)),[1,1.5],'.-','color',[0,0,1],'markersize',10,'linewidth',3);


% plot([0,0],[0.9,1.6],'m-.','linewidth',1);
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
saveas(gcf,['plots/fig_7B'],'pdf');

