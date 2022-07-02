for i=1:18
load(['post_samps_v22/vardur/subj_samps_vardur_bfit_',num2str(i),'.mat']);


ss1=subj.phi_theta(get_samples_burnin_thin(subj.post_samps,100,0,46,size(subj.mus{1},1)),subj);
for jj=1:3
kpa_2=ss1(:,5+jj);
kpv_2=ss1(:,8+jj);
sig_ap_2=ss1(:,13).^2;
sig_a_2{i}(:,jj)=sig_ap_2./kpa_2;
sig_v_2{i}(:,jj)=sig_ap_2.*kpv_2;
end
i
end
%%
for i=1:18

mus(i,:)=mean(sqrt(sig_a_2{i}));
ses(i,:)=std(sqrt(sig_a_2{i}));
i
end
wts=1./(repmat(var(mus),18,1)+ses.^2);
wts=wts./repmat(sum(wts),18,1);
mu_eff=sum(wts.*mus);
se_eff=sqrt(1./sum(1./(repmat(var(mus),18,1)+ses.^2)));
for i=1:18
cis=abs(quantile((sqrt(sig_a_2{i})),[0.16,0.84])-mean((sqrt(sig_a_2{i}))));
lb(i,:)=cis(1,:);
ub(i,:)=cis(2,:);
end


for i=1:18
errorbar([1,2,3]+0.1*(rand(1,3)-0.5),mus(i,:),lb(i,:),ub(i,:),'.-','markersize',5,'color',[0.75,0.75,0.75],'linewidth',1,'capsize',0);
hold on
drawnow
end
set(gca,'YScale','log')
errorbar([1:3],mu_eff,se_eff,'k.-','markersize',10,'capsize',0,'linewidth',2)
[h,p]=ttest(mus(:,1),mus(:,2),'tail','right')
p<0.001
[h,p]=ttest(mus(:,2),mus(:,3),'tail','right')
p<0.001

[h,p]=signtest(mus(:,1),mus(:,2))
[h,p]=signtest(mus(:,2),mus(:,3))
normcdf(mu_eff(1)-mu_eff(2))
[h,p]=ttest(mus(:,1),mus(:,3),'tail','right')
yticks([1e-4,1e-3,1e-2,1e-1,1e0])
set(gca,'YMinorTick','off')
yticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'})
xticks([1,2,3])
xlabel('duration')
xticklabels({'100','300','1000'})
xlabel('duration (ms)')
ylabel('sensory uncertainty')
yticklabels({'0.01','0.1','1'})
yticks([1e-2,1e-1,1e0])
xlabel('')
ylabel('')
standardize_figure(1,[2,1.5])
xticks([1,2,3])
xlim([0.5,3.5])
saveas(gcf,['plots/fig_7C'],'pdf');











% set(gca,'YScale','log')
% errorbar([1:3],mu_eff,se_eff,'k.-','markersize',10,'capsize',0,'linewidth',2)
% [h,p]=ttest(mus(:,1),mus(:,2),'tail','right')
% p<0.001
% [h,p]=ttest(mus(:,2),mus(:,3),'tail','right')
% p<0.001
% normcdf(mu_eff(1)-mu_eff(2))
% [h,p]=ttest(mus(:,1),mus(:,3),'tail','right')
% yticks([1e-4,1e-3,1e-2,1e-1,1e0])
% set(gca,'YMinorTick','off')
% yticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'})
% xticks([1,2,3])
% xlabel('duration')
% xticklabels({'100','300','1000'})
% xlabel('duration (ms)')
% standardize_figure(1,[2,2])
% ylabel('sensory uncertainty')
% saveas(gcf,['../../plots/sig_dur'],'pdf');
% yticklabels({'0.01','0.1','1'})
% yticks([1e-2,1e-1,1e0])
% saveas(gcf,['../../plots/sig_dur'],'pdf');