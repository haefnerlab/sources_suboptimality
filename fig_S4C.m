for i=1:18
load(['post_samps_v22/vardur/subj_samps_vardur_bfit_',num2str(i),'.mat']);
ve(i)=get_feve_vardur(subj.phi_theta(subj.params_phi_map,subj),subj,subj.num_eff_params,0);

load(['post_samps_v22/vardur_null/subj_samps_vardur_bfit_',num2str(i),'.mat']);
ve_null(i)=get_feve_vardur(subj.phi_theta(subj.params_phi_map,subj),subj,subj.num_eff_params,0);

load(['post_samps_v22/vardur_null1/subj_samps_vardur_bfit_',num2str(i),'.mat']);
ve_null1(i)=get_feve_vardur(subj.phi_theta(subj.params_phi_map,subj),subj,subj.num_eff_params,0);


i
end

%%

ves=[ve(:),ve_null(:),ve_null1(:)];
figure(1);
for i=1:size(ves,1)
    
    rs(i,:)=(rand(1,3)-0.5);
plot(ves(i,1:3),[1,1.5,2]+0.2*rs(i,:),'.-','color',[0.75,0.75,0.75],'markersize',8,'linewidth',1.5);


hold on
drawnow
end

% 
% 
ids=[2];
cc=cbrewer('div','BrBG',11);
cc=cc([3],:);

for ii=1:numel(ids)

% plot(ves(ids(ii),1:3),[1,1.5,2]+0.2*rs(ids(ii),:),'.-','color',cc(ii,:),'markersize',8,'linewidth',1.5);
plot(ves(ids(ii),1:3),[1,1.5,2]+0.2*rs(ids(ii),:),'.-','color',[0,0,0],'markersize',8,'linewidth',1);
hold on
drawnow
end
% 


plot(mean(ves),[1,1.5,2],'.-','color',[0,0,0],'markersize',12,'linewidth',2);



xticks([0.875,0.9,0.95,1])
yticks([])
set(gca,'TickDir','out')
xticklabels({'0','0.9','','1'});
yticklabels('');
xlim([0.87,1.01])
ylim([0.8,2.2]);

% xlabel({'explainable variance explained','(adj. for # params)'});




% load('subjdata_all_agg_ve.mat')
% ves=[ve(:),ve_null1(:),ve_null(:)];
% plot(mean(ves),[1,1.5,2],'.-','color',[0,0,1],'markersize',12,'linewidth',3);




standardize_figure(1,[3,1.5])
saveas(gcf,['plots/fig_S4C'],'pdf');

