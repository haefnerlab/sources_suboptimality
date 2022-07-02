% load('subjdata_all_agg.mat')
% load(['post_samps_v22/fixed/subjdata_fixed_agg_v22_fit_theta.mat']);
% 
% for i=1:460
% model=fixed_model{i};
% model.design_matrix_pred{1}.eps_a_tone=linspace(-3,3,101);
% model.design_matrix_pred{1}.eps_v_right=0*linspace(-3,3,101);
% model.design_matrix_pred{1}.npts=101;
% model.design_matrix_pred{2}.eps_a_tone=linspace(-3,3,101);
% model.design_matrix_pred{2}.eps_v_right=abs(linspace(-3,3,101));
% model.design_matrix_pred{2}.npts=101;
% tmp=model.simulate_prob_resp_pred(model.phi_theta(model.params_phi_map,model),model);
% xs=model.design_matrix_pred{1}.eps_a_tone;
% tmp_cen(i,:)=tmp{1};
% tmp_mat(i,:)=tmp{2};
% i
% end
% 
% save('summary_agg_fixed','tmp_cen','tmp_mat','xs');


%%

load('summary_agg_fixed')
% load(['post_samps_v22/fixed/subjdata_fixed_agg_v22_fit_theta.mat']);
figure(1);
hold on

% plot([-3,3],[0,0],'b:','linewidth',1);
% plot([0,0],[-0.1,0.1],'b-.','linewidth',1);

ci=abs(quantile(tmp_mat-tmp_cen,[0.16,0.84])-repmat(mean(tmp_mat-tmp_cen),2,1));
shaded_errorbar(xs,mean(tmp_mat-tmp_cen),ci([2,1],:),'lineprop','-k','transparent',false);
plot(xs,mean(tmp_mat-tmp_cen),'k','linewidth',2);
cc=cbrewer('div','RdYlGn',11);
cc=cc([7,9,11],:);
plot([-3,3],[0,0],':','color','k','linewidth',2);

xlim([-3,3]);
ylim([-0.125,0.125]);
yticks([-0.1:0.1:0.1]);

% xlabel({'normalized tone position','(dimensionless)'});
% ylabel({'\Delta prob. ''right'' response','(matched-central)'});

standardize_figure(1,[2,1.5])
saveas(gcf,['plots/fig_6B'],'pdf');