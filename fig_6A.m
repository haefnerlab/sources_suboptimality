% load('subjdata_all_agg.mat')
% load('post_samps_v22/fixed/subjdata_fixed_agg_v22_fit.mat');
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

cols=[0,0,1;1,0,0];
cc=cbrewer('div','PiYG',11);
cols=cc([10,2],:);

figure(1);
hold on

% plot([-3,3],[0.5,0.5],'b-.','linewidth',1);
% plot([0,0],[0,1],'b-.','linewidth',1);

ci=abs(quantile(tmp_cen,[0.16,0.84])-repmat(mean(tmp_cen),2,1));
shaded_errorbar(xs,mean(tmp_cen),ci([2,1],:),'lineprop',{'-','color',cols(1,:)},'transparent',false);


ci=abs(quantile(tmp_mat,[0.16,0.84])-repmat(mean(tmp_mat),2,1));
shaded_errorbar(xs,mean(tmp_mat),ci([2,1],:),'lineprop',{'-','color',cols(2,:)},'transparent',false);


ys=normcdf(xs,0,1);
cc=cbrewer('div','RdYlGn',11);
cc=cc([7,9,11],:);
plot(xs,ys,':','color','k','linewidth',2)


plot(xs,mean(tmp_cen),'color',cols(1,:),'linewidth',2);
plot(xs,mean(tmp_mat),'color',cols(2,:),'linewidth',2);






% xlabel({'normalized tone position','(dimensionless)'});
% ylabel('prob. of ''right'' response');

standardize_figure(1,[2,1.5])
saveas(gcf,['plots/fig_6A'],'pdf');



