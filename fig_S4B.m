
% load(['post_samps_v22/fixed/subjdata_fixed_agg_v22_fit.mat']);
% for i=1:460
% model=fixed_model{i};
% ve(i)=get_feve(model.phi_theta(model.params_phi_map,model),model,model.num_eff_params,0);
% i
% end
% for i=1:460
% model=fixed_model_null{i};
% ve_null(i)=get_feve(model.phi_theta(model.params_phi_map,model),model,model.num_eff_params,0);
% i
% end
% for i=1:460
% model=fixed_model_null1{i};
% ve_null1(i)=get_feve(model.phi_theta(model.params_phi_map,model),model,model.num_eff_params,0);
% i
% end
% save('subjdata_all_agg_ve_v2','ve','ve_null','ve_null1');


load('subjdata_all_agg_ve_v2.mat')
ves=[ve(:),ve_null1(:),ve_null(:)];


for i=1:size(ves,1)
    
    rs(i,:)=(rand(1,3)-0.5);
plot(ves(i,1:3),[1,1.5,2]+0.2*rs(i,:),'-','color',[0,0,0,0.025],'markersize',8,'linewidth',0.5);


hold on
drawnow
end









plot(mean(ves),[1,1.5,2],'.-','color',[0,0,0],'markersize',12,'linewidth',2);

xticks([0.875,0.9,0.95,1])
yticks([1,1.5,2])
set(gca,'TickDir','out')
xticklabels({'0','0.9','','1'});
yticklabels('');
xlim([0.87,1.01])
ylim([0.8,2.2]);

% xlabel({'explainable variance explained','(adj. for # params)'});









standardize_figure(1,[2,1.5])
saveas(gcf,['plots/fig_S4B'],'pdf');

