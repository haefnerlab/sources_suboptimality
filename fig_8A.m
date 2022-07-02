% %
% 
% load('subjdata_all_agg_fit_v2.mat', 'vardur_model')

% load(['post_samps_v22/vardur/subjdata_vardur_agg_v22_fit_theta.mat'],'vardur_model');


% 
% for i=1:460
%     model=vardur_model{i};
%     for k=1:3
%         model.design_matrix_pred{k}.eps_a_tone=linspace(-2,2,101);
%         model.design_matrix_pred{k}.eps_v_right=0*model.design_matrix_pred{k}.eps_a_tone;
%         model.design_matrix_pred{3+k}.eps_a_tone=linspace(-2,2,101);
%         model.design_matrix_pred{3+k}.eps_v_right=abs(model.design_matrix_pred{3+k}.eps_a_tone);
%         xs=model.design_matrix_pred{k}.eps_a_tone;
%         tmp=model.simulate_prob_resp_pred(model.phi_theta(model.params_phi_map,model),model);
%         bb1{k}(i,:)=tmp{k};
%         bb2{k}(i,:)=tmp{3+k};
%     end
%     i
% end
% 
% save('summary_agg_vardur1','bb1','bb2','xs');





%%

load('summary_agg_vardur1');
cc=cbrewer('qual','Dark2',6);
cc=cc(4:6,:);

cc1=[0,0,0];


for kk=1:3
figure(kk);
hold on


% plot([0,0],[-0.05,0.05],'b-.','linewidth',1);
hold on

% ci=abs(quantile(bb2{kk}-bb1{kk},[0.16,0.84])-repmat(mean(bb2{kk}-bb1{kk}),2,1));
% shaded_errorbar(xs,mean(bb2{kk}-bb1{kk}),ci([2,1],:),'lineprops',{'-','color',cc(kk,:)},'transparent',false);

ci=abs(quantile(bb2{kk}-bb1{kk},[0.16,0.84])-repmat(mean(bb2{kk}-bb1{kk}),2,1));
shaded_errorbar(xs,mean(bb2{kk}-bb1{kk}),ci([2,1],:),'lineprops',{'-','color',cc1},'transparent',false);
xlim([-2,2]);
ylim([-0.06,0.06]);
yticks([-0.05:0.025:0.05]);
yticklabels({'-0.05','','0','','0.05'});
end
for kk=1:3
    
figure(kk);
    
% plot(xs,mean(bb2{kk}-bb1{kk}),'color',cc(kk,:),'linewidth',2);

plot(xs,mean(bb2{kk}-bb1{kk}),'color',cc1,'linewidth',2);
cc=cbrewer('div','RdYlGn',11);
cc=cc([7,9,11],:);
plot([-2,2],[0,0],':','linewidth',2,'color','k');
% if kk==3
% xlabel({'normalized tone position','(dimensionless)'});
% end
% if kk==2
% ylabel({'\Delta prob. ''right'' response','(matched-central)'});
% end
standardize_figure(kk,[2,1.5])
saveas(kk,['plots/fig_8A_',num2str(kk)],'pdf');
end