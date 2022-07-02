ii=20;
load(['post_samps_v22/fixed/subj_samps_',num2str(ii),'.mat']);
ps0=subj.phi_theta(subj.params_phi_map_theta,subj);
ps0=[-0.1,20,0.8,1.5,0.01,-0.6,0.1,0.05,0.3,1,1];

% ps0(4)=0.1;


ps1=[0,1,ps0(3:5),0,ps0(7:end)];
model=subj;
% model.exact_inference=1;
model.params_phi=model.theta_phi(ps1,model);
model=model.simulate_data(model,1e4);
model=rmfield(model,'params_phi_map');
model.plot_model_pred(model,subplot(1,1,1),[0,0,0]);

xlabel('');
ylabel('');
standardize_figure(1,[2,2])
saveas(1,['plots/fig_4A_1'],'pdf');


figure;
% close all

ps1=[0,ps0(2:5),ps0(6),ps0(7:end)];
model=subj;
model.exact_inference=1;
model.params_phi=model.theta_phi(ps1,model);
model=model.simulate_data(model,1e4);
model=rmfield(model,'params_phi_map');
model.plot_model_pred(model,subplot(1,1,1),[0,0,0]);

xlabel('');
ylabel('');
standardize_figure(2,[2,2])
saveas(2,['plots/fig_4A_2'],'pdf');
figure;
% close all
ps1=[ps0(1),ps0(2:5),0,ps0(7:end)];
model=subj;
model.exact_inference=1;
model.params_phi=model.theta_phi(ps1,model);
model=model.simulate_data(model,1e4);
model=rmfield(model,'params_phi_map');
model.plot_model_pred(model,subplot(1,1,1),[0,0,0]);

xlabel('');
ylabel('');
standardize_figure(3,[2,2])
saveas(3,['plots/fig_4A_3'],'pdf');
figure;
% close all

figure(4)
ps0=[-0.025,20,0.8,1.5,0.01,0.6,0.1,0.05,0.3,1,1];
gam=ps0(4)./(1+ps0(4));

ps1=[ps0(1),1,ps0(3:end)];
% ps1=[normcdf(ps0(6)*(1-sqrt(1-gam)))-normcdf(ps0(6)),1,ps0(3:end)];
% ps1=[0.5-normcdf(ps0(6)),1,ps0(3:end)];
model=subj;
model.exact_inference=1;
model.params_phi=model.theta_phi(ps1,model);
model=model.simulate_data(model,1e4);
model=rmfield(model,'params_phi_map');
model.plot_model_pred(model,subplot(1,1,1),[0,0,0]);



hold on


gam=ps0(4)./(1+ps0(4));
ps1=[normcdf(ps0(6)*(1-sqrt(1-gam)))-normcdf(ps0(6)),1,0,ps0(4:end)];
% ps1=[0.5-normcdf(ps0(6)),1,0,ps0(4:end)];
model=subj;
model.exact_inference=1;
model.params_phi=model.theta_phi(ps1,model);
tmp=model.simulate_prob_resp_pred(model.phi_theta(model.params_phi,model),model);
% plot(model.design_matrix_pred{1}.eps_a_tone,tmp{1},'k-','linewidth',1.5)





xlabel('');
ylabel('');
standardize_figure(4,[2,2])
saveas(4,['plots/fig_4A_4'],'pdf');




% figure;
% close all
% 
% ps1=ps0;
% model=subj;
% % model.exact_inference=1;
% model.params_phi=model.theta_phi(ps1,model);
% model=model.simulate_data(model,1e4);
% model=rmfield(model,'params_phi_map');
% model.plot_model_pred(model,subplot(1,1,1),[0,0,0]);
% 
% xlabel('');
% ylabel('');
% standardize_figure(5,[2,2])
% saveas(5,['../../plots/fig2h_5'],'pdf');



