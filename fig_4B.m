ids=[20,10,1,17];
for i=1:numel(ids)
    
ii=ids(i);
load(['post_samps_v22/fixed/subj_samps_',num2str(ii),'.mat']);
subj.params_phi_map=subj.params_phi_map_theta;
subj.plot_model_pred(subj,subplot(1,1,1),0);drawnow

xlabel('');
ylabel('');
standardize_figure(1,[2,2])

saveas(gcf,['plots/fig_4B_',num2str(i)],'pdf');
close all
end


