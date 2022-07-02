load('nsamp_post_100_3s_null_v22.mat', 'models1')
model=models1{21,1};
model.plot_model_pred(model,subplot(1,1,1),[0,0,0])
standardize_figure(1,[2.5,2])
saveas(gcf,['plots/fig_3D'],'pdf');