% load('post_samps_v22/vardur/subjdata_vardur_agg_v22_fit.mat')
% parfor i=1:460
% h=num_hess(vardur_model{i});
% covs(:,:,i)=0.5*(inv(h)+inv(h)');
% mus(i,:)=vardur_model{i}.params_phi_map;
% i
% end
% 
% 
% save('vardur_agg_summary_v22.mat')
load('vardur_agg_summary_v22.mat','mus','covs')






xs=[1:0.1:100];a=0.5*(xs(1:end-1)+xs(2:end));
clear posts

cc=cbrewer('div','RdGy',11);
cc=cc(8:10,:);
for kk=1:3
    for i=1:460
        post=get_post_nsamp_cov(mus(i,1+kk),covs(1+kk,1+kk,i),xs);
        posts{kk}(i,:)=post/mean(diff(xs));
        i
    end
    
    
    
    
end


%%

for kk=1:3
%     bar(a,mean(posts{kk},1),1,'facecolor',cc(kk,:),'facealpha',0.5);
    
    plot(a,mean(posts{kk},1),'color',cc(kk,:),'linewidth',2);
    hold on
end

% plot([1,100],[0.01,0.01],'m-.','linewidth',2);

xticks([1,25,50,75,100])
yticks([0.5,1:5,10]*(1/numel(a)))
ylim([0.4,15]*(1/numel(a)));
xlim([1,100]);
xticklabels({'1','25','','','100'});
% xlabel('num. of samples');

dd=max(xs)-min(xs);
yticks([1,2:2:8]*(1/dd));
ylim([0,8]*(1/dd));
yticklabels(round([1,2:2:8]*(1/dd),2));


% yticklabels({'0.5','1','','','','5','10'})

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'YScale','linear')
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')


standardize_figure(1,[2,1.5])
% xlabel('number of samples')
% ylabel('normalized likelihood')



saveas(gcf,['plots/fig_8C'],'pdf');


