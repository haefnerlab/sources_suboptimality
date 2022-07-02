load thresh_vardur

cc=cbrewer('div','PiYG',11);
cols=cc([10,2],:);


figure(4);
for i=1:18
mu_cen(i,:)=100*(1-nanmean((thresh_cen_exact{i}./thresh_cen{i}),1));
plot([1,2,3],mu_cen(i,:),'.-','color',cols(1,:)*0.25+[1,1,1]*0.75,'markersize',10);
hold on
drawnow
end
for i=1:18
mu_mat(i,:)=100*(1-nanmean((thresh_mat_exact{i}./thresh_mat{i}),1));
plot([1,2,3],mu_mat(i,:),'.-','color',cols(2,:)*0.25+[1,1,1]*0.75,'markersize',10);
hold on
drawnow
end
errorbar([1,2,3]-0.1,mean(mu_cen,1),std(mu_cen,[],1),'.-','color',cols(1,:),'markersize',20,'linewidth',1.5,'capsize',0);
errorbar([1,2,3]+0.1,mean(mu_mat,1),std(mu_mat,[],1),'.-','color',cols(2,:),'markersize',20,'linewidth',1.5,'capsize',0);

% sigstar({[1,2],[2,3]},[signtest(mu_cen(:,1),mu_cen(:,2)),signtest(mu_cen(:,2),mu_cen(:,3))])


xlim([0.5,3.5]);
xticks([1,2,3]);
xticklabels({'100','300','1000'});
ylim(100*[0,0.25]);
yticks([0:5:25]);
xtickangle(0);
% xlabel('duration (ms)');
% ylabel('rel. contribution of sens. noise');


standardize_figure(4,[1.5,1.5]);
saveas(gcf,['plots/fig_S2'],'pdf');

pvals_cen(4,:)=[signtest(mu_cen(:,1),mu_cen(:,2)),signtest(mu_cen(:,2),mu_cen(:,3))];
pvals_mat(4,:)=[signtest(mu_mat(:,1),mu_mat(:,2)),signtest(mu_mat(:,2),mu_mat(:,3))];