








xs=[1:0.1:100];a=0.5*(xs(1:end-1)+xs(2:end));

% cols=eye(3);


for ii=1:18
    for kk=1:3
%         figure(kk);
%         hold on
        load(['post_samps_v22/vardur/subj_samps_vardur_bfit_',num2str(ii),'.mat']);
        ss1=subj.phi_theta(get_samples_burnin_thin(subj.post_samps,100,0,46,size(subj.mus{1},1)),subj);
        [q1,q2,q3]=histcounts(ss1(:,1+kk),xs);
%         q1=q1/sum(q1*mean(diff(a)));
        q1=q1/sum(q1);
%         plot(a,q1,'color',[0.75,0.75,1],'linewidth',1);
        q1s{kk}(ii,:)=q1;
        hold on
        drawnow
        
    end
    ii
end


% ids=[4,7];
%

ids=[2];
cc=cbrewer('div','BrBG',11);
cc=cc([3],:);

%%

for ii=1:18
    for kk=1:3
        figure(kk);
        plot(a,smooth(q1s{kk}(ii,:)./mean(diff(xs)),10),'color',[0.75,0.75,0.75],'linewidth',1);
        
        hold on
        drawnow
        
    end

end




for ii=1:numel(ids)
for kk=1:3
    figure(kk);
% plot(a,q1s{kk}(ids(ii),:),'color',cc(ii,:),'linewidth',1.5);
plot(a,smooth(q1s{kk}(ids(ii),:)./mean(diff(xs)),10),'color',[0,0,0],'linewidth',1);
hold on
drawnow
end
end
for kk=1:3
    figure(kk);
plot(a,smooth(mean(q1s{kk},1)./mean(diff(xs)),10),'color','k','linewidth',3);
end


for kk=1:3
    figure(kk);
% plot([1,100],[0.01,0.01],'m-.','linewidth',2);


xticks([1,25,50,75,100])
yticks([0.5,1:5]*(1/numel(a)))
ylim([0.4,6]*(1/numel(a)));
xlim([1,100]);
xticklabels({'1','25','','','100'});

dd=max(xs)-min(xs);
yticks([1,2:2:8]*(1/dd));
ylim([0,8]*(1/dd));
yticklabels(round([1,2:2:8]*(1/dd),2));


% yticklabels({'0.5','1','','','','5'})

set(gca,'XScale','log')
% set(gca,'YScale','log')
set(gca,'YScale','linear')
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')


standardize_figure(kk,[2,1.5])
% xlabel('num. of samples')
% ylabel({'likelihood function','(normalized)'})









% % xticks([1,25,50,75,100])
% % yticks([0:0.01:0.05])
% % 
% % ylim([0,0.06]);
% % 
% % if kk==3
% % xticklabels({'1','25','','','100'});
% % xlabel('num. of samples');
% % 
% % else
% % xticklabels('')
% % end
% % yticklabels({'','0.01','','','','0.05'})
% % 
% % set(gca,'XScale','log')
% % set(gca,'YScale','log')
% % set(gca,'XMinorTick','off')
% % set(gca,'YMinorTick','off')
% % 
% % 
% % standardize_figure(kk,[2.5,2])
% % % if kk==2
% % ylabel('likelihood fn.')
% % % end
saveas(gcf,['plots/fig_7D_',num2str(kk)],'pdf');
end




