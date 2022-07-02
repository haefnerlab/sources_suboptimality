cc=cbrewer('div','RdYlGn',11);
cc=cc([3,9],:);


figure(1);


hold on



xs=[1:0.1:100];a=0.5*(xs(1:end-1)+xs(2:end));

for ii=1:20
    load(['post_samps_v22/fixed/subj_samps_',num2str(ii),'.mat']);
    ss1=subj.phi_theta(get_samples_burnin_thin(subj.post_samps,500,0,46,size(subj.mus{1},1)),subj);
    [q1,q2,q3]=histcounts(ss1(:,2),xs);
    q1=q1/(sum(q1)*mean(diff(xs)));
    plot(a,smooth(q1),'color',[0.75,0.75,0.75],'linewidth',1);
    q1s(ii,:)=q1;
    hold on
    drawnow
end
%%

for ii=1:20
    
    plot(a,(q1s(ii,:)),'color',[0.75,0.75,0.75],'linewidth',1);
    
    hold on
    drawnow
end


ids=[4,17];

for ii=1:numel(ids)

% plot(a,q1s(ids(ii),:),'color',cc(ii,:),'linewidth',1.5);
plot(a,(q1s(ids(ii),:)),'color',[0,0,0],'linewidth',1);
hold on
drawnow
end

plot(a,(mean(q1s,1)),'color',[0,0,0],'linewidth',3);


% plot([1,100],[0.01,0.01],'m-.','linewidth',2);

dd=max(xs)-min(xs);

xticks([1,25,50,75,100])
yticks([1,2:2:8]*(1/dd));
ylim([0,8]*(1/dd));
xlim([1,100]);
xticklabels({'1','25','','','100'});





% yticklabels({'0.2','0.5','1','','','','5'})
yticklabels(round([1,2:2:8]*(1/dd),2));

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'XMinorTick','off')
set(gca,'YMinorTick','off')

set(gca,'YScale','linear')


standardize_figure(1,[2,1.5])

saveas(gcf,['plots/fig_5C'],'pdf');

