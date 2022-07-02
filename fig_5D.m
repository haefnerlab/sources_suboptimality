load('bias_cov','xs','ys','mus');
for i=1:20
    %     load(['post_samps_v22/fixed/subj_samps_',num2str(i),'.mat'],'subj');
    %     ss1=subj.phi_theta(get_samples_burnin_thin(subj.post_samps,1000,0,46,size(subj.mus{1},1)),subj);
    %     bs=[normcdf(ss1(:,6)),ss1(:,1)];
    %     mu=mean(bs);
    %     c=cov(bs);
    %     [x,y]=get_cov_ellipse(mu,c,0.68,100);
    x=xs{i};
    y=ys{i};
    plot(x,y,'linewidth',1,'color',[0.75,0.75,0.75]);
    if i==4 || i==17
        plot(x,y,'linewidth',1,'color',[0,0,0]);
    end
    
    hold on
    %     plot(mu(1),mu(2),'.','markersize',10,'color',[0.5,0.5,0.5]);
    %     if i==4 || i==7
    %         plot(mu(1),mu(2),'.','markersize',10,'color',[0,0,0]);
    %     end
    %     mus(i,:)=mu;
    drawnow
    xs{i}=x;
    ys{i}=y;
end










for i=1:20
    plot(mus(i,1),mus(i,2),'.','markersize',10,'color',[0.5,0.5,0.5]);
    if i==4 || i==7
        plot(mus(i,1),mus(i,2),'.','markersize',10,'color',[0,0,0]);
    end
end


for i=[4,17]
    %     load(['post_samps_v22/fixed/subj_samps_',num2str(i),'.mat'],'subj');
    %     ss1=subj.phi_theta(get_samples_burnin_thin(subj.post_samps,1000,0,46,size(subj.mus{1},1)),subj);
    %     bs=[normcdf(ss1(:,6)),ss1(:,1)];
    %     mu=mean(bs);
    %     c=cov(bs);
    %     [x,y]=get_cov_ellipse(mu,c,0.68,100);
    x=xs{i};
    y=ys{i};
    plot(x,y,'linewidth',1,'color',[0.75,0.75,0.75]);
    if i==4 || i==17
        plot(x,y,'linewidth',1,'color',[0,0,0]);
    end
    
    hold on
    drawnow
end


for i=1:20
    
    if i==4 || i==17
        plot(mus(i,1),mus(i,2),'.','markersize',10,'color',[0,0,0]);
    end
end

% plot(mean(mus(:,1)),mean(mus(:,2)),'.','markersize',10,'color',[0,0,0]);
xlim([0,1]);
ylim([-1,1]);
% plot([0,1],[0,0],'k:')
% plot([0.5,0.5],[-1,1],'k:')
load(['post_samps_v22/fixed/subj_samps_',num2str(i),'.mat'],'subj');
model=subj;
ss0=model.phi_theta(normrnd(0,1,[100000,11]),model);
bs=[normcdf(ss0(:,6)),ss0(:,1)];
mu=mean(bs);
c=cov(bs);
[x,y]=get_cov_ellipse(mu,c,0.68,100);
plot(x,y,':','linewidth',1,'color',[0,0,0]);
% xlabel('categorical bias implied by perceptual bias')
% ylabel('calibration bias')
standardize_figure(1,[2,1.5])
saveas(gcf,['plots/fig_5D'],'pdf');






% cc=cbrewer('div','RdYlGn',11);
% cc=cc([3,9],:);
%
% figure(1);
%
%
% hold on
%
%
%
% a=[-1-(0.02/2):0.02:1+(0.02/2)];
% xs=0.5*(a(1:end-1)+a(2:end));
%
% for ii=1:20
%     load(['post_samps_v22/fixed/subj_samps_',num2str(ii),'.mat']);
%     ss1=subj.phi_theta(get_samples_burnin_thin(subj.post_samps,500,0,46,size(subj.mus{1},1)),subj);
%     q1=normpdf(xs,mean(ss1(:,1)),std(ss1(:,1)));
%
%
%
%
% %     [q1,q2,q3]=histcounts(ss1(:,1),a);
% %     q1=q1/sum(q1*mean(diff(xs)));
%     plot(xs,q1,'color',[0.75,0.75,0.75],'linewidth',1);
%     q1s(ii,:)=q1;
%     hold on
%     drawnow
% end
%
%
%
% ids=[4,17];
%
% for ii=1:numel(ids)
%
% % plot(xs,q1s(ids(ii),:),'color',cc(ii,:),'linewidth',2);
% plot(xs,q1s(ids(ii),:),'color',[0,0,0],'linewidth',1);
% hold on
% drawnow
% end
%
%
% % plot(xs,normpdf(xs,0,0.25),'m-.','linewidth',2)
%
% plot(xs,mean(q1s,1),'color',[0,0,0],'linewidth',3);
%
%
%
%
%
% plot([0,0],[0,4],'k:','linewidth',1.5);
%
%
%
%
%
% xticks([-1:0.5:1])
% yticks([0:2:4])
%
% yticklabels([])
%
%
% standardize_figure(1,[2,1.5])
% % xlabel('choice bias')
% % ylabel('posterior probability')
%
%
% saveas(gcf,['../../plots/fig4d1'],'pdf');
