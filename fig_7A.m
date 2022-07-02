ii=2;
load(['post_samps_v22/vardur/subj_samps_vardur_bfit_',num2str(ii),'.mat']);
subj.params_phi_map=subj.params_phi_map_theta;
for kk=1:6
subj.design_matrix_pred{kk}.npts=1001;
subj.design_matrix_pred{kk}.eps_v_right=linspace(min(subj.design_matrix_pred{kk}.eps_v_right),max(subj.design_matrix_pred{kk}.eps_v_right),1001);
subj.design_matrix_pred{kk}.eps_a_tone=linspace(min(subj.design_matrix_pred{kk}.eps_a_tone),max(subj.design_matrix_pred{kk}.eps_a_tone),1001);
end
model=subj;


cols=[repmat([0,0,1],3,1);repmat([1,0,0],3,1)];
cc=cbrewer('div','PiYG',11);
cols=cc([10,2],:);
cols=[repmat(cols(1,:),3,1);repmat(cols(2,:),3,1)];


model.num_dur=3;
num_pts=10;
tmp=[];
for i=1:numel(model.design_matrix{1}.eps_a_tone)
    tmp=[tmp;model.design_matrix{1}.eps_a_tone(:)];
end

if isfield(model,'params_phi_mle')
    
    pred=model.simulate_prob_resp_pred(model.phi_theta(model.params_phi_mle,model),model);

elseif isfield(model,'params_phi_map')
    pred=model.simulate_prob_resp_pred(model.phi_theta(model.params_phi_map,model),model);

else
    pred=model.simulate_prob_resp_pred(model.phi_theta(model.params_phi,model),model);
end


for i=1:2
    for j=1:3
        figure(j);

        es=[model.design_matrix{j}.eps_a_tone(:);model.design_matrix{3+j}.eps_a_tone(:)];
        [f,x]=ecdf(es);
        bin_width=0.1;
        clear id
        for k=1:((1/bin_width)-1)
            id(k)=find(f>=bin_width*k,1);
        end
        id=[1,id,numel(f)];
        ed=x(id);
        es=model.design_matrix{(i-1)*3+j}.eps_a_tone(:);        
        [q1,q2,q3]=histcounts(es,ed);
        ns=model.data{(i-1)*3+j}.num_repeats(:);
        ks=model.data{(i-1)*3+j}.num_ch1(:);
        
        
        q1=q1(:);

        clear xss nss kss
        uq3=unique(q3);
        for k=1:numel(uq3)
            xss(k)=0.5*(q2(uq3(k))+q2(uq3(k)+1));
            nss(k)=sum(ns(q3==uq3(k)));
            kss(k)=sum(ks(q3==uq3(k)));
        end
        
        
        
        mu=kss./nss;
        se=sqrt(mu.*(1-mu)./nss);
        
        errorbar(xss,mu,se,'.','Color',cols((i-1)*3+j,:),'linewidth',1,'markersize',20)
        hold on
        plot(model.design_matrix_pred{(i-1)*3+j}.eps_a_tone,pred{(i-1)*3+j},'Color',cols((i-1)*3+j,:),'linewidth',2);
    end
end











load(['post_samps_v22/vardur_null1/subj_samps_vardur_bfit_',num2str(ii),'.mat']);
subj.params_phi_map=subj.params_phi_map_theta;
model=subj;


cols=[repmat([0,0,1],3,1);repmat([1,0,0],3,1)];
cc=cbrewer('div','PiYG',11);
cols=cc([10,2],:);
cols=[repmat(cols(1,:),3,1);repmat(cols(2,:),3,1)];


model.num_dur=3;
num_pts=10;
tmp=[];
for i=1:numel(model.design_matrix{1}.eps_a_tone)
    tmp=[tmp;model.design_matrix{1}.eps_a_tone(:)];
end

if isfield(model,'params_phi_mle')
    
    pred=model.simulate_prob_resp_pred(model.phi_theta(model.params_phi_mle,model),model);

elseif isfield(model,'params_phi_map')
    pred=model.simulate_prob_resp_pred(model.phi_theta(model.params_phi_map,model),model);

else
    pred=model.simulate_prob_resp_pred(model.phi_theta(model.params_phi,model),model);
end


for i=1:2
    for j=1:3
        figure(j);

        es=[model.design_matrix{j}.eps_a_tone(:);model.design_matrix{3+j}.eps_a_tone(:)];
        [f,x]=ecdf(es);
        bin_width=0.1;
        clear id
        for k=1:((1/bin_width)-1)
            id(k)=find(f>=bin_width*k,1);
        end
        id=[1,id,numel(f)];
        ed=x(id);
        es=model.design_matrix{(i-1)*3+j}.eps_a_tone(:);        
        [q1,q2,q3]=histcounts(es,ed);
        ns=model.data{(i-1)*3+j}.num_repeats(:);
        ks=model.data{(i-1)*3+j}.num_ch1(:);
        
        q1=q1(:);

        clear xss nss kss
        uq3=unique(q3);
        for k=1:numel(uq3)
            xss(k)=0.5*(q2(uq3(k))+q2(uq3(k)+1));
            nss(k)=sum(ns(q3==uq3(k)));
            kss(k)=sum(ks(q3==uq3(k)));
        end
        mu=kss./nss;
        se=sqrt(mu.*(1-mu)./nss);
        
%         errorbar(xss,mu,se,'.','Color',cols((i-1)*3+j,:),'linewidth',1,'markersize',20)
        hold on
        plot(model.design_matrix_pred{(i-1)*3+j}.eps_a_tone,pred{(i-1)*3+j},'Color',cols((i-1)*3+j,:),'linewidth',2,'linestyle',':');
    end
end












for j=1:3
    figure(j);

%     xlabel('tone position (deg)');
%     ylabel('prop. ''right'' response');
    
    
    
    
    
    
    
    set(gca,'fontsize',12,'fontweight','normal','fontname','Helvetica Neue')
    set(gcf,'color','white')
    set(gca,'linewidth',2)
    set(gca,'TickDir','out');
    
    box off
    xlim([-25,25]);
    ylim([0,1]);
    xticks([-25:12.5:25]);
    yticks([0:0.5:1]);
%     xticklabels('');
%     if j>1
%         yticks([]);
%         ylabel('');
%         xlabel('');
%     end
    
    
    axis fill
    plot([-25,25],[0.5,0.5],'k:','linewidth',1.5);
    plot([0,0],[0,1],'k:','linewidth',1.5);
    
    
    standardize_figure(j,[2,1.5])

saveas(gcf,['plots/fig_7A_',num2str(j)],'pdf');
    
end