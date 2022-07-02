function model = generate_model_struct_v22_vardur_null_exact(id,sub_init,incomplete_dir_save)
model.params_names = {...   
'Choice prior',...
'Number of samples (s)',...
'Number of samples (m)',...
'Number of samples (l)',...
'Prior combination prob.',...
'Aud. sensor noise (std;s)',...
'Aud. sensor noise (std;m)',...
'Aud. sensor noise (std;l)',...
'Vis. sensor noise (std;s)',...
'Vis. sensor noise (std;m)',...
'Vis. sensor noise (std;l)',...
'Aud. prior mean',...
'Aud. prior std',...
'Vis. prior mean',...
'Vis. prior std',...
'Lapse rate',...
'Lapse prob',...
'eps0',...
'd'...
};


model.id = id;

if exist('incomplete_dir_save','var')
    model.incomplete_dir_save=incomplete_dir_save;
end

% Parameter family:
% Beta(a,b)=1
% Gamma(k,theta)=2
% Normal(mu,sig)=3
% Lognormal(mu,sig)=4
% Beta-prime(mu,sig)=5




model.min_sig=1e-4;
model.min_lapse=1e-4;

model.params_family = [3 , 1, 5,5,5, 5,5,5, 3,4,3,3,3,3    1,1,1,1,1,1,        1,1];
model.params_additive_offset = [0,         0,     0,0,0,    0,0,0, 0,model.min_sig,0,0,0,0,       model.min_lapse,model.min_lapse,model.min_lapse,model.min_lapse,model.min_lapse,model.min_lapse,    0,0.5];
model.params_multiplicative_offset = [1, 1, 1,1,1,   1,1,1, 1,1,1,1,1,1,        1-model.min_lapse,1-model.min_lapse,1-model.min_lapse,1-model.min_lapse,1-model.min_lapse,1-model.min_lapse,   1,1];
model.params_hyperprior_params{1} = [0,       1.25,  1,1,1,   1,1,1,   0,0,0,0,0,0,      1,1.25,1,1.25,1,1.25,     1,2];
model.params_hyperprior_params{2} = [0.25,       1.25,  1,1,1,    1,1,1,   0.5,2.35,0.5*2,2.35*2,0.5*2,2.35*2,  5,1.25,5,1.25,5,1.25,    1,2];%3.525,2];



model.exact_inference=1;
model.normprior_std = 1;
model.max_phi=6*model.normprior_std;
model.ns_inf=100;
model.eps_max=atand((0.5*60.96)/50); %Maximum eccentricity

model.num_params = length(model.params_family);
model.num_eff_params = model.num_params;

model.phi_theta = @(phi, model) phi_theta(phi, model);
model.theta_phi = @(theta, model) theta_phi(theta, model);

model.npts_quad = 20;
[model.pts_leg, model.wts_leg] = GaussLegendre(model.npts_quad);

model.num_dur=numel(sub_init.eps_a_raw_central);
for kk=1:numel(sub_init.eps_a_raw_central)
    model.design_matrix{kk}.eps_a_tone=0.5*sub_init.eps_a_raw_central{kk};
    model.design_matrix{kk}.eps_v_right=0*model.design_matrix{kk}.eps_a_tone;
    model.data{kk}.num_ch1=sub_init.resp_raw_central{kk};
    model.data{kk}.num_repeats=ones(size(sub_init.resp_raw_central{kk}));
    model.design_matrix{kk}.npts=length(model.design_matrix{kk}.eps_a_tone);
    model.design_matrix_pred{kk}.eps_a_tone=0.5*linspace(min(sub_init.eps_a_raw_central{kk}),max(sub_init.eps_a_raw_central{kk}),101);
    model.design_matrix_pred{kk}.eps_v_right=0*model.design_matrix_pred{kk}.eps_a_tone;
    model.design_matrix_pred{kk}.npts=length(model.design_matrix_pred{kk}.eps_a_tone);
end
for kk=1:numel(sub_init.eps_a_raw_matched)
    model.design_matrix{numel(sub_init.eps_a_raw_central)+kk}.eps_a_tone=0.5*sub_init.eps_a_raw_matched{kk};
    model.design_matrix{numel(sub_init.eps_a_raw_central)+kk}.eps_v_right=abs(model.design_matrix{numel(sub_init.eps_a_raw_central)+kk}.eps_a_tone);
    model.data{numel(sub_init.eps_a_raw_central)+kk}.num_ch1=sub_init.resp_raw_matched{kk};
    model.data{numel(sub_init.eps_a_raw_central)+kk}.num_repeats=ones(size(sub_init.resp_raw_matched{kk}));
    model.design_matrix{numel(sub_init.eps_a_raw_central)+kk}.npts=length(model.design_matrix{numel(sub_init.eps_a_raw_central)+kk}.eps_a_tone);
    model.design_matrix_pred{numel(sub_init.eps_a_raw_central)+kk}.eps_a_tone=0.5*linspace(min(sub_init.eps_a_raw_matched{kk}),max(sub_init.eps_a_raw_matched{kk}),101);
    model.design_matrix_pred{numel(sub_init.eps_a_raw_central)+kk}.eps_v_right=abs(model.design_matrix_pred{numel(sub_init.eps_a_raw_central)+kk}.eps_a_tone);
    model.design_matrix_pred{numel(sub_init.eps_a_raw_central)+kk}.npts=length(model.design_matrix_pred{numel(sub_init.eps_a_raw_central)+kk}.eps_a_tone);
end





model.set_default = [11:14];
model.default_values = [zeros(1,4)];
[~, ids_sort] = sort(model.set_default, 'descend');
model.num_eff_params = model.num_params - numel(model.set_default);
model.set_default = model.set_default(ids_sort);
model.default_values = model.default_values(ids_sort);

model.params_family_eff = model.params_family;
model.params_additive_offset_eff = model.params_additive_offset;
model.params_multiplicative_offset_eff = model.params_multiplicative_offset;
model.params_hyperprior_params_eff = model.params_hyperprior_params;
model.params_family_eff(model.set_default) = [];
model.params_additive_offset_eff(model.set_default) = [];
model.params_multiplicative_offset_eff(model.set_default) = [];
for jj = 1:numel(model.params_hyperprior_params)
    model.params_hyperprior_params_eff{jj}(model.set_default) = [];
end
model.params_theta = rand(1, model.num_eff_params);
model.params_theta(model.params_family_eff == 1) = betainv(model.params_theta(model.params_family_eff == 1), ...
    model.params_hyperprior_params_eff{1}(model.params_family_eff == 1), ...
    model.params_hyperprior_params_eff{2}(model.params_family_eff == 1));
model.params_theta(model.params_family_eff == 2) = gaminv(model.params_theta(model.params_family_eff == 2), ...
    model.params_hyperprior_params_eff{1}(model.params_family_eff == 2), ...
    model.params_hyperprior_params_eff{2}(model.params_family_eff == 2));
model.params_theta(model.params_family_eff == 3) = norminv(model.params_theta(model.params_family_eff == 3), ...
    model.params_hyperprior_params_eff{1}(model.params_family_eff == 3), ...
    sqrt(model.params_hyperprior_params_eff{2}(model.params_family_eff == 3)));
model.params_theta = model.params_multiplicative_offset_eff .* model.params_theta + model.params_additive_offset_eff;
model.params_phi = model.theta_phi(model.params_theta, model);

model.simulate_data = @(model, num_repeats) simulate_data(model, num_repeats);
model.simulate_prob_resp = @(params, model) simulate_prob_resp(params, model);
model.simulate_prob_resp_pred = @(params, model) simulate_prob_resp_pred(params, model);
model.simulate_prob_resp_pred_specific= @(params,eps_a,eps_v,dur,model) simulate_prob_resp_pred_specific(params,eps_a,eps_v,dur, model);
% model.plot_model_pred=@(model,figid,showlegend) plot_model_pred(model,figid,showlegend);
model.plot_model_pred = @(model, plt_sz,plt_ids, showlegend) plot_model_pred(model, plt_sz,plt_ids, showlegend);

model.get_samples_posterior = @(model, numchains, num_samples, num_samples_save) get_samples_posterior(model, numchains, num_samples, num_samples_save);
model.get_samples_posterior_parallel = @(model, numchains, num_samples, num_samples_save) get_samples_posterior_parallel(model, numchains, num_samples, num_samples_save);

model.log_unnorm_post = @(params, model) log_unnorm_post(params, model);
model.log_likelihood = @(params, model) log_likelihood(params, model);
model.log_prior_theta = @(theta, model) log_prior_theta(theta, model);
model.log_prior_phi = @(phi, model) log_prior_phi(phi, model);

model.get_map = @(model, nstart, gethessian, useparallel) get_map(model, nstart, gethessian, useparallel);
model.get_map_theta = @(model, nstart, gethessian, useparallel) get_map_theta(model, nstart, gethessian, useparallel);
model.get_mle = @(model, nstart, gethessian, useparallel) get_mle(model, nstart, gethessian, useparallel);

model.get_bootstraps_map = @(model, nboot) get_bootstraps_map(model, nboot);
model.get_bootstraps_mle = @(model, nboot) get_bootstraps_mle(model, nboot);
end

function theta = phi_theta(phi, model)
phi=min(model.max_phi,max(-1*model.max_phi,phi));
theta = normcdf(phi, 0, model.normprior_std);
if sum(model.params_family_eff == 1)>0
    theta(:,model.params_family_eff == 1) = betainv(theta(:,model.params_family_eff == 1), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 1),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 1),size(phi,1),1));
end
if sum(model.params_family_eff == 2)>0
    theta(:,model.params_family_eff == 2) = gaminv(theta(:,model.params_family_eff == 2), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 2),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 2),size(phi,1),1));
end
if sum(model.params_family_eff == 3)>0
    theta(:,model.params_family_eff == 3) = norminv(theta(:,model.params_family_eff == 3), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 3),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 3),size(phi,1),1));
end
if sum(model.params_family_eff == 4)>0
    theta(:,model.params_family_eff == 4) = logninv(theta(:,model.params_family_eff == 4), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 4),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 4),size(phi,1),1));
end
if sum(model.params_family_eff == 5)>0
    theta(:,model.params_family_eff == 5) = betaprinv(theta(:,model.params_family_eff == 5), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 5),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 5),size(phi,1),1));
end



theta = theta.*repmat(model.params_multiplicative_offset_eff,size(theta,1),1) + repmat(model.params_additive_offset_eff,size(theta,1),1);
end

function phi = theta_phi(theta, model)

phi = (theta - repmat(model.params_additive_offset_eff,size(theta,1),1)) ./ repmat(model.params_multiplicative_offset_eff,size(theta,1),1);
if sum(model.params_family_eff == 1)>0
    phi(:,model.params_family_eff == 1) = betacdf(phi(:,model.params_family_eff == 1), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 1),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 1),size(phi,1),1));
end
if sum(model.params_family_eff == 2)>0
    phi(:,model.params_family_eff == 2) = gamcdf(phi(:,model.params_family_eff == 2), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 2),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 2),size(phi,1),1));
end
if sum(model.params_family_eff == 3)>0
    phi(:,model.params_family_eff == 3) = normcdf(phi(:,model.params_family_eff == 3), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 3),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 3),size(phi,1),1));
end
if sum(model.params_family_eff == 4)>0
    phi(:,model.params_family_eff == 4) = logncdf(phi(:,model.params_family_eff == 4), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 4),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 4),size(phi,1),1));
end
if sum(model.params_family_eff == 5)>0
    phi(:,model.params_family_eff == 5) = betaprcdf(phi(:,model.params_family_eff == 5), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 5),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 5),size(phi,1),1));
end

phi = norminv(phi, 0, model.normprior_std);
phi=min(model.max_phi,max(-1*model.max_phi,phi));
end

function plot_model_pred(model, plt_sz,plt_ids, showlegend)
cols=[repmat([0,0,0],3,1);repmat([1,0,0],3,1)];
% cc=cbrewer('div','RdGy',11);
% cols=cc([2,3,4,10,9,8],:);

model.num_dur=3;
num_pts=10;
tmp=[];
for i=1:numel(model.design_matrix{1}.eps_a_tone)
    tmp=[tmp;model.design_matrix{1}.eps_a_tone(:)];
end

if isfield(model,'params_phi_mle')
    %     pt=get_full_params(model.phi_theta(model.params_phi_mle,model),model);
    pt=model.phi_theta(model.params_phi_mle,model);
    pred=simulate_prob_resp_pred(pt,model);
    
elseif isfield(model,'params_phi_map')
    %     pt=get_full_params(model.phi_theta(model.params_phi_map,model),model);
    pt=model.phi_theta(model.params_phi_map,model);
    pred=simulate_prob_resp_pred(pt,model);
    
else
    %     pt=get_full_params(model.phi_theta(model.params_phi,model),model);
    pt=model.phi_theta(model.params_phi,model);
    pred=simulate_prob_resp_pred(pt,model);
end


for i=1:2
    for j=1:3
        subplot(plt_sz(1),plt_sz(2),plt_ids(j));
% figure(j)
        % %         [q1,q2,q3]=histcounts(model.design_matrix{(i-1)*3+j}.eps_a_tone,num_pts);
        % %         for k=1:num_pts
        % %             xs(k)=0.5*(q2(k)+q2(k+1));
        % %             ys(k)=sum(model.data{(i-1)*3+j}.num_ch1(q3==k))./sum(q3==k);
        % %             se(k)=sqrt((ys(k)*(1-ys(k)))/sum(q3==k));
        % %         end
        
        es=model.design_matrix{(i-1)*3+j}.eps_a_tone(:);
        ns=model.data{(i-1)*3+j}.num_repeats(:);
        ks=model.data{(i-1)*3+j}.num_ch1(:);
        [f,x]=ecdf(es);
        bin_width=0.1;
        clear id
        for k=1:((1/bin_width)-1)
            id(k)=find(f>=bin_width*k,1);
        end
        id=[1,id,numel(f)];
        ed=x(id);
        [q1,q2,q3]=histcounts(es,ed);
        q1=q1(:);
%         id_del=diff(q2)>2;%  | q1<1;
%         for k1=1:numel(q1)
%             if id_del(k1)==1
%                 es(q3==k1)=[];
%                 ns(q3==k1)=[];    ks(q3==k1)=[];
%                 q3(q3==k1)=[];
%             end
%         end
        clear xss nss kss
        uq3=unique(q3);
        for k=1:numel(uq3)
            xss(k)=0.5*(q2(uq3(k))+q2(uq3(k)+1));
            nss(k)=sum(ns(q3==uq3(k)));
            kss(k)=sum(ks(q3==uq3(k)));
        end
        mu=kss./nss;
        se=sqrt(mu.*(1-mu)./nss);
        
        errorbar(xss,mu,se,'.','Color',cols((i-1)*3+j,:),'linewidth',2,'markersize',30)
        hold on
        plot(model.design_matrix_pred{(i-1)*3+j}.eps_a_tone,pred{(i-1)*3+j},'Color',cols((i-1)*3+j,:),'linewidth',2);
    end
end

for j=1:3
    subplot(plt_sz(1),plt_sz(2),plt_ids(j));
% figure(j)
    if showlegend==1
        legend(pl1,{'Central pred.','Matched pred.'},'Location', 'best');
    end
    xlabel('Tone position');
    ylabel('Probability of responding right');
    set(gca,'fontsize',12,'fontweight','bold');
    set(gcf,'color','white')
    set(gca,'linewidth',4)
    set(gca,'TickDir','out');
    if showlegend==1
        legend boxoff
    end
    box off
    xlim([-1*max(abs(model.design_matrix_pred{j}.eps_a_tone)),max(abs(model.design_matrix_pred{j}.eps_a_tone))]);
    ylim([0,1]);
    hline(0.5,'b--');
    vline(0,'b--');
%     axis square
end

end

function model = get_samples_posterior(model, numchains, num_samples, num_samples_save)
model = gess_model(model, numchains, num_samples, num_samples_save, 0, 0);
end

function model = get_samples_posterior_parallel(model, numchains, num_samples, num_samples_save)



if isfield(model, 'incomplete_dir_save')
    
    model = gess_parallel_v2(model.incomplete_dir_save,model, numchains, num_samples, num_samples_save, 0, 0);
else
    model = gess_parallel(model, numchains, num_samples, num_samples_save, 0, 0);
end
end

function model1 = get_bootstraps_map(model, nboot)


options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e4,'MaxIterations',1e4);

for i = 1:nboot
    model1{i} = model;
    model1{i}.data{1}.num_ch1 = binornd(model1{i}.data{1}.num_repeats, model1{i}.data{1}.num_ch1 ./ model1{i}.data{1}.num_repeats);
end

parfor i = 1:nboot
    i
    
    param0 = normrnd(0, model.normprior_std, [1, model1{i}.num_eff_params]);
    boots(i, :) = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x, model), model)), param0, options);
end
model.boots = boots;
end

function model1 = get_bootstraps_mle(model, nboot)


options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e4,'MaxIterations',1e4);

for i = 1:nboot
    model1{i} = model;
    model1{i}.data{1}.num_ch1 = binornd(model1{i}.data{1}.num_repeats, model1{i}.data{1}.num_ch1 ./ model1{i}.data{1}.num_repeats);
end

parfor i = 1:nboot
    i
    
    param0 = normrnd(0, model.normprior_std, [1, model1{i}.num_eff_params]);
    boots(i, :) = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
end
model.boots = boots;
end







function model = get_map_theta(model, nstart, gethessian, useparallel)

options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e5,'MaxIterations',1e5);
if useparallel == 1
    parfor i = 1:nstart
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        
        if gethessian == 1
            
%             [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x,model), model)), param0, options);
            
            
            
            
        else
%             [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x,model), model)), param0, options);
            %             [fit(i,:),neg_lups(i)]=fminunc(@(x) -1*model.log_unnorm_post(x,model),param0,options);
            %         [fit(i,:),neg_lups(i)]=bads(@(x) -1*model.log_unnorm_post(x,model),param0,-100*ones(1,model.num_params),100*ones(1,model.num_params));
        end
    end
else
    
    for i = 1:nstart
%         param0 = model.params_phi_map;
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        if gethessian == 1
%             [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x,model), model)), param0, options);
        else
%             [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x,model), model)), param0, options);
        end
    end
end

model.params_phi_map = fit(find(neg_lups == min(neg_lups), 1), :);
model.map_allparams=fit;
model.map_neg_lups=neg_lups;
if gethessian == 1
    model.params_phi_hessian = hes{find(neg_lups == min(neg_lups), 1)};
    
end
end



function model = get_map(model, nstart, gethessian, useparallel)

options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e5,'MaxIterations',1e5);
if useparallel == 1
    parfor i = 1:nstart
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        
        if gethessian == 1
            
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
%             [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x,model), model)), param0, options);
            
            
            
            
        else
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
%             [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x,model), model)), param0, options);
            %             [fit(i,:),neg_lups(i)]=fminunc(@(x) -1*model.log_unnorm_post(x,model),param0,options);
            %         [fit(i,:),neg_lups(i)]=bads(@(x) -1*model.log_unnorm_post(x,model),param0,-100*ones(1,model.num_params),100*ones(1,model.num_params));
        end
    end
else
    
    for i = 1:nstart
%         param0 = model.params_phi_map;
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        if gethessian == 1
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
%             [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x,model), model)), param0, options);
        else
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
%             [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x,model), model)), param0, options);
        end
    end
end

model.params_phi_map = fit(find(neg_lups == min(neg_lups), 1), :);
model.map_allparams=fit;
model.map_neg_lups=neg_lups;
if gethessian == 1
    model.params_phi_hessian = hes{find(neg_lups == min(neg_lups), 1)};
    
end
end

function model = get_mle(model, nstart, gethessian, useparallel)

options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e4,'MaxIterations',1e4);
if useparallel == 1
    parfor i = 1:nstart
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        if gethessian == 1
            
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
        else
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
        end
    end
else
    
    for i = 1:nstart
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        if gethessian == 1
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
        else
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
        end
    end
end

model.params_phi_mle = fit(find(neg_lups == min(neg_lups), 1), :);
model.mle_allparams=fit;
model.mle_neg_ll=neg_lups;
if gethessian == 1
    model.params_phi_hessian = hes{find(neg_lups == min(neg_lups), 1)};
    
end
end

function lp_theta = log_prior_theta(theta, model)
lp_theta = 0;
theta = (theta - repmat(model.params_additive_offset_eff,size(theta,1),1)) ./ repmat(model.params_multiplicative_offset_eff,size(theta,1),1);

if sum(model.params_family_eff == 1)>0
    lp_theta = lp_theta + sum(log(betapdf(theta(:,model.params_family_eff == 1), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 1),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 1),size(theta,1),1))),2);
    
end
if sum(model.params_family_eff == 2)>0
    lp_theta = lp_theta + sum(log(gampdf(theta(:,model.params_family_eff == 2), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 2),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 2),size(theta,1),1))),2);
    
end
if sum(model.params_family_eff == 3)>0
    lp_theta = lp_theta + sum(lognormpdf(theta(:,model.params_family_eff == 3), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 3),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 3),size(theta,1),1)),2);
    
end

if sum(model.params_family_eff == 4)>0
    lp_theta = lp_theta + sum(log(lognpdf(theta(:,model.params_family_eff == 4), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 4),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 4),size(theta,1),1))),2);
    
    
end

if sum(model.params_family_eff == 5)>0
    lp_theta = lp_theta + sum(logbetaprpdf(theta(:,model.params_family_eff == 5), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 5),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 5),size(theta,1),1)),2);
    
end


if any(isinf(lp_theta) | isnan(lp_theta))
    warning('inf or nan log prior');
    lp_theta=-1e4;
    %keyboard
end
lp_theta= lp_theta - sum(log(repmat(model.params_multiplicative_offset_eff,size(theta,1),1)),2);

end

function lp_phi = log_prior_phi(phi, model)


lp_phi=sum(lognormpdf(phi,0,model.normprior_std),2);


end

function ll = log_likelihood(params_phi, model)

p_ch=simulate_prob_resp(model.phi_theta(params_phi,model),model);
ll=0;
for jj=[1:6]%1:numel(model.design_matrix)
    if any(p_ch{jj}<=0) || any(p_ch{jj}>=1)
        %keyboard
        warning('0 or 1 probability of choice');
        p_ch{jj}=min(1-eps,max(eps,p_ch{jj}));
        
    end
    
    
    ll=ll+sum((model.data{jj}.num_ch1.*log(p_ch{jj}))+((model.data{jj}.num_repeats-model.data{jj}.num_ch1).*log(1-p_ch{jj})),2);
    
end

end
function lup = log_unnorm_post(params_phi, model)
lup = model.log_likelihood(params_phi, model) + model.log_prior_phi(params_phi,model);
end

function model = simulate_data(model, num_repeats)
p_ch1=simulate_prob_resp(model.phi_theta(model.params_phi,model),model);
for jj=1:numel(model.design_matrix)
    model.data{jj}.num_repeats=num_repeats*ones(1,model.design_matrix{jj}.npts);
    model.data{jj}.num_ch1=binornd(num_repeats,p_ch1{jj});
end
end

function theta = get_full_params(theta0, model)

theta = zeros(1, model.num_params);

theta(model.set_default) = model.default_values;
theta(setdiff([1:model.num_params], model.set_default)) = theta0;

end

function p_ch1 = simulate_prob_resp(params, model)

params=get_full_params(params,model);



alp0=params(21);
d=params(22);


for jj=1:numel(model.design_matrix)
    
    if jj>model.num_dur
        
        
        lapse_rate=params(14+2*jj-1-2*model.num_dur);
        lapse_prob=params(14+2*jj-2*model.num_dur);
        kpa_2=params(2+jj-model.num_dur);
        kpv_2=params(5+jj-model.num_dur);
        sig_ap_2=params(10).^2;
        sig_a_2=sig_ap_2./kpa_2;
        sig_v_2=sig_a_2.*kpv_2;
        

        
    else
        lapse_rate=params(14+2*jj-1);
        lapse_prob=params(14+2*jj);
        kpa_2=params(2+jj);
        kpv_2=params(5+jj);
        sig_ap_2=params(10).^2;
        sig_a_2=sig_ap_2./kpa_2;
        sig_v_2=sig_a_2.*kpv_2;
        
    end
    
    
    
    
    num_inputs=model.design_matrix{jj}.npts;
    
    
    model.design_matrix{jj}.eps_a_tone=map_stim_normstim(model.design_matrix{jj}.eps_a_tone,alp0,d,model.eps_max);
    model.design_matrix{jj}.eps_v_right=map_stim_normstim(model.design_matrix{jj}.eps_v_right,alp0,d,model.eps_max);
    
    
    if model.exact_inference==0
        X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
        X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
        if jj>model.num_dur
            tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj-model.num_dur);
        else
            tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj);
        end
        inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
        wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
        p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
    else
        
        try
            nsig=3;
            ax(1).val=linspace(min(model.design_matrix{jj}.eps_a_tone)-nsig*sqrt(sig_a_2),max(model.design_matrix{jj}.eps_a_tone)+nsig*sqrt(sig_a_2),10);
            ax(2).val=linspace(min(model.design_matrix{jj}.eps_v_right)-nsig*sqrt(sig_v_2),max(model.design_matrix{jj}.eps_v_right)+nsig*sqrt(sig_v_2),10);
            Niteration=3;
            
            if jj>model.num_dur
                evalc('mdbm_sol=mdbm(ax,@(ax1) p_ch_x_exact(params,ax1(1,:),ax1(2,:),jj-model.num_dur),Niteration)');
            else
                evalc('mdbm_sol=mdbm(ax,@(ax1) p_ch_x_exact(params,ax1(1,:),ax1(2,:),jj),Niteration)');
            end
            
            
            xa=mdbm_sol.posinterp(1,:)';
            xv=mdbm_sol.posinterp(2,:)';
            [xv,idd]=sort(xv);
            xa=xa(idd);
            pa1=repmat(normcdf(max(ax(1).val),model.design_matrix{jj}.eps_a_tone(:)',sqrt(sig_a_2)),size(xa,1),1)-normcdf(repmat(xa,1,num_inputs),repmat(model.design_matrix{jj}.eps_a_tone(:)',size(xa,1),1),sqrt(sig_a_2));
            pa2=normcdf(repmat(xa,1,num_inputs),repmat(model.design_matrix{jj}.eps_a_tone(:)',size(xa,1),1),sqrt(sig_a_2))-repmat(normcdf(min(ax(1).val),model.design_matrix{jj}.eps_a_tone(:)',sqrt(sig_a_2)),size(xa,1),1);
            pv=normcdf(repmat(xv,1,num_inputs),repmat(model.design_matrix{jj}.eps_v_right(:)',size(xv,1),1),sqrt(sig_v_2));
            
            a1=sum(0.5*(pa1(1:end-1,:)+pa1(2:end,:)).*abs((pv(1:end-1,:)-pv(2:end,:))),1)./abs(pv(1,:)-pv(end,:));
            a2=sum(0.5*(pa2(1:end-1,:)+pa2(2:end,:)).*abs((pv(1:end-1,:)-pv(2:end,:))),1)./abs(pv(1,:)-pv(end,:));
            p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,a1./(a1+a2));
            p_ch1{jj}(isnan(a1) | isnan(a2))=NaN;
        catch er
            X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
            X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
            if jj>model.num_dur
                tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj-model.num_dur);
            else
                tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj);
            end
            inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
            wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
            p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
        end
    end
    
    
end
end
function p_ch1 = simulate_prob_resp_pred(params, model)

params=get_full_params(params,model);


model.num_dur=3;




alp0=params(21);
d=params(22);









for jj=1:numel(model.design_matrix_pred)
    if jj>model.num_dur
        
        
        lapse_rate=params(14+2*jj-1-2*model.num_dur);
        lapse_prob=params(14+2*jj-2*model.num_dur);
        kpa_2=params(2+jj-model.num_dur);
        kpv_2=params(5+jj-model.num_dur);
        sig_ap_2=params(10).^2;
        sig_a_2=sig_ap_2./kpa_2;
        sig_v_2=sig_a_2.*kpv_2;
        

        
    else
        lapse_rate=params(14+2*jj-1);
        lapse_prob=params(14+2*jj);
        kpa_2=params(2+jj);
        kpv_2=params(5+jj);
        sig_ap_2=params(10).^2;
        sig_a_2=sig_ap_2./kpa_2;
        sig_v_2=sig_a_2.*kpv_2;
        
    end
    num_inputs=model.design_matrix_pred{jj}.npts;
    
    model.design_matrix_pred{jj}.eps_a_tone=map_stim_normstim(model.design_matrix_pred{jj}.eps_a_tone,alp0,d,model.eps_max);
    model.design_matrix_pred{jj}.eps_v_right=map_stim_normstim(model.design_matrix_pred{jj}.eps_v_right,alp0,d,model.eps_max);
    
    
    
    if model.exact_inference==0
        X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
        X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
        if jj>model.num_dur
            tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj-model.num_dur);
        else
            tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj);
        end
        inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
        wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
        p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
    else
        
        try
            nsig=3;
            ax(1).val=linspace(min(model.design_matrix_pred{jj}.eps_a_tone)-nsig*sqrt(sig_a_2),max(model.design_matrix_pred{jj}.eps_a_tone)+nsig*sqrt(sig_a_2),10);
            ax(2).val=linspace(min(model.design_matrix_pred{jj}.eps_v_right)-nsig*sqrt(sig_v_2),max(model.design_matrix_pred{jj}.eps_v_right)+nsig*sqrt(sig_v_2),10);
            Niteration=3;
            
            if jj>model.num_dur
                evalc('mdbm_sol=mdbm(ax,@(ax1) p_ch_x_exact(params,ax1(1,:),ax1(2,:),jj-model.num_dur),Niteration)');
            else
                evalc('mdbm_sol=mdbm(ax,@(ax1) p_ch_x_exact(params,ax1(1,:),ax1(2,:),jj),Niteration)');
            end
            
            
            xa=mdbm_sol.posinterp(1,:)';
            xv=mdbm_sol.posinterp(2,:)';
            [xv,idd]=sort(xv);
            xa=xa(idd);
            pa1=repmat(normcdf(max(ax(1).val),model.design_matrix_pred{jj}.eps_a_tone(:)',sqrt(sig_a_2)),size(xa,1),1)-normcdf(repmat(xa,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',size(xa,1),1),sqrt(sig_a_2));
            pa2=normcdf(repmat(xa,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',size(xa,1),1),sqrt(sig_a_2))-repmat(normcdf(min(ax(1).val),model.design_matrix_pred{jj}.eps_a_tone(:)',sqrt(sig_a_2)),size(xa,1),1);
            pv=normcdf(repmat(xv,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_v_right(:)',size(xv,1),1),sqrt(sig_v_2));
            
            a1=sum(0.5*(pa1(1:end-1,:)+pa1(2:end,:)).*abs((pv(1:end-1,:)-pv(2:end,:))),1)./abs(pv(1,:)-pv(end,:));
            a2=sum(0.5*(pa2(1:end-1,:)+pa2(2:end,:)).*abs((pv(1:end-1,:)-pv(2:end,:))),1)./abs(pv(1,:)-pv(end,:));
            p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,a1./(a1+a2));
            p_ch1{jj}(isnan(a1) | isnan(a2))=NaN;
        catch er
            X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
            X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
            if jj>model.num_dur
                tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj-model.num_dur);
            else
                tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj);
            end
            inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
            wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
            p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
        end
    end
    
    
end
end
function p_ch1 = simulate_prob_resp_pred_specific(params,eps_a,eps_v,dur, model)

params=get_full_params(params,model);



alp0=params(21);
d=params(22);
jj=dur;
lapse_rate=params(14+2*jj-1);
lapse_prob=params(14+2*jj);
kpa_2=params(2+jj);
kpv_2=params(5+jj);
sig_ap_2=params(10).^2;
sig_a_2=sig_ap_2./kpa_2;
sig_v_2=sig_a_2.*kpv_2;











% sig_rightleft_true_2=params(8+jj).^2;

num_inputs=numel(eps_a);

model.design_matrix_pred{jj}.eps_a_tone=map_stim_normstim(eps_a,alp0,d,model.eps_max);
    model.design_matrix_pred{jj}.eps_v_right=map_stim_normstim(eps_v,alp0,d,model.eps_max);





if model.exact_inference==0
    X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
    X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
    if jj>model.num_dur
        tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj-model.num_dur);
    else
        tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj);
    end
    inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
    wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
    p_ch1=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
else
    
    try
        nsig=3;
        ax(1).val=linspace(min(model.design_matrix_pred{jj}.eps_a_tone)-nsig*sqrt(sig_a_2),max(model.design_matrix_pred{jj}.eps_a_tone)+nsig*sqrt(sig_a_2),10);
        ax(2).val=linspace(min(model.design_matrix_pred{jj}.eps_v_right)-nsig*sqrt(sig_v_2),max(model.design_matrix_pred{jj}.eps_v_right)+nsig*sqrt(sig_v_2),10);
        Niteration=3;
        
        if jj>model.num_dur
            evalc('mdbm_sol=mdbm(ax,@(ax1) p_ch_x_exact(params,ax1(1,:),ax1(2,:),jj-model.num_dur),Niteration)');
        else
            evalc('mdbm_sol=mdbm(ax,@(ax1) p_ch_x_exact(params,ax1(1,:),ax1(2,:),jj),Niteration)');
        end
        
        
        xa=mdbm_sol.posinterp(1,:)';
        xv=mdbm_sol.posinterp(2,:)';
        [xv,idd]=sort(xv);
        xa=xa(idd);
        pa1=repmat(normcdf(max(ax(1).val),model.design_matrix_pred{jj}.eps_a_tone(:)',sqrt(sig_a_2)),size(xa,1),1)-normcdf(repmat(xa,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',size(xa,1),1),sqrt(sig_a_2));
        pa2=normcdf(repmat(xa,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',size(xa,1),1),sqrt(sig_a_2))-repmat(normcdf(min(ax(1).val),model.design_matrix_pred{jj}.eps_a_tone(:)',sqrt(sig_a_2)),size(xa,1),1);
        pv=normcdf(repmat(xv,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_v_right(:)',size(xv,1),1),sqrt(sig_v_2));
        
        a1=sum(0.5*(pa1(1:end-1,:)+pa1(2:end,:)).*abs((pv(1:end-1,:)-pv(2:end,:))),1)./abs(pv(1,:)-pv(end,:));
        a2=sum(0.5*(pa2(1:end-1,:)+pa2(2:end,:)).*abs((pv(1:end-1,:)-pv(2:end,:))),1)./abs(pv(1,:)-pv(end,:));
        p_ch1=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,a1./(a1+a2));
    catch
        X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
        X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
        if jj>model.num_dur
            tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj-model.num_dur);
        else
            tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj);
        end
        inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
        wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
        p_ch1=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
    end
end



end

function p_ch = p_ch_x(params,X_a_tonenoise,X_v_rightleft,jj)


num_inputs_quad=size(X_a_tonenoise,1);
num_inputs=size(X_a_tonenoise,2);


inds=cartprod([1:num_inputs_quad],[1:num_inputs_quad]);
numel_sz=size(inds,1);

eta_R=100;%params(2);
gamma=params(2);
kpa_2=params(2+jj);
kpv_2=params(5+jj);
sig_ap_2=params(10).^2;
mu_a=params(9)*sqrt(sig_ap_2);
sig_a_2=sig_ap_2./kpa_2;
beta=normcdf(mu_a./sqrt(sig_ap_2))+params(1);
beta(beta<0 | beta>1)=betarnd(1,1,sum(beta<0 | beta>1),1);
% beta=sigmoid(log(normcdf(mu_a./sqrt(sig_ap_2))./normcdf(-mu_a./sqrt(sig_ap_2)))+params(1));
sig_v_2=sig_a_2.*kpv_2;
sig_vp_2=exp(0.5*log(sig_ap_2)+params(12)).^2;
mu_v=(params(3)+params(11))*sqrt(sig_vp_2);
sig_avp_2=exp(0.5*log(sig_ap_2)+params(14)).^2;
mu_av=(params(3)+params(13))*sqrt(sig_avp_2);





alp_a=sig_ap_2./(sig_ap_2+sig_a_2);
sig_aap_2=sig_a_2.*alp_a;

alp_v=sig_vp_2./(sig_vp_2+sig_v_2);
sig_vvp_2=sig_v_2.*alp_v;

alp_av=sig_v_2./(sig_a_2+sig_v_2);
sig_av_2=sig_a_2.*alp_av;

alp_avavp=sig_avp_2./(sig_avp_2+sig_av_2);
sig_avavp_2=sig_av_2.*alp_avavp;



lp_xr_c0=zeros(numel_sz,num_inputs,2);
lp_xr_c1=zeros(numel_sz,num_inputs,2);
lp_xc=zeros(numel_sz,num_inputs,2);


X_a_tonenoise_toneprior=X_a_tonenoise(inds(:,1),:).*alp_a+mu_a.*(1-alp_a);

lp_xr_c0(:,:,1)=(lognormcdf(-1*X_a_tonenoise_toneprior,0,sqrt(sig_aap_2))...
    -lognormcdf(-1*mu_a,0,sqrt(sig_ap_2))...
    +log(1-beta)...
    );
lp_xr_c0(:,:,2)=(lognormcdf(1*X_a_tonenoise_toneprior,0,sqrt(sig_aap_2))...
    -lognormcdf(1*mu_a,0,sqrt(sig_ap_2))...
    +log(beta)...
    );



X_v_rightleft_rightprior=X_v_rightleft(inds(:,2),:).*alp_v+mu_v.*(1-alp_v);

% lp_x_c0=lognormcdf(X_v_rightleft_rightprior,0,sqrt(sig_vvp_2))...
%     -lognormcdf(mu_v,0,sqrt(sig_vp_2))...
%     +logsumexp(lp_xr_c0,3);

lp_x_c0=lognormcdf(X_v_rightleft_rightprior,0,sqrt(sig_vvp_2))...
    -lognormcdf(mu_v,0,sqrt(sig_vp_2))...
    +lognormpdf(X_a_tonenoise(inds(:,1),:),mu_a,sqrt(sig_a_2+sig_ap_2))...
    +lognormpdf(X_v_rightleft(inds(:,2),:),mu_v,sqrt(sig_v_2+sig_vp_2))...
    +logsumexp(lp_xr_c0,3);




X_av_m1=X_a_tonenoise(inds(:,1),:)*alp_av-X_v_rightleft(inds(:,2),:)*(1-alp_av);
X_av_1=X_a_tonenoise(inds(:,1),:)*alp_av+X_v_rightleft(inds(:,2),:)*(1-alp_av);


X_avavp_m1=X_av_m1*alp_avavp+mu_av*(1-alp_avavp);
X_avavp_1=X_av_1*alp_avavp+mu_av*(1-alp_avavp);


lp_xr_c1(:,:,1)=lognormcdf(-1*X_avavp_m1,0,sqrt(sig_avavp_2))...
    -lognormcdf(-1*mu_av,0,sqrt(sig_avp_2))...
    +lognormpdf(X_a_tonenoise(inds(:,1),:),-1*X_v_rightleft(inds(:,2),:),sqrt(sig_a_2+sig_v_2))...
    +lognormpdf(X_av_m1,mu_av,sqrt(sig_av_2+sig_avp_2))...
    +log(1-beta);


lp_xr_c1(:,:,2)=lognormcdf(X_avavp_1,0,sqrt(sig_avavp_2))...
    -lognormcdf(mu_av,0,sqrt(sig_avp_2))...
    +lognormpdf(X_a_tonenoise(inds(:,1),:),X_v_rightleft(inds(:,2),:),sqrt(sig_a_2+sig_v_2))...
    +lognormpdf(X_av_1,mu_av,sqrt(sig_av_2+sig_avp_2))...
    +log(beta);



lp_x_c1=squeeze(logsumexp(lp_xr_c1,3));


lmu_01=squeeze(lp_xr_c0(:,:,2))-squeeze(logsumexp(lp_xr_c0,3));

lmu_11=squeeze(lp_xr_c1(:,:,2))-squeeze(logsumexp(lp_xr_c1,3));

lp_xc(:,:,1)=log(1-gamma)+lp_x_c0;
lp_xc(:,:,2)=log(gamma)+lp_x_c1;



lmu_c0=squeeze(lp_xc(:,:,1))-squeeze(logsumexp(lp_xc,3));
lmu_c1=squeeze(lp_xc(:,:,2))-squeeze(logsumexp(lp_xc,3));

lp_xr_c1(:,:,1)=lmu_c0+lmu_01;
lp_xr_c1(:,:,2)=lmu_c1+lmu_11;

mu_R=exp(squeeze(logsumexp(lp_xr_c1,3)));

try
    p_ch = betainc(min(1,max(0,mu_R)),0.5*(eta_R+1),0.5*(eta_R+1));
    
catch er
%          keyboard
    'check3'
end




end



function mu_R = p_ch_x_exact(params,X_a_tonenoise,X_v_rightleft,jj)

if numel(X_a_tonenoise)>1 && numel(X_v_rightleft)==1
    sz=size(X_a_tonenoise);
    X_v_rightleft=X_v_rightleft*ones(sz);
elseif numel(X_v_rightleft)>1 && numel(X_a_tonenoise)==1
    sz=size(X_v_rightleft);
    X_a_tonenoise=X_a_tonenoise*ones(sz);
else
    sz=size(X_v_rightleft);
end



gamma=params(2);
kpa_2=params(2+jj);
kpv_2=params(5+jj);
sig_ap_2=params(10).^2;
mu_a=params(9)*sqrt(sig_ap_2);
sig_a_2=sig_ap_2./kpa_2;
beta=normcdf(mu_a./sqrt(sig_ap_2))+params(1);
beta(beta<0 | beta>1)=betarnd(1,1,sum(beta<0 | beta>1),1);
% beta=sigmoid(log(normcdf(mu_a./sqrt(sig_ap_2))./normcdf(-mu_a./sqrt(sig_ap_2)))+params(1));
sig_v_2=sig_a_2.*kpv_2;
sig_vp_2=exp(0.5*log(sig_ap_2)+params(12)).^2;
mu_v=(params(3)+params(11))*sqrt(sig_vp_2);
sig_avp_2=exp(0.5*log(sig_ap_2)+params(14)).^2;
mu_av=(params(3)+params(13))*sqrt(sig_avp_2);





alp_a=sig_ap_2./(sig_ap_2+sig_a_2);
sig_aap_2=sig_a_2.*alp_a;

alp_v=sig_vp_2./(sig_vp_2+sig_v_2);
sig_vvp_2=sig_v_2.*alp_v;

alp_av=sig_v_2./(sig_a_2+sig_v_2);
sig_av_2=sig_a_2.*alp_av;

alp_avavp=sig_avp_2./(sig_avp_2+sig_av_2);
sig_avavp_2=sig_av_2.*alp_avavp;



% lp_xr_c0=zeros(numel_sz,num_inputs,2);
% lp_xr_c1=zeros(numel_sz,num_inputs,2);
% lp_xc=zeros(numel_sz,num_inputs,2);


X_a_tonenoise_toneprior=X_a_tonenoise.*alp_a+mu_a.*(1-alp_a);

lp_xr_c0(:,:,1)=(lognormcdf(-1*X_a_tonenoise_toneprior,0,sqrt(sig_aap_2))...
    -lognormcdf(-1*mu_a,0,sqrt(sig_ap_2))...
    +log(1-beta)...
    );
lp_xr_c0(:,:,2)=(lognormcdf(1*X_a_tonenoise_toneprior,0,sqrt(sig_aap_2))...
    -lognormcdf(1*mu_a,0,sqrt(sig_ap_2))...
    +log(beta)...
    );



X_v_rightleft_rightprior=X_v_rightleft.*alp_v+mu_v.*(1-alp_v);

% lp_x_c0=lognormcdf(X_v_rightleft_rightprior,0,sqrt(sig_vvp_2))...
%     -lognormcdf(mu_v,0,sqrt(sig_vp_2))...
%     +logsumexp(lp_xr_c0,3);

lp_x_c0=lognormcdf(X_v_rightleft_rightprior,0,sqrt(sig_vvp_2))...
    -lognormcdf(mu_v,0,sqrt(sig_vp_2))...
    +lognormpdf(X_a_tonenoise,mu_a,sqrt(sig_a_2+sig_ap_2))...
    +lognormpdf(X_v_rightleft,mu_v,sqrt(sig_v_2+sig_vp_2))...
    +logsumexp(lp_xr_c0,3);




X_av_m1=X_a_tonenoise*alp_av-X_v_rightleft*(1-alp_av);
X_av_1=X_a_tonenoise*alp_av+X_v_rightleft*(1-alp_av);


X_avavp_m1=X_av_m1*alp_avavp+mu_av*(1-alp_avavp);
X_avavp_1=X_av_1*alp_avavp+mu_av*(1-alp_avavp);


lp_xr_c1(:,:,1)=lognormcdf(-1*X_avavp_m1,0,sqrt(sig_avavp_2))...
    -lognormcdf(-1*mu_av,0,sqrt(sig_avp_2))...
    +lognormpdf(X_a_tonenoise,-1*X_v_rightleft,sqrt(sig_a_2+sig_v_2))...
    +lognormpdf(X_av_m1,mu_av,sqrt(sig_av_2+sig_avp_2))...
    +log(1-beta);


lp_xr_c1(:,:,2)=lognormcdf(X_avavp_1,0,sqrt(sig_avavp_2))...
    -lognormcdf(mu_av,0,sqrt(sig_avp_2))...
    +lognormpdf(X_a_tonenoise,X_v_rightleft,sqrt(sig_a_2+sig_v_2))...
    +lognormpdf(X_av_1,mu_av,sqrt(sig_av_2+sig_avp_2))...
    +log(beta);



lp_x_c1=squeeze(logsumexp(lp_xr_c1,3));


lmu_01=squeeze(lp_xr_c0(:,:,2))-squeeze(logsumexp(lp_xr_c0,3));

lmu_11=squeeze(lp_xr_c1(:,:,2))-squeeze(logsumexp(lp_xr_c1,3));

lp_xc(:,:,1)=log(1-gamma)+lp_x_c0;
lp_xc(:,:,2)=log(gamma)+lp_x_c1;



lmu_c0=squeeze(lp_xc(:,:,1))-squeeze(logsumexp(lp_xc,3));
lmu_c1=squeeze(lp_xc(:,:,2))-squeeze(logsumexp(lp_xc,3));

lp_xr_c1(:,:,1)=lmu_c0+lmu_01;
lp_xr_c1(:,:,2)=lmu_c1+lmu_11;

mu_R=exp(squeeze(logsumexp(lp_xr_c1,3)))-0.5;


end