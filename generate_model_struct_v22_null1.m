function model = generate_model_struct_v22_null1(id,incomplete_dir_save)
model.params_names = {...
    'Choice prior',...
    'Number of samples',...
    'Prior combination prob.',...
    'Aud. sensor noise (std)',...
    'Vis. sensor noise (std)',...
    'Aud. prior mean',...
    'Aud. prior std',...
    'Vis. prior mean',...
    'Vis. prior std',...
    'Lapse rate',...
    'Lapse prob',...
    'alp0',...
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

model.params_family = [3 , 1, 5, 5, 3,4,3,3,3,3    1,1,1,1];
model.params_additive_offset = [0,0,0,0 ,0,model.min_sig,0,0,0,0,       model.min_lapse,model.min_lapse,0,0.5];
model.params_multiplicative_offset = [1, 1, 1,1, 1,1,1,1,1,1,        1-model.min_lapse,1-model.min_lapse,1,1];
model.params_hyperprior_params{1} = [0,1.25,  1,1,   0,0,0,0,0,0,      1,1.25,1,2];
model.params_hyperprior_params{2} = [0.25,1.25,  1,1,   0.5,2.35,0.5*2,2.35*2,0.5*2,2.35*2,  5,1.25,1,2];%3.525,2];
% model.params_hyperprior_params{2} = [0.25,1,1.25,  1,1,   0.5,2.35,0.01,0.01,0.01,0.01,  5,1.25,3.525,2];

% mu chosen such that 95 percent mass between 0.05 and 0.95
% sig_pr chosen such that 95 percent mass between 1e-2 and 1e2



%quantile(lognrnd(-2.3,2,[1e5,1]),[0.125,0.875])=0.01,1


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


tp=0.5*[-20,-10,-5,-2.5,-1.25,1.25,2.5,5,10,20];
% tp=linspace(-10,10,100);
model.design_matrix{1}.eps_a_tone=tp;
model.design_matrix{1}.eps_v_right=0*model.design_matrix{1}.eps_a_tone;
model.design_matrix{2}.eps_a_tone=tp;
model.design_matrix{2}.eps_v_right=abs(model.design_matrix{2}.eps_a_tone);
model.design_matrix{1}.npts=length(model.design_matrix{1}.eps_a_tone);
model.design_matrix{2}.npts=length(model.design_matrix{1}.eps_a_tone);

model.design_matrix_pred{1}.eps_a_tone=0.5*linspace(-20,20,101);
model.design_matrix_pred{1}.eps_v_right=0*model.design_matrix_pred{1}.eps_a_tone;
model.design_matrix_pred{2}.eps_a_tone=0.5*linspace(-20,20,101);
model.design_matrix_pred{2}.eps_v_right=abs(model.design_matrix_pred{2}.eps_a_tone);
model.design_matrix_pred{1}.npts=length(model.design_matrix_pred{1}.eps_a_tone);
model.design_matrix_pred{2}.npts=length(model.design_matrix_pred{2}.eps_a_tone);

model.set_default = [7:10];
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

% model.plot_model_pred=@(model,figid,showlegend) plot_model_pred(model,figid,showlegend);
model.plot_model_pred = @(model, fig_handle, showlegend) plot_model_pred(model, fig_handle, showlegend);

model.get_samples_posterior = @(model, numchains, num_samples, num_samples_save) get_samples_posterior(model, numchains, num_samples, num_samples_save);
model.get_samples_posterior_parallel = @(model, numchains, num_samples, num_samples_save) get_samples_posterior_parallel(model, numchains, num_samples, num_samples_save);

model.log_unnorm_post = @(params, model) log_unnorm_post(params, model);
model.log_likelihood = @(params, model) log_likelihood(params, model);
model.log_prior_theta = @(theta, model) log_prior_theta(theta, model);
model.log_prior_phi = @(phi, model) log_prior_phi(phi, model);

model.get_map = @(model, nstart, gethessian, useparallel) get_map(model, nstart, gethessian, useparallel);
model.get_map_theta = @(model, nstart, gethessian, useparallel) get_map_theta(model, nstart, gethessian, useparallel);
model.get_mle = @(model, nstart, gethessian, useparallel) get_mle(model, nstart, gethessian, useparallel);
model.get_full_params= @(theta0, model) get_full_params(theta0, model);
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

function plot_model_pred(model, fig_handle, showlegend)
cols=[0,0,1;1,0,0];


% axes(fig_handle);





if isfield(model,'params_phi_mle')
    pt=model.phi_theta(model.params_phi_mle,model);
    ch1=simulate_prob_resp_pred(pt,model);
    
elseif isfield(model,'params_phi_map')
    pt=model.phi_theta(model.params_phi_map,model);
    ch1=simulate_prob_resp_pred(pt,model);
    
else
    pt=model.phi_theta(model.params_phi,model);
    ch1=simulate_prob_resp_pred(pt,model);
    
    
end

% subplot(2,2,1)
% hold on
% pl1(1)=plot(model.design_matrix_pred{1}.eps_a_tone,ch1{1},'linewidth',2);
%
% subplot(2,2,2)
% hold on
% pl1(2)=plot(model.design_matrix_pred{1}.eps_a_tone,ch1{2},'linewidth',2);
plot(fig_handle,[-10,10],[0.5,0.5],'k:','linewidth',1.5);
plot(fig_handle,[0,0],[0,1],'k:','linewidth',1.5);

for jj=1:numel(model.design_matrix)
    pl1(jj)=plot(fig_handle,model.design_matrix_pred{jj}.eps_a_tone,ch1{jj},'Color',cols(jj,:),'linewidth',2,'linestyle',':');
% pl1(jj)=plot(fig_handle,model.design_matrix_pred{1}.eps_a_tone,ch1{jj},'linewidth',2);
    hold on
end







% hline(0.5,'b--');
% vline(0,'b--');
xlim([min(model.design_matrix_pred{1}.eps_a_tone)-1,max(model.design_matrix_pred{1}.eps_a_tone)+1]);
ylim([0,1]);


if isfield(model,'data')
    for jj=1:numel(model.design_matrix)
        pr_mu=model.data{jj}.num_ch1./model.data{jj}.num_repeats;
        pr_se=sqrt((pr_mu.*(1-pr_mu))./model.data{jj}.num_repeats);
        pl2(jj)=errorbar(fig_handle,model.design_matrix{jj}.eps_a_tone,pr_mu,pr_se,'.','Color',cols(jj,:),'linewidth',2,'markersize',12,'capsize',0);
        %         pl2(jj)=errorbar(fig_handle,model.design_matrix{jj}.eps_a_tone,pr_mu,pr_se,'o','linewidth',2);
        
        
    end
    if showlegend==1
        legend([pl2,pl1],{'Central data','Matched data','Central pred.','Matched pred.'},'Location', 'best');
    end
else
    if showlegend==1
        legend(pl1,{'Central pred.','Matched pred.'},'Location', 'best');
    end
end



xlabel('tone position (deg)');
ylabel('prob. of ''right'' response');
set(gca,'fontsize',20,'fontweight','bold');
set(gcf,'color','white')
set(gca,'linewidth',4)
set(gca,'TickDir','out');

if showlegend==1
    legend boxoff
end
box off
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

options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e5,'MaxIterations',1e5);
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
% p_ch1=simulate_prob_resp(model.phi_theta(model.params_phi,model),model);

ll=0;
for jj=1:numel(model.design_matrix)
    if any(p_ch{jj}<=0) || any(p_ch{jj}>=1)
        warning('0 or 1 probability of choice');
        %         keyboard
        p_ch{jj}=min(1-eps,max(eps,p_ch{jj}));
    end
    
    ll=ll+sum((model.data{jj}.num_ch1.*log(p_ch{jj}))+((model.data{jj}.num_repeats-model.data{jj}.num_ch1).*log(1-p_ch{jj})),2);

%     model.data{jj}.num_repeats=1e6;
%     ll=ll+sum((model.data{jj}.num_repeats.*p_ch1{jj}.*log(p_ch{jj}))+(model.data{jj}.num_repeats.*(1-p_ch1{jj}).*log(1-p_ch{jj})),2);
%     ll=ll-sum(((p_ch1{jj}-p_ch{jj}).^2)./(p_ch1{jj}.*(1-p_ch1{jj})),2);
%     ll=ll+sum(lognormpdf(p_ch{jj},p_ch1{jj},sqrt(p_ch1{jj}.*(1-p_ch1{jj}))),2);

%     ll=ll+sum(exp(2*log(abs(p_ch1{jj}-p_ch{jj}))-log(p_ch{jj})-log(1-p_ch{jj})),2);
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



kpa_2=params(3);
kpv_2=params(4);
sig_ap_2=params(6).^2;
sig_a_2=sig_ap_2./kpa_2;
sig_v_2=sig_a_2.*kpv_2;
lapse_rate=params(11);
lapse_prob=params(12);
alp0=params(13);
d=params(14);



for jj=1:numel(model.design_matrix)
    
    
    model.design_matrix{jj}.eps_a_tone=map_stim_normstim(model.design_matrix{jj}.eps_a_tone,alp0,d,model.eps_max);
    model.design_matrix{jj}.eps_v_right=map_stim_normstim(model.design_matrix{jj}.eps_v_right,alp0,d,model.eps_max);
    
    
    num_inputs=model.design_matrix{jj}.npts;
    
    
    
    
    
    if model.exact_inference==0
        X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
        X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
        
        tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft);
        inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
        wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
        p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
    else
        try
            nsig=3;
            ax(1).val=linspace(min(model.design_matrix{jj}.eps_a_tone)-nsig*sqrt(sig_a_2),max(model.design_matrix{jj}.eps_a_tone)+nsig*sqrt(sig_a_2),20);
            ax(2).val=linspace(min(model.design_matrix{jj}.eps_v_right)-nsig*sqrt(sig_v_2),max(model.design_matrix{jj}.eps_v_right)+nsig*sqrt(sig_v_2),20);
            Niteration=3;
            evalc('mdbm_sol=mdbm(ax,@(ax1) p_ch_x_exact(params,ax1(1,:),ax1(2,:)),Niteration)');
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
        catch er
            X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
            X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
            
            tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft);
            inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
            wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
            p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
        end
        
        %         plot(mdbm_sol.posinterp(1,:),mdbm_sol.posinterp(2,:),'o-','linewidth',2)
    end
    
    
end
end
function p_ch1 = simulate_prob_resp_pred(params, model)

params=get_full_params(params,model);


kpa_2=params(3);
kpv_2=params(4);
sig_ap_2=params(6).^2;
sig_a_2=sig_ap_2./kpa_2;
sig_v_2=sig_a_2.*kpv_2;
lapse_rate=params(11);
lapse_prob=params(12);
alp0=params(13);
d=params(14);


for jj=1:numel(model.design_matrix_pred)
    
    
    model.design_matrix_pred{jj}.eps_a_tone=map_stim_normstim(model.design_matrix_pred{jj}.eps_a_tone,alp0,d,model.eps_max);
    model.design_matrix_pred{jj}.eps_v_right=map_stim_normstim(model.design_matrix_pred{jj}.eps_v_right,alp0,d,model.eps_max);
    
    num_inputs=model.design_matrix_pred{jj}.npts;
    
    
    
    
    
    if model.exact_inference==0
        X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
        X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
        
        tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft);
        inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
        wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
        p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
    else
        try
            nsig=3;
            ax(1).val=linspace(min(model.design_matrix_pred{jj}.eps_a_tone)-nsig*sqrt(sig_a_2),max(model.design_matrix_pred{jj}.eps_a_tone)+nsig*sqrt(sig_a_2),10);
            ax(2).val=linspace(min(model.design_matrix_pred{jj}.eps_v_right)-nsig*sqrt(sig_v_2),max(model.design_matrix_pred{jj}.eps_v_right)+nsig*sqrt(sig_v_2),10);
            Niteration=3;
            evalc('mdbm_sol=mdbm(ax,@(ax1) p_ch_x_exact(params,ax1(1,:),ax1(2,:)),Niteration)');
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
        catch er
            X_a_tonenoise=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_a_tone(:)',model.npts_quad,1),sqrt(sig_a_2));
            X_v_rightleft=norminv(repmat(0.5*model.pts_leg(:)+0.5,1,num_inputs),repmat(model.design_matrix_pred{jj}.eps_v_right(:)',model.npts_quad,1),sqrt(sig_v_2));
            
            tmp=p_ch_x(params,X_a_tonenoise,X_v_rightleft);
            inds=cartprod([1:model.npts_quad],[1:model.npts_quad]);
            wp=prod([model.wts_leg(inds(:,1)),model.wts_leg(inds(:,2))],2);
            p_ch1{jj}=(lapse_rate*lapse_prob)+(1-lapse_rate)*min(1,sum((0.5*0.5*repmat(wp,1,num_inputs)).*tmp,1));
        end
        %         plot(mdbm_sol.posinterp(1,:),mdbm_sol.posinterp(2,:),'o-','linewidth',2)
    end
    
    
end
end

function p_ch = p_ch_x(params,X_a_tonenoise,X_v_rightleft)

num_inputs_quad=size(X_a_tonenoise,1);
num_inputs=size(X_a_tonenoise,2);


inds=cartprod([1:num_inputs_quad],[1:num_inputs_quad]);
numel_sz=size(inds,1);

eta_R=100;%params(2);
gamma=params(2);
kpa_2=params(3);
kpv_2=params(4);
sig_ap_2=params(6).^2;
mu_a=params(5)*sqrt(sig_ap_2);
sig_a_2=sig_ap_2./kpa_2;
beta=normcdf(mu_a./sqrt(sig_ap_2))+params(1);
beta(beta<0 | beta>1)=betarnd(1,1,sum(beta<0 | beta>1),1);
% beta=sigmoid(log(normcdf(mu_a./sqrt(sig_ap_2))./normcdf(-mu_a./sqrt(sig_ap_2)))+params(1));
sig_v_2=sig_a_2.*kpv_2;
sig_vp_2=exp(0.5*log(sig_ap_2)+params(8)).^2;
mu_v=(params(5)+params(7))*sqrt(sig_vp_2);
sig_avp_2=exp(0.5*log(sig_ap_2)+params(10)).^2;
mu_av=(params(5)+params(9))*sqrt(sig_avp_2);





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



function mu_R = p_ch_x_exact(params,X_a_tonenoise,X_v_rightleft)

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
kpa_2=params(3);
kpv_2=params(4);
sig_ap_2=params(6).^2;
mu_a=params(5)*sqrt(sig_ap_2);
sig_a_2=sig_ap_2./kpa_2;
beta=normcdf(mu_a./sqrt(sig_ap_2))+params(1);
beta(beta<0 | beta>1)=betarnd(1,1,sum(beta<0 | beta>1),1);
% beta=sigmoid(log(normcdf(mu_a./sqrt(sig_ap_2))./normcdf(-mu_a./sqrt(sig_ap_2)))+params(1));
sig_v_2=sig_a_2.*kpv_2;
sig_vp_2=exp(0.5*log(sig_ap_2)+params(8)).^2;
mu_v=(params(5)+params(7))*sqrt(sig_vp_2);
sig_avp_2=exp(0.5*log(sig_ap_2)+params(10)).^2;
mu_av=(params(5)+params(9))*sqrt(sig_avp_2);




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



mu_R=exp(logsumexp(lp_xr_c1,3))-0.5;

end