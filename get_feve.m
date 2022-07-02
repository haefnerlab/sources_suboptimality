function ve=get_feve(params_theta,model,nparam,compensate_bias)
data_mu=[];
data_var=[];
pred_mu=[];
nreps_data=[];

ch=model.simulate_prob_resp(params_theta,model);
for i=1:numel(model.data)
   pred_mu=[pred_mu;ch{i}(:)]; 
   alp=model.data{i}.num_ch1(:)+0.5; 
   bet=model.data{i}.num_repeats(:)-model.data{i}.num_ch1(:)+0.5;
   data_mu=[data_mu;alp./(alp+bet)];
   data_var=[data_var;(alp.*bet)./(((alp+bet).^2).*(alp+bet+1))];
   nreps_data=[nreps_data;model.data{i}.num_repeats(:)];
end
ve=feve2(data_mu,pred_mu,data_var,nreps_data,nparam,compensate_bias);

end