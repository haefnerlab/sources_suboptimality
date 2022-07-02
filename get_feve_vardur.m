function ve=get_feve_vardur(params_theta,model,nparam,compensate_bias)


data_mu=[];
data_var=[];
pred_mu=[];
nreps_data=[];
num_pts=10;
ch=model.simulate_prob_resp(params_theta,model);
for i=1:numel(model.data)
    
    [q1,q2,q3]=histcounts(model.design_matrix{i}.eps_a_tone,10);
    for k=1:num_pts
        
        alp(k)=sum(model.data{i}.num_ch1(q3==k))+0.5;
        bet(k)=sum(q3==k)-sum(model.data{i}.num_ch1(q3==k))+0.5;
        alp1(k)=sum(ch{i}(q3==k))+0.5;
        bet1(k)=sum(q3==k)-sum(ch{i}(q3==k))+0.5;
        nr(k)=sum(q3==k);
    end
    
    
    pred_mu=[pred_mu;alp1(:)./(alp1(:)+bet1(:))];
    alp=alp(:);
    bet=bet(:);
    data_mu=[data_mu;alp./(alp+bet)];
    data_var=[data_var;(alp.*bet)./(((alp+bet).^2).*(alp+bet+1))];
    nreps_data=[nreps_data;nr(:)];
end
ve=feve2(data_mu,pred_mu,data_var,nreps_data,nparam,compensate_bias);

end