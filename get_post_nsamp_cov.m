function post=get_post_nsamp_cov(mu,var,bins)


post=zeros(1,numel(bins)-1);
for i=1:numel(post)
    
    post(i)=normcdf(norminv((bins(i+1)-1)/99,0,1),mu,sqrt(var))-normcdf(norminv((bins(i)-1)/99,0,1),mu,sqrt(var));
    
    
end



end