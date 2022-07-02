cc=cbrewer('div','RdBu',11);
cc=cc([end,end-1, end-2],:);

emax=atand((0.5*60.96)/50);

mu_deg=[-10,-5,-2.5,0,2.5,5,10];



mu_internal_linear=mu_deg./emax;
std_internal=0.05*ones(size(mu_internal_linear));
std_deg_linear=(emax).*std_internal;



m1=map_stim_normstim(mu_deg,0,1,emax);
ms=sign(m1);
m1=abs(m1);
m=log(1+emax)*m1;

s_2=(log(1+emax).*std_internal).^2;
mu_internal_log=(exp(m+0.5*s_2)-1).*ms;
std_deg_log=sqrt((exp(s_2)-1).*exp(2*m+s_2));
m1=m1.*ms;



m2=map_stim_normstim(mu_deg,0.5,1,emax);


for j=1:numel(m2)
j
x=m2(j)+normrnd(0,std_internal(j),1e3,1);
for i=1:1e3
y(i)=map_normstim_stim(x(i),0.5,1,emax);
end
y=y(:);
std_deg_inter(j)=std(y);
end




errorbar(mu_deg,mu_internal_linear,std_internal,std_internal,std_deg_linear,std_deg_linear,'.-','linewidth',1.5,'markersize',20,'color',cc(3,:));
hold on
errorbar(mu_deg,m2,std_internal,std_internal,std_deg_inter,std_deg_inter,'.-','linewidth',1.5,'markersize',20,'color',cc(2,:));
errorbar(mu_deg,m1,std_internal,std_internal,std_deg_log,std_deg_log,'.-','linewidth',1.5,'markersize',20,'color',cc(1,:));


axis square
set(gca,'fontsize',60,'fontweight','bold')
set(gca,'linewidth',8)
box off
set(gca,'TickDir','out')
set(gcf,'color','white')
ylim([-1,1]);
xlim([-15,15])
xticklabels([])
yticklabels([])
