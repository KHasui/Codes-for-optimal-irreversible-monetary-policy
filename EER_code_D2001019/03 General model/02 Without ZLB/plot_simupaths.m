% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.

close all
clear
clc
disp('Deriving simulated paths with policy reversal aversion:')

% 1. Generating pseudo demand shocks:
varg     = (3.72/16)^2;
rhog     = 0.6;
nshocks  = 13;
[eg,prb] = qnwnorm(nshocks,0,varg); % recall qnwnorm from CompEcon Toolbox
nburn    = 1000; % # of burning periods (default setting is 200)
T        = 10000 + nburn + 1;
eg_vec   = eg(discrand(T,prb));
Tupper   = T -1;

%
% Model 1: Policy reversal aversion (lambda_ir = 0.0016)
%

% prepare memories for control variables
i_path_pd   = zeros(1,T);
di_path_pd  = zeros(1,T);
dlt_path_pd = zeros(1,T);

i_path_ir10   = zeros(1,T);
di_path_ir10  = zeros(1,T);
dlt_path_ir10 = zeros(1,T);

i_path_ir0016   = zeros(1,T);
di_path_ir0016  = zeros(1,T);
dlt_path_ir0016 = zeros(1,T);

g_path      = zeros(1,T+1);


load result_irrev_indlamdaf1_zlb_sig100_kpp24 gmat ilagmat dltlagmat imat
Fi_pd = griddedInterpolant(gmat,ilagmat,dltlagmat,imat,'spline');

load result_irrev_indlamdaf2_zlb_sig100_kpp24 gmat ilagmat dltlagmat imat
Fi_0016 = griddedInterpolant(gmat,ilagmat,dltlagmat,imat,'spline');

load result_irrev_indlamdaf3_zlb_sig100_kpp24 gmat ilagmat dltlagmat imat rss
Fi_10 = griddedInterpolant(gmat,ilagmat,dltlagmat,imat,'spline');

% initial condition of state variables
dltlag_pd = 0;
ilag_pd   = 0;
dltlag_ir10 = 0;
ilag_ir10   = 0;
dltlag_ir0016 = 0;
ilag_ir0016   = 0;

g_t       = eg_vec(1);

for t = 1:Tupper
    
    % derive responses of i(t), pi(t), and x(t) using interporation
    i_pd   = max(-rss,Fi_pd(g_t,ilag_pd,dltlag_pd)); 
    di_pd  = i_pd-ilag_pd;          % delta
    dlt_pd = sign(round(di_pd,3)); % sign of delta
    
    i_ir10   = max(-rss,Fi_10(g_t,ilag_ir10,dltlag_ir10)); 
    di_ir10  = i_ir10-ilag_ir10;          % delta
    dlt_ir10 = sign(round(di_ir10,3)); % sign of delta
    
    i_ir0016   = max(-rss,Fi_0016(g_t,ilag_ir0016,dltlag_ir0016)); 
    di_ir0016  = i_ir0016-ilag_ir0016;          % delta
    dlt_ir0016 = sign(round(di_ir0016,3)); % sign of delta
    
    % insert impulse responses
    i_path_pd(t)   = i_pd;
    di_path_pd(t)  = di_pd;
    dlt_path_pd(t) = dlt_pd;
    
    i_path_ir10(t)   = i_ir10;
    di_path_ir10(t)  = di_ir10;
    dlt_path_ir10(t) = dlt_ir10;
    
    i_path_ir0016(t)   = i_ir0016;
    di_path_ir0016(t)  = di_ir0016;
    dlt_path_ir0016(t) = dlt_ir0016;
    
    g_path(t)   = g_t;
    
    % update state variable
    ilag_pd   = i_pd;
    dltlag_pd = dlt_pd;
    
    ilag_ir10   = i_ir10;
    dltlag_ir10 = dlt_ir10;
    
    ilag_ir0016   = i_ir0016;
    dltlag_ir0016 = dlt_ir0016;
    
    g_t      = rhog*g_t + eg_vec(t+1);
    
end

ipaths_pd     = i_path_pd(nburn+1:T) + rss;
ipaths_ir10   = i_path_ir10(nburn+1:T) + rss;
ipaths_ir0016 = i_path_ir0016(nburn+1:T) + rss;
tgrids        = 0:size(ipaths_pd,2)-1;

gfonts = 15;
xfonts = 22;
yfonts = 22;

figure(1)
    h = plot(tgrids,4*ipaths_pd,tgrids,4*ipaths_ir10,tgrids,4*ipaths_ir0016,tgrids,1.5*ones(size(tgrids)));
    xlim([0 40])
    
    h_axes = gca;
    h_axes.XAxis.FontSize = gfonts;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = gfonts;
    h_axes.YAxis.FontName = 'Times';
    
    xlabel('Quarters','FontName','Times','Fontsize',xfonts)
    ylabel('Interest rate','FontName','Times','Fontsize',yfonts)
    
    set(h(1),'color','black','LineWidth',1,'LineStyle','-')
    set(h(2),'color','red','LineWidth',1,'LineStyle','--')
    set(h(3),'color','blue','LineWidth',2,'LineStyle',':')
    set(h(4),'color','black','LineWidth',1,'LineStyle','--')
    
    legend('Pure discretion','$\lambda_{\rm ir} = 10$','$\lambda_{\rm ir} = 0.0016$','Steady State','Interpreter','Latex','EdgeColor',[1,1,1])





