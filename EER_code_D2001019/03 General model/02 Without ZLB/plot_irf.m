% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.

close all
clear
clc
disp('Deriving impulse responses with policy reversal aversion:')

z = 0;  % Temporal shock for 1 period
%z = 1;  % Temporal shock for 2 periods

shock  = -3/4; % shock valuen (quarterly)
Tupper = 6;    % for simulated periods 
eg_vec = zeros(Tupper+2,1);
for j = 2:2+z
    eg_vec(j) = shock;
end

%
% Model 1: Policy reversal aversion (lambda_ir = 0.0016)
%

lamdaf_irr = 0.0016;
lamdafs    = [0,0.0016,10];

indlamdaf = find(lamdafs == lamdaf_irr);

eval(['load result_irrev_indlamdaf',num2str(indlamdaf),'_nozlb_sig100_kpp24'])


% Derive Simulated Paths

% initial condition of state variables
dltlag_t = 0;
ilag_t   = 0;
g_t      = eg_vec(1);

% prepare memories for control variables
i_imp_irr = zeros(1,Tupper+1);
p_imp_irr = zeros(1,Tupper+1);
y_imp_irr = zeros(1,Tupper+1);
r_imp_irr = zeros(1,Tupper+1);
g_imp_irr = zeros(1,Tupper+1);
di_imp    = zeros(1,Tupper+1);
dlt_imp   = zeros(1,Tupper+1);

Fi = griddedInterpolant(gmat,ilagmat,dltlagmat,imat,'spline');

for t = 1:Tupper+1
    
    % derive responses of i(t), pi(t), and x(t) using interporation
    i_t   = max(-rss,Fi(g_t,ilag_t,dltlag_t));
    p_t   = Fp(g_t,ilag_t,dltlag_t);
    y_t   = Fy(g_t,ilag_t,dltlag_t);
    
    ep_t  = (p_t-kpp*y_t)/bet; % expected inflation
    r_t   = i_t - ep_t;        % real interest rate
    
    di_t  = i_t-ilag_t;          % delta
    dlt_t = sign(round(di_t,3)); % sign of delta
    
    % insert impulse responses
    p_imp_irr(t)   = p_t;
    y_imp_irr(t)   = y_t;
    i_imp_irr(t)   = i_t;
    r_imp_irr(t)   = r_t;
    g_imp_irr(t)   = g_t;
    di_imp(t)      = di_t;
    dlt_imp(t)     = dlt_t;
    
    % update state variable
    ilag_t   = i_t;
    dltlag_t = dlt_t;
    g_t      = eg_vec(t+1);
    
end


%
% Model 2: Pure discretion
%

lamdaf_pd = 0;
indlamdaf = find(lamdafs == lamdaf_pd);

eval(['load result_irrev_indlamdaf',num2str(indlamdaf),'_nozlb_sig100_kpp24'])


% Derive Simulated Paths

% initial condition of state variables
dltlag_t_pd = 0;
ilag_t_pd   = 0;
g_t_pd      = eg_vec(1);

% prepare memories for control variables
i_imp_pd   = zeros(1,Tupper+1);
p_imp_pd   = zeros(1,Tupper+1);
y_imp_pd   = zeros(1,Tupper+1);
r_imp_pd   = zeros(1,Tupper+1);
g_imp_pd   = zeros(1,Tupper+1);
di_imp_pd  = zeros(1,Tupper+1);
dlt_imp_pd = zeros(1,Tupper+1);

Fi = griddedInterpolant(gmat,ilagmat,dltlagmat,imat,'spline');

for t = 1:Tupper+1
    
    % derive responses of i(t), pi(t), and x(t) using interporation
    i_t_pd   = max(-rss+0.01/4,Fi(g_t_pd,ilag_t_pd,dltlag_t_pd));
    p_t_pd   = Fp(g_t_pd,ilag_t_pd,dltlag_t_pd);
    y_t_pd   = Fy(g_t_pd,ilag_t_pd,dltlag_t_pd);
    
    ep_t_pd  = (p_t_pd-kpp*y_t_pd)/bet; % expected inflation
    r_t_pd   = i_t_pd - ep_t_pd;        % real interest rate
    
    di_t_pd  = i_t_pd-ilag_t_pd;          % delta
    dlt_t_pd = sign(round(di_t_pd,3)); % sign of delta
    
    % insert impulse responses
    p_imp_pd(t)   = p_t_pd;
    y_imp_pd(t)   = y_t_pd;
    i_imp_pd(t)   = i_t_pd;
    r_imp_pd(t)   = r_t_pd;
    g_imp_pd(t)   = g_t_pd;
    di_imp_pd(t)  = di_t_pd;
    dlt_imp_pd(t) = dlt_t_pd;
    
    % update state variable
    ilag_t_pd   = i_t_pd;
    dltlag_t_pd = dlt_t_pd;
    g_t_pd      = eg_vec(t+1);
    
end


% 3: PLOT FIGURE:

tgrid = (0:Tupper)'; % time grid
    
MarkerSize_irr = 6; 
MarkerSize_pd  = 6; 
TitleFontSize  = 15;
    
figure(1)

subplot(2,2,1) % output gap
    h=plot(tgrid,y_imp_pd,tgrid,y_imp_irr,tgrid,zeros(size(tgrid)));
    set(h(1),'color','black','linewidth',1,'Marker','^','MarkerSize',MarkerSize_pd);
    set(h(2),'color','blue','linewidth',1,'Marker','o','MarkerSize',MarkerSize_irr);
    set(h(3),'color','black','linestyle','--');
    legend('Pure discretion',['$\lambda_{\rm ir} = ',num2str(round(lamdaf_irr,4)),'$'],'Steady State','Interpreter','Latex','fontsize',11,'EdgeColor','n')
    xticks(tgrid)
    title('(a) Output gap','fontsize',TitleFontSize+2,'FontName','Times')
    
    h_axes = gca;
    h_axes.XAxis.FontSize = TitleFontSize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = TitleFontSize;
    h_axes.YAxis.FontName = 'Times';
    
subplot(2,2,2) % inflation rate
    h=plot(tgrid,4*p_imp_pd,tgrid,4*p_imp_irr,tgrid,zeros(size(tgrid)));
    set(h(1),'color','black','linewidth',1,'Marker','^','MarkerSize',MarkerSize_pd);
    set(h(2),'color','blue','linewidth',1,'Marker','o','MarkerSize',MarkerSize_irr);
    set(h(3),'color','black','linestyle','--');
    title('(b) Inflation rate','fontsize',TitleFontSize+2,'FontName','Times')
    xticks(tgrid)
    
    h_axes = gca;
    h_axes.XAxis.FontSize = TitleFontSize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = TitleFontSize;
    h_axes.YAxis.FontName = 'Times';
    
subplot(2,2,3) % Nominal interest rate
    h=plot(tgrid,4*i_imp_pd+4*rss,tgrid,4*i_imp_irr+4*rss,tgrid,4*rss*ones(size(tgrid)));
    set(h(1),'color','black','linewidth',1,'Marker','^','MarkerSize',MarkerSize_pd);
    set(h(2),'color','blue','linewidth',1,'Marker','o','MarkerSize',MarkerSize_irr);
    set(h(3),'color','black','linestyle','--');
    title('(c) Nominal interest rate','fontsize',TitleFontSize+2,'FontName','Times')
    ylim([0 2])
    xticks(tgrid)
    
    h_axes = gca;
    h_axes.XAxis.FontSize = TitleFontSize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = TitleFontSize;
    h_axes.YAxis.FontName = 'Times';
    
subplot(2,2,4) % Real interest rate
    h=plot(tgrid,4*g_imp_irr+4*rss,tgrid,4*r_imp_pd+4*rss,tgrid,4*r_imp_irr+4*rss,tgrid,4*rss*ones(size(tgrid)));
    set(h(1),'color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',1.8)
    set(h(2),'color','black','linewidth',1,'Marker','^','MarkerSize',MarkerSize_pd);
    set(h(3),'color','blue','linewidth',1,'Marker','o','MarkerSize',MarkerSize_irr);
    set(h(4),'color','black','linestyle','--');
    title('(d) Real interest rate','fontsize',TitleFontSize+2,'FontName','Times')
    legend(h(1),'Natural rate','Interpreter','Latex','fontsize',11,'EdgeColor','n')
    xticks(tgrid)
    
    h_axes = gca;
    h_axes.XAxis.FontSize = TitleFontSize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = TitleFontSize;
    h_axes.YAxis.FontName = 'Times';
    
    annotation(figure(1),'textbox',[0.451 0.0190476190476191 0.106142857142857 0.0500000000000002],...
    'String',{'Quarters'},'LineStyle','none','FontSize',TitleFontSize,'FontName','Times New Roman','FitBoxToText','off');