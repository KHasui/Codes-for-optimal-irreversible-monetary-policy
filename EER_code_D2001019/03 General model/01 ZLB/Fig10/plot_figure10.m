% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.
close all
clear
clc
disp('Analyzing distribution of binding for degree of reversal aversions')
disp(' NOTE: This code takes at least 70 minutes to complete (if you stop, push "ctrl + c")')

rss    = 1.5/4;
rhog   = .6;
imin   = -rss;
varg   = (3.72/16)^2;
stdg   = 3.72/16;
nburn  = 1000; % # of burning periods (default setting is 200)
T      = 100000 + nburn;
npath  = 100; % # of paths
eg_mat = stdg*randn(npath,T);


lamdafs = 0:0.002:0.01; % grids of lamda_ir
nlf     = size(lamdafs,2);

loss              = zeros(npath,nlf);
brate             = zeros(npath,nlf);
bnum              = zeros(npath,nlf);
interval_zlb      = zeros(npath,nlf);
periods_zlbpolicy = zeros(npath,nlf);
periods_zlb_mat   = zeros(T,npath,nlf);
num_liftoff_mat   = zeros(npath,nlf);
num_ret_zlb_mat   = zeros(npath,nlf);

ind    = 1;
nind   = nlf*npath*T;
secnum = 20;
backspaces = '  ';

for ilf  = 1:nlf

    for ipath = 1:npath
        
        %initial value of state variables
        g_t      = 0;
        ilag_t   = 0;
        dltlag_t = 0;
        
        %initial value of state indeces
        ind_period         = 1;
        ind_bind_lag       = 0;
        ind_bindp_down_lag = 0;
        ind_rest_4q_lag    = 0;
        
        %
        bindrate         = zeros(T,1);
        dum_periods_zlb  = zeros(T,1);
        periods_nozlb    = zeros(T,1);
        y_impmat_irr     = zeros(T,1);
        p_impmat_irr     = zeros(T,1);
    
        eval(['load result_irrev_zlb_indlamdaf',num2str(ilf),'_sig100_kpp24 gmat ilagmat dltlagmat imat Fp Fy'])
    
        Fi = griddedInterpolant(gmat,ilagmat,dltlagmat,imat,'spline');
    
        for t = 1:T
            i_t                = max(imin,Fi(g_t,ilag_t,dltlag_t));
            p_impmat_irr(t,:)  = Fp(g_t,ilag_t,dltlag_t);
            y_impmat_irr(t,:)  = Fy(g_t,ilag_t,dltlag_t);
            
            ind_bind = i_t == imin;
            
            ind_liftoff = (ind_bind_lag-ind_bind)*ind_bind_lag;
            ind_rest_4q = (4+1)*ind_liftoff + ind_rest_4q_lag*(1-ind_liftoff);
            ind_in_4q   = ind_rest_4q >= 1;
            ind_ret_zlb = ind_bind*ind_in_4q*(1-ind_liftoff);
            
            dltlag_t       = round(sign(i_t-ilag_t));
            bindrate(t,:)  = ind_bind;
            ind_bindp_down = round(abs(ind_bind*(1-ind_bind_lag)));
            ind_bindp_up   = round(abs(ind_bind_lag*(1-ind_bind)));
                       
            ind_period = ind_period + round(ind_bindp_up);
            
            periods_zlb_mat(ind_period,ipath,ilf) = periods_zlb_mat(ind_period,ipath,ilf) + round(ind_bind);
            dum_periods_zlb(ind_period,:)         = round(ind_bind);
            periods_nozlb(ind_period,:)           = periods_nozlb(ind_period,:) + round(abs(1-ind_bind));
            num_liftoff_mat(ipath,ilf)            = num_liftoff_mat(ipath,ilf) + ind_liftoff;
            num_ret_zlb_mat(ipath,ilf)            = num_ret_zlb_mat(ipath,ilf) + ind_ret_zlb;
            
            %update state variables
            ilag_t = i_t;
            g_t    = rhog*g_t + eg_mat(ipath,t);
            ind_bind_lag       = ind_bind;
            ind_bindp_down_lag = ind_bindp_down;
            ind_rest_4q_lag = max(ind_rest_4q - 1,0)*(1-ind_bind);
            
                    % display progress
                    percentage = ind/nind*100;
                    percstr = sprintf('Iteration: %3.1f percent completed...',percentage);
                    fprintf([backspaces, percstr]);
                    backspaces = repmat(sprintf('\b'), 1, length(percstr));
            
                    ind = ind+1;
            
        end %for t = 1:T
    
        ndump = sum(dum_periods_zlb);
        nburn = 100; % # of burning periods
        periods_zlb  =  sum((periods_zlb_mat(:,ipath,ilf)))/ndump;
        periods_zlbpolicy(ipath,ilf) = periods_zlb;
        periods_nozlb  =  mean(periods_nozlb(1:ndump,:));
        interval_zlb(ipath,ilf) = periods_nozlb;

    end %for ipath = 1:npath

end %for ilf  = 1:nlf
disp(' iteration finished.')

average_zlb_periods = mean(interval_zlb);
[i1,i2] = find(isnan(periods_zlbpolicy)==1);
periods_zlbpolicy(i1,i2) = 0;
average_periods_zlbpolicy = mean(periods_zlbpolicy);

prob_liftoff_retzlb_mat = num_ret_zlb_mat./num_liftoff_mat;
prob_liftoff_retzlb = mean(prob_liftoff_retzlb_mat);

%
% PLOT FIGURE
% 

left_color  = [0,0,0];
right_color = [0,0,0];

lamdairs = lamdafs;


titlefontsize  = 22;
xlabelfontsize = 25;
legendfontsize = 14;
axisfontsize   = 15;

set(figure(1),'defaultAxesColorOrder',[left_color; right_color]);
subplot(1,2,1)
yyaxis left % duration
    bar(lamdairs,average_zlb_periods,'FaceColor',[0.5,0.6,1])
    
    h_axes = gca;
    h_axes.YAxis(1).FontSize = axisfontsize;
    h_axes.YAxis(1).FontName = 'Times';
    
    ylim([0 140]);
    xlim([-2e-003 12e-003]);
    xticks(lamdairs)
    ylabel('Number of periods','FontName','Times','FontSize',titlefontsize)
    
yyaxis right % interval
    h = plot(lamdairs,average_periods_zlbpolicy);
    
    h_axes = gca;
    h_axes.XAxis.FontSize = axisfontsize;
    h_axes.XAxis.FontName = 'Times';
    
    h_axes.YAxis(2).FontSize = axisfontsize;
    h_axes.YAxis(2).FontName = 'Times';
    
    set(h(1),'color','black','Marker','n','LineStyle','-','LineWidth',1.5)
    xlabel('$\lambda_{\rm ir}$','Interpreter','Latex','FontSize',xlabelfontsize);
    ylim([1 1.7]);
    xlim([-2e-003 12e-003]);
    xticks(lamdairs)
    ylabel('Number of periods','FontName','Times','FontSize',titlefontsize);
    legend('Average interval between ZIRPs (left scale)','Average duration of a ZIRP (right scale)','EdgeColor',[1 1 1],'Interpreter','Latex','FontSize',legendfontsize);
    title('  ','FontSize',35)
    
    annotation(figure(1),'textbox',[0.0753610899247043 0.935714285714293 0.043973544973545 0.0666666666666705],'LineStyle','none',...
        'String',{'(a)'},'FontSize',30,'FontName','Arial','FitBoxToText','off');
    
    
subplot(1,2,2)
    bar(lamdairs,prob_liftoff_retzlb,'FaceColor',[0.5,1,0.6])
    
    h_axes = gca;
    h_axes.XAxis.FontSize = axisfontsize;
    h_axes.XAxis.FontName = 'Times';
    
    h_axes.YAxis.FontSize = axisfontsize;
    h_axes.YAxis.FontName = 'Times';
    
    xlabel('$\lambda_{\rm ir}$','Interpreter','Latex','FontSize',xlabelfontsize);
    ylabel('Probability of returning to ZLB','FontName','Times','FontSize',titlefontsize)
    title('  ','FontSize',35)
    xlim([-2e-003 12e-003]);
    xticks(lamdairs)
    
    annotation(figure(1),'textbox',[0.511111112710083 0.93809523809525 0.043973544973545 0.0666666666666706],'LineStyle','none',...
        'String','(b)','FontSize',30,'FontName','Arial','FitBoxToText','off');
    
    
