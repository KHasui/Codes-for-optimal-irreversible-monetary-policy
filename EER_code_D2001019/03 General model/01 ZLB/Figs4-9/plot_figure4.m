% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.

close all
clear

%
% 1. Load result
%

% interpolation nodes

rhog = .6;
gmax = 4*sqrt(0.7^2/(1-rhog^2));
gmin = -gmax;
imax = 8;
rss  = 1.5/4;
sig    = 1;

gridg         = linspace(gmin,gmax,55);
gridilag      = 0*ones(size(gridg));
griddltlagneg = -1*ones(size(gridg));
griddltlag0   = 0*ones(size(gridg));
griddltlag1   = 1*ones(size(gridg));
lamdafs       = [0,0.0016,10];

% Load pure discretion
lamda_ir = 0;
lind = find(lamdafs == lamda_ir);
eval(['load result_irrev_indlamdaf',num2str(lind),'_zlb_sig100_kpp24 gmat ilagmat dltlagmat imat pmat ymat'])
Fipd = griddedInterpolant(gmat,ilagmat,dltlagmat,imat);
Fppd = griddedInterpolant(gmat,ilagmat,dltlagmat,pmat);
Fypd = griddedInterpolant(gmat,ilagmat,dltlagmat,ymat);
clear gmat ilagmat dltlagmat imat pmat ymat

% Load policy irreversibility (lambda_ir = 10)
lamda_ir10 = 10;
lind = find(lamdafs == lamda_ir10);
eval(['load result_irrev_indlamdaf',num2str(lind),'_zlb_sig100_kpp24 gmat ilagmat dltlagmat imat pmat ymat'])
Fi10 = griddedInterpolant(gmat,ilagmat,dltlagmat,imat,'spline');
Fp10 = griddedInterpolant(gmat,ilagmat,dltlagmat,pmat,'spline');
Fy10 = griddedInterpolant(gmat,ilagmat,dltlagmat,ymat,'spline');
clear gmat ilagmat dltlagmat imat pmat ymat

% Load policy irreversibility (lambda_ir = 0.0016)
lamda_ir0016 = 0.0016;
lind = find(lamdafs == lamda_ir0016);
eval(['load result_irrev_indlamdaf',num2str(lind),'_zlb_sig100_kpp24 gmat ilagmat dltlagmat imat pmat ymat'])
Fi016 = griddedInterpolant(gmat,ilagmat,dltlagmat,imat);
Fp016 = griddedInterpolant(gmat,ilagmat,dltlagmat,pmat);
Fy016 = griddedInterpolant(gmat,ilagmat,dltlagmat,ymat);
%clear gmat ilagmat dltlagmat imat pmat ymat

% check extrapolation for nominal rate
delta_i_upper = imax - max(max(max(imat)));
if delta_i_upper < 0
    warning('Extrapolation is dectected for nominal interest rate!')
end
%min(min(min(imat)))

%
% 2. derive policy function of nominal interest rate (interpolation)
%

% pure discretion
i_intp_pd = 4*max(-rss,Fipd(gridg,gridilag,griddltlag0));
p_intp_pd = 4*Fppd(gridg,gridilag,griddltlag0);
y_intp_pd = Fypd(gridg,gridilag,griddltlag0);

% lamda_ir = 0.0016
%   \delta_{t-1} = -1
i016_intpneg = 4*max(-rss,Fi016(gridg,gridilag,griddltlagneg));
p016_intpneg = 4*Fp016(gridg,gridilag,griddltlagneg);
y016_intpneg = Fy016(gridg,gridilag,griddltlagneg);
%   \delta_{t-1} = 0
i016_intp0   = 4*max(-rss,Fi016(gridg,gridilag,griddltlag0));
p016_intp0   = 4*Fp016(gridg,gridilag,griddltlag0);
y016_intp0   = Fy016(gridg,gridilag,griddltlag0);
%   \delta_{t-1} = 1
i016_intp1   = 4*max(-rss,Fi016(gridg,gridilag,griddltlag1));
p016_intp1   = 4*Fp016(gridg,gridilag,griddltlag1);
y016_intp1   = Fy016(gridg,gridilag,griddltlag1);


% lamda_ir = 10
%   \delta_{t-1} = -1
i10_intpneg = 4*max(-rss,Fi10(gridg,gridilag,griddltlagneg));
p10_intpneg = 4*Fp10(gridg,gridilag,griddltlagneg);
y10_intpneg = Fy10(gridg,gridilag,griddltlagneg);
%   \delta_{t-1} = 0
i10_intp0   = 4*max(-rss,Fi10(gridg,gridilag,griddltlag0));
p10_intp0   = 4*Fp10(gridg,gridilag,griddltlag0);
y10_intp0   = Fy10(gridg,gridilag,griddltlag0);
%   \delta_{t-1} = 1
i10_intp1   = 4*max(-rss,Fi10(gridg,gridilag,griddltlag1));
p10_intp1   = 4*Fp10(gridg,gridilag,griddltlag1);
y10_intp1   = Fy10(gridg,gridilag,griddltlag1);


%
% 3. plot figures
%

gridg = 4*(gridg/sig + rss) ;

xfont = 20;
yfont = 25;
tfont = 30;
gridfontsize = 19;
lsize = 1.3;


figure(1)

%(a) lambda_ir = 10      
subplot(1,2,1)
h=plot(gridg,i_intp_pd+4*rss,gridg,i10_intpneg+4*rss,gridg,i10_intp1+4*rss,gridg,i10_intp0+4*rss,[4*rss 4*rss],[-100 100]);
    h_axes = gca;
    h_axes.XAxis.FontSize = gridfontsize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = gridfontsize;
    h_axes.YAxis.FontName = 'Times';

    legend([h(2) h(4) h(3) h(1) h(5)],...
        '$\delta_{t-1} = -1$','$\delta_{t-1} = 0$','$\delta_{t-1} = 1$','Pure discretion','$4\times i^*$',...
        'Interpreter','Latex','EdgeColor',[1 1 1],...
        'Position',[0.141043572077454 0.588809522844497 0.134901884168141 0.242857146263123],'FontSize',16);
    %legend([h(2),h(1)],'a','b')
    ylabel('Nominal interest rate','FontName','Times','FontSize',yfont)
    
    set(h(1),'color','black','linestyle','-','marker','n','linewidth',lsize)
    set(h(2),'color','blue','linestyle','n','marker','o','linewidth',lsize,'MarkerSize',10)
    set(h(3),'color','red','linestyle','n','marker','d','linewidth',lsize,'MarkerSize',10)
    set(h(4),'color',[0,.5,0],'linestyle','--','marker','n','linewidth',2,'MarkerSize',10)
    set(h(5),'color','black','linestyle',':','marker','n','linewidth',lsize)
    xlim([-6+4*rss 6+4*rss]);xticks(-4:1:8);
    ylim([-1.5 9])%;yticks([-1.5:2:0,1.5,3:2:9])
    title(' ','FontName','Times','FontSize',tfont)


     annotation(figure(1),'textbox',...
    [0.083117247864793 0.864285714285717 0.0441189499589827 0.119047619047625],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontSize',30,...
    'FontName','Arial',...
    'FitBoxToText','off');

    annotation(figure(1),'textbox',...
    [0.122145001978779 0.873809523809525 0.350284175642088 0.104761904761907],...
    'String',{['$\lambda_{\rm ir} =',num2str(lamda_ir10),'$']},...
    'LineStyle','none',...
    'Interpreter','Latex',...
    'FontSize',30,...
    'FontName','.AppleSystemUIFont',...
    'FitBoxToText','off');

%(a) lambda_ir = 0.0016
subplot(1,2,2)
h=plot(gridg,i_intp_pd+4*rss,gridg,i016_intpneg+4*rss,gridg,i016_intp1+4*rss,gridg,i016_intp0+4*rss,[4*rss 4*rss],[-100 100]);
    h_axes = gca;
    h_axes.XAxis.FontSize = gridfontsize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = gridfontsize;
    h_axes.YAxis.FontName = 'Times';
    
    xlabel(' ','FontName','Times','FontSize',40)
    title('@','FontName','Times','FontSize',40)
    
    annotation(figure(1),'textbox',...
    [0.522649606793412 0.86428571428572 0.0441189499589826 0.119047619047624],...
    'String','(b)',...
    'LineStyle','none',...
    'FontSize',30,...
    'FontName','Arial',...
    'FitBoxToText','off');
    
    annotation(figure(1),'textbox',...
    [0.561996350161188 0.873809523809529 0.350284175642088 0.104761904761907],...
    'String',{['$\lambda_{\rm ir} =',num2str(lamda_ir0016),'$']},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',30,...
    'FontName','.AppleSystemUIFont',...
    'FitBoxToText','off');
    

    annotation(figure(1),'textbox',...
    [0.430184549356223 0.0333333333333334 0.340250233252473 0.0666666666666671],...
    'String',{'Natural interest rate'},...
    'LineStyle','none',...
    'FontSize',25,...
    'FontName','Times',...
    'FitBoxToText','off');

    
    set(h(1),'color','black','linestyle','-','marker','n','linewidth',lsize)
    set(h(2),'color','blue','linestyle','n','marker','o','linewidth',lsize,'MarkerSize',10)
    set(h(3),'color','red','linestyle','n','marker','d','linewidth',lsize,'MarkerSize',10)
    set(h(4),'color',[0,.5,0],'linestyle','--','marker','n','linewidth',2,'MarkerSize',10)
    set(h(5),'color','black','linestyle',':','marker','n','linewidth',lsize)
    xlim([-6+4*rss 6+4*rss]);xticks(-4:1:8);
    ylim([-1.5 9])%;yticks([-1.5:2:0,1.5,3:2:9])
 
    
    
    
    
figure(2)

%(a) lambda_ir = 10      
subplot(2,2,1)
h=plot(gridg,p_intp_pd,gridg,p10_intpneg,gridg,p10_intp1,gridg,p10_intp0,[4*rss 4*rss],[-100 100]);
    h_axes = gca;
    h_axes.XAxis.FontSize = gridfontsize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = gridfontsize;
    h_axes.YAxis.FontName = 'Times';

    legend([h(2) h(4) h(3) h(1) h(5)],...
        '$\delta_{t-1} = -1$','$\delta_{t-1} = 0$','$\delta_{t-1} = 1$','Pure discretion','$4\times i^*$',...
        'Interpreter','Latex','EdgeColor',[1 1 1],...
        'Position',[0.311731431688016 0.5923207901334 0.148364874624437 0.156202146141672],'FontSize',16);
    %legend([h(2),h(1)],'a','b')
    ylabel('Inflation rate','FontName','Times','FontSize',yfont)
    
    set(h(1),'color','black','linestyle','-','marker','n','linewidth',lsize)
    set(h(2),'color','blue','linestyle','n','marker','o','linewidth',lsize,'MarkerSize',10)
    set(h(3),'color','red','linestyle','n','marker','d','linewidth',lsize,'MarkerSize',10)
    set(h(4),'color',[0,.5,0],'linestyle','--','marker','n','linewidth',2,'MarkerSize',10)
    set(h(5),'color','black','linestyle',':','marker','n','linewidth',lsize)
    xlim([-6+4*rss 6+4*rss]);xticks(-4:1:8);
    ylim([-.4 .2])%;yticks([-1.5:2:0,1.5,3:2:9])
    title(' ','FontName','Times','FontSize',tfont)


     annotation(figure(2),'textbox',...
    [0.083117247864793 0.864285714285717 0.0441189499589827 0.119047619047625],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontSize',30,...
    'FontName','Arial',...
    'FitBoxToText','off');

    annotation(figure(2),'textbox',...
    [0.122145001978779 0.873809523809525 0.350284175642088 0.104761904761907],...
    'String',{['$\lambda_{\rm ir} =',num2str(lamda_ir10),'$']},...
    'LineStyle','none',...
    'Interpreter','Latex',...
    'FontSize',30,...
    'FontName','.AppleSystemUIFont',...
    'FitBoxToText','off');

%(a) lambda_ir = 0.0016
subplot(2,2,2)
h=plot(gridg,p_intp_pd,gridg,p016_intpneg,gridg,p016_intp1,gridg,p016_intp0,[4*rss 4*rss],[-100 100]);
    h_axes = gca;
    h_axes.XAxis.FontSize = gridfontsize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = gridfontsize;
    h_axes.YAxis.FontName = 'Times';
    
    title('@','FontName','Times','FontSize',40)
    
    annotation(figure(2),'textbox',...
    [0.522649606793412 0.86428571428572 0.0441189499589826 0.119047619047624],...
    'String','(b)',...
    'LineStyle','none',...
    'FontSize',30,...
    'FontName','Arial',...
    'FitBoxToText','off');
    
    annotation(figure(2),'textbox',...
    [0.561996350161188 0.873809523809529 0.350284175642088 0.104761904761907],...
    'String',{['$\lambda_{\rm ir} =',num2str(lamda_ir0016),'$']},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',30,...
    'FontName','.AppleSystemUIFont',...
    'FitBoxToText','off');
    

    annotation(figure(2),'textbox',...
    [0.414055517098159 0.00331588132635246 0.340250233252473 0.0666666666666671],...
    'String',{'Natural interest rate'},...
    'LineStyle','none',...
    'FontSize',25,...
    'FontName','Times',...
    'FitBoxToText','off');

    
    set(h(1),'color','black','linestyle','-','marker','n','linewidth',lsize)
    set(h(2),'color','blue','linestyle','n','marker','o','linewidth',lsize,'MarkerSize',10)
    set(h(3),'color','red','linestyle','n','marker','d','linewidth',lsize,'MarkerSize',10)
    set(h(4),'color',[0,.5,0],'linestyle','--','marker','n','linewidth',2,'MarkerSize',10)
    set(h(5),'color','black','linestyle',':','marker','n','linewidth',lsize)
    xlim([-6+4*rss 6+4*rss]);xticks(-4:1:8);
    ylim([-.4 .2])%;yticks([-1.5:2:0,1.5,3:2:9])
    

    
    
 subplot(2,2,3)
h=plot(gridg,y_intp_pd,gridg,y10_intpneg,gridg,y10_intp1,gridg,y10_intp0,[4*rss 4*rss],[-100 100]);
    h_axes = gca;
    h_axes.XAxis.FontSize = gridfontsize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = gridfontsize;
    h_axes.YAxis.FontName = 'Times';

    ylabel('Output gap','FontName','Times','FontSize',yfont)
    
    set(h(1),'color','black','linestyle','-','marker','n','linewidth',lsize)
    set(h(2),'color','blue','linestyle','n','marker','o','linewidth',lsize,'MarkerSize',10)
    set(h(3),'color','red','linestyle','n','marker','d','linewidth',lsize,'MarkerSize',10)
    set(h(4),'color',[0,.5,0],'linestyle','--','marker','n','linewidth',2,'MarkerSize',10)
    set(h(5),'color','black','linestyle',':','marker','n','linewidth',lsize)
    xlim([-6+4*rss 6+4*rss]);xticks(-4:1:8);
    ylim([-3 2])
    title(' ','FontName','Times','FontSize',tfont)



%(a) lambda_ir = 0.0016
subplot(2,2,4)
h=plot(gridg,y_intp_pd,gridg,y016_intpneg,gridg,y016_intp1,gridg,y016_intp0,[4*rss 4*rss],[-100 100]);
    h_axes = gca;
    h_axes.XAxis.FontSize = gridfontsize;
    h_axes.XAxis.FontName = 'Times';
    h_axes.YAxis.FontSize = gridfontsize;
    h_axes.YAxis.FontName = 'Times';
    
    xlabel(' ','FontName','Times','FontSize',40)
    title('@','FontName','Times','FontSize',40)
    
    
    set(h(1),'color','black','linestyle','-','marker','n','linewidth',lsize)
    set(h(2),'color','blue','linestyle','n','marker','o','linewidth',lsize,'MarkerSize',10)
    set(h(3),'color','red','linestyle','n','marker','d','linewidth',lsize,'MarkerSize',10)
    set(h(4),'color',[0,.5,0],'linestyle','--','marker','n','linewidth',2,'MarkerSize',10)
    set(h(5),'color','black','linestyle',':','marker','n','linewidth',lsize)
    xlim([-6+4*rss 6+4*rss]);xticks(-4:1:8);
    ylim([-3 2])