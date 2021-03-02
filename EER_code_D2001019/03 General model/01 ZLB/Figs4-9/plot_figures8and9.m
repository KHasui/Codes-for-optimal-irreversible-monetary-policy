% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.
% More detailed calculation of welfare (Figure 8) is available upon request (Kohei Hasui: khasui@g.matsuyama-u.ac.jp).

close all
clear;
clc
%disp('Analyzing relationship between welfare gain from reversal aversion and the frequency of hitting ZLB under pure discretion.')
disp('NOTE: This code takes at least 30 minutes to complete (If you stop, push "ctrl + c.")')

% 1. Generating pseudo demand shocks:
varg     = (3.72/16)^2;
nshocks  = 13;
[eg,prb] = qnwnorm(nshocks,0,varg); % recall qnwnorm from CompEcon Toolbox
nburn    = 1000; % # of burning periods (default setting is 200)
T        = 10000 + nburn;
npath    = 1000; % # of paths
eg_mat   = zeros(npath,T);
    for k = 1:npath % recall discrand.m from CompEcon Toolbox
        eg_mat(k,:) = (eg(discrand(T,prb)))';
    end
    %save data_shock_process eg_mat T npath


lamdafs = [0,0.0016];
nlf     = size(lamdafs,2);

loss  = zeros(npath,nlf); 
brate = zeros(npath,nlf);
bnum  = zeros(npath,nlf);

% 2. Calculate binding rate and welfare loss:
for ilf  = 1:nlf % loop for lambda_{ir} = 0 and 0.0016.

    for ipath = 1:npath % 1000 paths
    
        g_t      = 0;
        ilag_t   = 0;
        dltlag_t = 0;
    
        y_impmat_irr  = zeros(T,1);
        p_impmat_irr  = zeros(T,1);
        bindrate      = zeros(T,1);
    
        % loading result lamda_ir= 0 and 0.0016:
        eval(['load result_irrev_indlamdaf',num2str(ilf),'_zlb_sig100_kpp24'])
    
        Fi = griddedInterpolant(gmat,ilagmat,dltlagmat,imat,'spline');
        
        % counting binding
        for t = 1:T
            i_t = max(-rss,Fi(g_t,ilag_t,dltlag_t));
            p_impmat_irr(t,:)   = Fp(g_t,ilag_t,dltlag_t);
            y_impmat_irr(t,:)   = Fy(g_t,ilag_t,dltlag_t);
            
            % counting bind of interest rate:
            bindrate(t,:) = max(imin,i_t) == imin;
            
            % state transition
            dltlag_t = sign(i_t-ilag_t);
            ilag_t   = i_t;
            g_t      = rhog*g_t + eg_mat(ipath,t);
        end 
        
        %  discarded initial periods (default setting is 200 periods)
        numbind = sum(bindrate(nburn+1:T,:));
        freqzlb = numbind/(T-nburn)*100;
        nmat    = floor(sqrt(T-nburn)); % finally, reshape the simulated path, n_mat by n_mat matrix
        
        yresh = y_impmat_irr(nburn+1:nburn+nmat^2,:); 
        presh = p_impmat_irr(nburn+1:nburn+nmat^2,:); 
        yresh_mat = reshape(yresh,nmat,nmat);         % ypath: n_mat by n_mat matrix
        presh_mat = reshape(presh,nmat,nmat);         % ppath: n_mat by n_mat matrix
        losse_zlb_mat = zeros(1,nmat); % prepare box of loss
        
        % calculating welfare loss
        for iy=1:nmat
            for ip=1:nmat
                losse_zlb_mat(iy)=losse_zlb_mat(iy)+...
                bet^(ip-1)*((presh_mat(ip,iy))^2+lamdax*(yresh_mat(ip,iy))^2);
            end
        end
    
        L = mean(losse_zlb_mat); % insert average of loss
        loss(ipath,ilf)  = L;
        brate(ipath,ilf) = freqzlb;
        bnum(ipath,ilf)  = numbind;
        
        % Display progress:
        if ipath/npath*100 == 25
            disp(['  25% completed (lambda_ir = ',num2str(lamdafs(ilf)),')...'])
        elseif ipath/npath*100 == 50
            disp(['  50% completed (lambda_ir = ',num2str(lamdafs(ilf)),')...'])
        elseif ipath/npath*100 == 75
            disp(['  75% completed (lambda_ir = ',num2str(lamdafs(ilf)),')...'])
        elseif ipath/npath*100 == 100
            disp(['  100% completed (lambda_ir = ',num2str(lamdafs(ilf)),').'])
            disp('  --')
        end

    end


end

% 3. Sorting bindings and gains

lamdaf = 0.0016;
indlamdaf = find(lamdafs == lamdaf);
inddisc   = find(lamdafs == 0);

bnumdisc = bnum(:,inddisc); % extract # of binding under pure disc
[bnumdisc_sort,ind_bnumdisc] = sort(bnumdisc); % reproduce index of binding

% extract welfare
welfare_pd   = loss(:,inddisc);   % pure discretion
welfare_irr  = loss(:,indlamdaf); % lambda_ir = 0.0016

% sorting welfare
welfare_pd_sort  = welfare_pd(ind_bnumdisc);
welfare_irr_sort = welfare_irr(ind_bnumdisc);

% gains
relgains = (welfare_pd_sort./welfare_irr_sort-1)*100;

nbt = size(bnumdisc_sort,1);

sumrelgains = zeros(nbt,1);
brate_disc  = zeros(nbt,1);
brate_dum   = zeros(nbt,1);
brate_n     = zeros(nbt,1);

indb = 0;
bnum_i_old = 0;

for i = 1:nbt
    
    % extract bindrate
    bnum_i = bnumdisc_sort(i);
    
    % extract relative gains
    relgain_i = relgains(i);
    
    % check bindrate increases or remains
    ind_shift = sign(bnum_i-bnum_i_old);
    indb = indb + ind_shift;
    
    % insert values
    brate_disc(indb) = bnum_i;
    sumrelgains(indb) = sumrelgains(indb) + relgain_i;
    
    brate_dum(indb) = 1;
    brate_n(indb) = brate_n(indb) + 1;
    
    bnum_i_old = bnum_i;
    
end


nb = sum(brate_dum);

brate_disc  = brate_disc(1:nb);
sumrelgains = sumrelgains(1:nb);
brate_n     = brate_n(1:nb);

mean_elgains = sumrelgains./brate_n;

% Plotting figure: Analyzing relationship between welfare gain from reversal aversion and the frequency of hitting ZLB under pure discretion

figure(1)
h = scatter(brate_disc,mean_elgains,65,'MarkerFaceColor',[0.8,0.9,1],'MarkerEdgeColor',[0,0,0.9],'LineWidth',1);
    ylabel('Welfare gain (%)','Fontsize',22,'FontName','Times');
    xlabel(['Number of periods hitting ZLB';'under pure discretion        '],'Fontsize',21,'FontName','Times')
    %title(['$\lambda_{\rm{ir}} = ',num2str(lamdafs(2)),'$'],'Interpreter','latex','Fontsize',27)
    ylim([5 20]);xlim([770 1030])
    
    xticks([800 900 1000])
    
    h_axes = gca;
    h_axes.XAxis.FontSize = 19;
    h_axes.XAxis.FontName = 'Times';
    
    h_axes = gca;
    h_axes.YAxis.FontSize = 19;
    h_axes.YAxis.FontName = 'Times';

    
    
load result_irrev_indlambdaf_zlb_nozlb_sig100_kpp24


left_color  = [0,0,0];right_color = [0,0,0];
    
set(figure(2),'defaultAxesColorOrder',[left_color; right_color]);
    h1 = plot(lamdafs,gain_zlb,lamdafs_nozlb,gain_nozlb,[0.001635 0.001635],[-30 30]);
    
    h_axes = gca;
    h_axes.XAxis.FontSize = 19;
    h_axes.XAxis.FontName = 'Times';
    
    h_axes = gca;
    h_axes.YAxis.FontSize = 19;
    h_axes.YAxis.FontName = 'Times';
    
    set(h1(1),'color','blue','LineWidth',3);
    set(h1(2),'color','red','LineWidth',3,'linestyle',':');
    set(h1(3),'color','black','LineWidth',2,'linestyle',':');
    %set(h1(2),'color','red','marker','o')
    %set(h1(2),'color','red','linestyle','--');
    %xlim([0 100]);
    %ylim([-4 4]);yticks(-3:1:3)
    xlim([0 0.008])
    ylim([-30 30])
    annotation(figure(2),'textbox',...
    [0.281357142857142 0.176190476190476 0.0704285714285719 0.0880952380952399],...
    'String',{'Å©'},...
    'LineStyle','none',...
    'FontSize',24,...
    'FontName','Times',...
    'FitBoxToText','off');
    annotation(figure(2),'textbox',...
    [0.334928571428569 0.166666666666667 0.0704285714285719 0.08809523809524],...
    'String',{'$\lambda_{\rm ir}^*$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',28,...
    'FontName','Times',...
    'FitBoxToText','off');

    ylabel('Welfare gain (%)','FontName','Times','FontSize',24)
    xlabel('$\lambda_{\rm ir}$','Interpreter','Latex','FontSize',28);
    lgd = legend('Welfare gain with ZLB','Welfare gain without ZLB');
    lgd.FontSize = 20;
    lgd.FontName = 'Times';
    lgd.EdgeColor = 'white';
    
    xticks([0 0.002 0.004 0.006 0.008])
    xticklabels([0 0.002 0.004 0.006 0.008])    
    
disp('Fin.')



