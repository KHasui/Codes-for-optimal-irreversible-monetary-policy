% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.

close all
clear
clc

disp('Main file: model with the policy irreversibility (without the ZLB).')
disp(' NOTE: This code takes at least 30 hours to complete. If you stop the process, push Keys "Ctrl + c."')

%
% 1. Parameters
%

% Model params
rss    = 1.5/4;
bet    = 1/(1+rss/100);
sig    = 1;
alfa   = .66;
tht    = 7.66;
omg    = .47;
kpp    = 0.024;
lamdax = 0.003;
lamdaf = 10; % degree of penalty 
lamdai = 0.003;
varg   = (0.25*3.72/4)^2;
rhog   = 0.6;

% grid of lambda_f
lamdafs = [0,0.0016,10];

%for lamdaf = lamdafs

    indlamdaf = find(lamdafs == lamdaf);
    strict_irreversi = lamdaf == 10;

    % pack params
    params = [rss; bet; sig; alfa; tht; omg; kpp; lamdax; lamdaf; lamdai; varg; rhog];

    % params for computation
    tol   = 1.48e-006;
    maxit = 4000;
    
%
% 2. State space
%

    % define state space
    nshocks = 13;
    [eg, prb]= qnwnorm(nshocks,0,varg); % Gaussian Quadrature (CompEcon Toolbox)
    
    % make grids of g
    gmax = 4*sqrt(0.7^2/(1-rhog^2));
    gmin = -gmax;
    ng   = 30;
    sg   = linspace(gmin,gmax,ng)';  % grid of g

    % make grids of lagged interest rate (NOTE: grids are slightly changed for the case of lambda_ir = 10)
    imax  = 8*(1-strict_irreversi) + 11*strict_irreversi;
    imin  = -imax;
    inti  = 0.5/4*(1-strict_irreversi) + 0.5/2.7*strict_irreversi;
    ni    = round((imax-imin)/inti+1);
    silag = linspace(imin,imax,ni)'; % grids of ilag

    % sign of shint in interes rate
    ndlt = 3;
    sdlt  = [-1;0;1];

    % state space
    [gmat,ilagmat,dltlagmat] = ndgrid(sg,silag,sdlt);

    g       = gmat(:);
    ilag    = ilagmat(:);
    dltlag  = dltlagmat(:);

    % the size of state space
    ns = ng*ni*ndlt;
    nd = [ng,ni,ndlt];

    expands = ones(1,nshocks);
    gnext = rhog*g*expands + ones(ns,1)*eg';
    %disp(['gnextmax = ',num2str(max(max(gnext)))]) 
    %disp(['gmax = ',num2str(gmax)])

    % check extrapolation for g
    if gmax - max(max(gnext)) < 0
        warning('Extrapolation is detected in upper bound of g!')
        disp(['gmax = ',num2str(gmax),'  upper of gnext = ',num2str( max(max(gnext)))])
    end
    if gmin - min(min(gnext)) > 0
        warning('Extrapolation is detected in lower bound of g!')
        disp(['gmin = ',num2str(gmin),'  lower of gnext = ',num2str( min(min(gnext)))])
    end

%
% 3. Solve the model
%

    vmat = 0*ones(nd);
    pmat = 0*ones(nd);
    ymat = 0*ones(nd);

    v = 0*ones(ns,1);
    p = 0*ones(ns,1);
    y = 0*ones(ns,1);
    i = 0*ones(ns,1);

    v_insert = 0*ones(ns,1);
    p_insert = 0*ones(ns,1);
    y_insert = 0*ones(ns,1);
    i_insert = 0*ones(ns,1);
    
    x = 0;
    options = optimset('Display','off');

    % search area for nominal interest rate
    b = imax;
    a = imin;


    for it = 1:maxit
    
    vold = v;
    
    % function approximation
    Fv = griddedInterpolant(gmat,ilagmat,dltlagmat,vmat,'spline');
    Fp = griddedInterpolant(gmat,ilagmat,dltlagmat,pmat,'spline');
    Fy = griddedInterpolant(gmat,ilagmat,dltlagmat,ymat,'spline');
    
        for i_s = 1:ns
            g_is     = g(i_s);
            gnext_is = gnext(i_s,:);
            ilag_is  = ilag(i_s,:);
            dltlag_is   = dltlag(i_s,:);
        
            % recall fminbnd
            [i_is,v_is] = fminbnd(@objective_irrev,a,b,options,g_is,ilag_is,dltlag_is,gnext_is,Fv,Fy,Fp,params,expands,prb);
            [y_is,p_is] = nkmodel_irrev(i_is,g_is,ilag_is,gnext_is,Fy,Fp,params,expands,prb);
        
            v_insert(i_s,:) = v_is;
            p_insert(i_s,:) = p_is;
            y_insert(i_s,:) = y_is;
            i_insert(i_s,:) = i_is;
        
        end
    
        v = v_insert;
    
        vmat = reshape(v_insert,nd);
        pmat = reshape(p_insert,nd);
        ymat = reshape(y_insert,nd);
        imat = reshape(i_insert,nd);
    
        % check convergence
        residv = abs((v-vold)./vold);
        maxresidv = max(residv);
        if maxresidv < tol
            break
        end
    
    end%for it = 1:maxit
    


eval(['save result_irrev_indlamdaf',num2str(indlamdaf),'_nozlb_sig',num2str(round(sig*100)),'_kpp',num2str(round(kpp*1000))])

