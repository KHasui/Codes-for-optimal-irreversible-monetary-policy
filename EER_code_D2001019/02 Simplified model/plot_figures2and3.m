% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.

close all
clear


%
% 1: Impulses
%


% Probabilities: probH denotes Prob(rnt=rL|rn_{t-1}=rH).
%probH = 0.01; % Case (a)
probH = 0.4;  % Case (b) 

probL = 0.875; % q in article

% model params
sig      = 1;
kappa    = 0.024;
rL       = -1/4;
rH       = 1.5/4;
bet      = 1/(1+rH/100);
lambda   = 0.003;
lambdai  = 0.01;


% Horizon and state space
nperiod = 4; % t = 0,1,2,3

imin = 0;
imax = 1.5;
inti = .5;
i_grid = (imin:inti:imax)/4;
d_grid = [-1,0,1];

[imat,dmat] = ndgrid(i_grid,d_grid);

si = imat(:);
sd = dmat(:);

St = [dmat(:),imat(:)];

[~,amat] = ndgrid(i_grid,[0,0,1]);
[~,bmat] = ndgrid(i_grid,[1,0,0]);
[~,cmat] = ndgrid(i_grid,[0,1,1]);

a = amat(:).*St(:,2);
b = bmat(:).*St(:,2) + imax/4*cmat(:);

ns = size(St,1);



%
% Model 1: Policy irreversibility
%

% Policy function (state rH)
ppol_irr_H = zeros(ns,nperiod);
xpol_irr_H = zeros(ns,nperiod);
ipol_irr_H = zeros(ns,nperiod);
val_irr_H  = zeros(ns,nperiod);

% Policy function (state rL)
ppol_irr_L = zeros(ns,nperiod);
xpol_irr_L = zeros(ns,nperiod);
ipol_irr_L = zeros(ns,nperiod);
val_irr_L  = zeros(ns,nperiod);

% --Solve model backward--
% t = 3
t = nperiod;
for i_s = 1:ns
    
    a_irr_is = a(i_s); b_irr_is = b(i_s); i4_vec = a_irr_is:inti/4:b_irr_is;
    
    x_irr_4H_vec = -sig*(i4_vec-rH);
    p_irr_4H_vec = kappa*x_irr_4H_vec;
    [vH,indv_irr_H] = min(0.5*p_irr_4H_vec.^2 + lambda*0.5*x_irr_4H_vec.^2);
    
    ppol_irr_H(i_s,t) = p_irr_4H_vec(indv_irr_H);
    xpol_irr_H(i_s,t) = x_irr_4H_vec(indv_irr_H);
    ipol_irr_H(i_s,t) = i4_vec(indv_irr_H);
    val_irr_H(i_s,t) = vH;
    
    x4L_vec = -sig*(i4_vec-rL);
    p4L_vec = kappa*x4L_vec;
    [vL,indv_irr_L] = min(0.5*p4L_vec.^2 + lambda*0.5*x4L_vec.^2);
    ppol_irr_L(i_s,t) = p4L_vec(indv_irr_L);
    xpol_irr_L(i_s,t) = x4L_vec(indv_irr_L);
    ipol_irr_L(i_s,t) = i4_vec(indv_irr_L);
    val_irr_L(i_s,t) = vL;
    
end
% t = 0,1,2
for t = nperiod-1:-1:1

    for i_s = 1:ns
    
    
        ilag_irr = si(i_s);
        a_irr_is = a(i_s); b_irr_is = b(i_s); i_irr_vec = a_irr_is:inti/4:b_irr_is;
    
        ni_irr_vec = size(i_irr_vec,2);
    
        Ep_irr_H = zeros(1,ni_irr_vec);
        Ex_irr_H = zeros(1,ni_irr_vec);
        Ev_irr_H = zeros(1,ni_irr_vec);
    
        Ep_irr_L = zeros(1,ni_irr_vec);
        Ex_irr_L = zeros(1,ni_irr_vec);
        Ev_irr_L = zeros(1,ni_irr_vec);
    
        for ii = 1:ni_irr_vec
        
            it_irr = i_irr_vec(ii);  dt_irr = sign(it_irr-ilag_irr);
            hit_ind_d_irr = St(:,1)==dt_irr;
            hit_ind_i_irr = St(:,2)==it_irr;
            hit_ind_irr = hit_ind_i_irr+hit_ind_d_irr;
            ind_2_irr = find(hit_ind_irr == 2);
        
            Ep_irr_H(ii) = probH*ppol_irr_L(ind_2_irr,t+1) + (1-probH)*ppol_irr_H(ind_2_irr,t+1);
            Ex_irr_H(ii) = probH*xpol_irr_L(ind_2_irr,t+1) + (1-probH)*xpol_irr_H(ind_2_irr,t+1);
            Ev_irr_H(ii) = probH*val_irr_L(ind_2_irr,t+1)  + (1-probH)*val_irr_H(ind_2_irr,t+1);
        
            Ep_irr_L(ii) = probL*ppol_irr_L(ind_2_irr,t+1) + (1-probL)*ppol_irr_H(ind_2_irr,t+1);
            Ex_irr_L(ii) = probL*xpol_irr_L(ind_2_irr,t+1) + (1-probL)*xpol_irr_H(ind_2_irr,t+1);
            Ev_irr_L(ii) = probL*val_irr_L(ind_2_irr,t+1)  + (1-probL)*val_irr_H(ind_2_irr,t+1);
        
        
        end
    
        xt_irr_H_vec = Ex_irr_H-sig*(i_irr_vec-Ep_irr_H-rH);
        pt_irr_H_vec = bet*Ep_irr_H +kappa*xt_irr_H_vec;
        [vt_irr_H,indv_irr_H] = min(0.5*pt_irr_H_vec.^2 + lambda*0.5*xt_irr_H_vec.^2 + bet*Ev_irr_H);
    
        ppol_irr_H(i_s,t) = pt_irr_H_vec(indv_irr_H);
        xpol_irr_H(i_s,t) = xt_irr_H_vec(indv_irr_H);
        ipol_irr_H(i_s,t) = i_irr_vec(indv_irr_H);
        val_irr_H(i_s,t) = vt_irr_H;
    
        xt_irr_L_vec = Ex_irr_L -sig*(i_irr_vec-Ep_irr_L-rL);
        pt_irr_L_vec = bet*Ep_irr_L + kappa*xt_irr_L_vec;
        [vt_irr_L,indv_irr_L] = min(0.5*pt_irr_L_vec.^2 + lambda*0.5*xt_irr_L_vec.^2 + bet*Ev_irr_L);
    
        ppol_irr_L(i_s,t) = pt_irr_L_vec(indv_irr_L);
        xpol_irr_L(i_s,t) = xt_irr_L_vec(indv_irr_L);
        ipol_irr_L(i_s,t) = i_irr_vec(indv_irr_L);
        val_irr_L(i_s,t) = vt_irr_L;
    
    end

end


% --Simulate path: r^n = rL rL rH rH--
vpath_irr = zeros(nperiod,1);
ppath_irr = zeros(nperiod,1);
xpath_irr = zeros(nperiod,1);
ipath_irr = zeros(nperiod,1);
rpath_irr = zeros(nperiod,1);
dpath_irr = zeros(nperiod,1);
i_0 = 0;
d_0 = 0;

rnpath = rH*ones(1,nperiod);
rnpath(1:1) = rL;
drL    = rnpath == rL;

% t = 0; r^n = rL
    T = 1;

    indd = St(:,1) == d_0;
    indi = St(:,2) == i_0;
    indv1 =find(indd+indi == 2);
    
    v1_irr = val_irr_L(indv1,T);
    i_1 = ipol_irr_L(indv1,T);
    d_1 = sign(i_1-i_0);
    p_1 = ppol_irr_L(indv1,T);
    x_1 = xpol_irr_L(indv1,T);
    Ep_2 = (p_1-kappa*x_1)/bet;
    r_1 = i_1 - Ep_2;
    vpath_irr(T) = v1_irr; 
    ipath_irr(T) = i_1; 
    ppath_irr(T) = p_1; 
    xpath_irr(T) = x_1; 
    dpath_irr(T) = d_1;
    rpath_irr(T) = r_1; 
    
% t = 1,2,3: r^n = rL
for T = 2:nperiod
    
    drLT = drL(T);
    
    i_lag = ipath_irr(T-1);
    d_lag = dpath_irr(T-1);
    indd = St(:,1) == d_lag;
    indi = St(:,2) == i_lag;
    indv_2 =find(indd+indi == 2);
    v_t = drLT*val_irr_L(indv_2,T) + (1-drLT)*val_irr_H(indv_2,T);
    vpath_irr(T) = v_t;
    i_t = drLT*ipol_irr_L(indv_2,T) + (1-drLT)*ipol_irr_H(indv_2,T);
    d_t = sign(i_t-i_lag);
    p_t = drLT*ppol_irr_L(indv_2,T) + (1-drLT)*ppol_irr_H(indv_2,T);
    x_t = drLT*xpol_irr_L(indv_2,T) + (1-drLT)*xpol_irr_H(indv_2,T);
    Ep_t = (p_t-kappa*x_t)/bet;
    r_t = i_t - Ep_t;
    
    ipath_irr(T)  = i_t; 
    ppath_irr(T)  = p_t; 
    xpath_irr(T)  = x_t;
    dpath_irr(T)  = d_t;
    rpath_irr(T)  = r_t;
    
end


%
% Model 2: Interest rate smoothing
%

si = i_grid';

St = si;

a = min(i_grid);
b = max(i_grid);
ns = size(St,1);

ppol_irs_H = zeros(ns,nperiod);
xpol_irs_H = zeros(ns,nperiod);
ipol_irs_H = zeros(ns,nperiod);
val_irs_H  = zeros(ns,nperiod);

ppol_irs_L = zeros(ns,nperiod);
xpol_irs_L = zeros(ns,nperiod);
ipol_irs_L = zeros(ns,nperiod);
val_irs_L  = zeros(ns,nperiod);

% --Solve model backward--
% t = 3
t = nperiod;
for i_s = 1:ns
    
    i4_irs_vec = i_grid;
    i3_irs = i_grid(i_s);
    
    x4_irs_H_vec = -sig*(i4_irs_vec-rH);
    p4_irs_H_vec = kappa*x4_irs_H_vec;
    [v_irs_H,indv_irs_H] = min(0.5*(p4_irs_H_vec.^2 + lambda*x4_irs_H_vec.^2 + lambdai*(i4_irs_vec-i3_irs).^2));
    
    ppol_irs_H(i_s,t) = p4_irs_H_vec(indv_irs_H);
    xpol_irs_H(i_s,t) = x4_irs_H_vec(indv_irs_H);
    ipol_irs_H(i_s,t) = i4_irs_vec(indv_irs_H);
    val_irs_H(i_s,t) = v_irs_H;
    
    x4_irs_L_vec = -sig*(i4_irs_vec-rL);
    p4_irs_L_vec = kappa*x4_irs_L_vec;
    [v_irs_L,indv_irs_L] = min(0.5*(p4_irs_L_vec.^2 + lambda*x4_irs_L_vec.^2 + lambdai*(i4_irs_vec-i3_irs).^2));
    ppol_irs_L(i_s,t) = p4_irs_L_vec(indv_irs_L);
    xpol_irs_L(i_s,t) = x4_irs_L_vec(indv_irs_L);
    ipol_irs_L(i_s,t) = i4_irs_vec(indv_irs_L);
    val_irs_L(i_s,t) = v_irs_L;
    
end
% t = 0,1,2
for t = nperiod-1:-1:1

    for i_s = 1:ns
    
    
        ilag_irs_ = si(i_s);
        i_irs__vec = i_grid;
    
        ni_irs_vec = size(i_irs__vec,2);
    
        Ep_irs_H = zeros(1,ni_irs_vec);
        Ex_irs_H = zeros(1,ni_irs_vec);
        Ev_irs_H = zeros(1,ni_irs_vec);
        Ep_irs_L = zeros(1,ni_irs_vec);
        Ex_irs_L = zeros(1,ni_irs_vec);
        Ev_irs_L = zeros(1,ni_irs_vec);
    
        for ii = 1:ni_irs_vec
        
            it = i_irs__vec(ii);
            ind_it = find(St == it);
            Ep_irs_H(ii) = probH*ppol_irs_L(ind_it,t+1) + (1-probH)*ppol_irs_H(ind_it,t+1);
            Ex_irs_H(ii) = probH*xpol_irs_L(ind_it,t+1) + (1-probH)*xpol_irs_H(ind_it,t+1);
            Ev_irs_H(ii) = probH*val_irs_L(ind_it,t+1)  + (1-probH)*val_irs_H(ind_it,t+1);
            
            Ep_irs_L(ii) = probL*ppol_irs_L(ind_it,t+1) + (1-probL)*ppol_irs_H(ind_it,t+1);
            Ex_irs_L(ii) = probL*xpol_irs_L(ind_it,t+1) + (1-probL)*xpol_irs_H(ind_it,t+1);
            Ev_irs_L(ii) = probL*val_irs_L(ind_it,t+1)  + (1-probL)*val_irs_H(ind_it,t+1);
        
        end
        
        xt_irs_H_vec = Ex_irs_H-sig*(i_irs__vec-Ep_irs_H-rH);
        pt_irs_H_vec = bet*Ep_irs_H +kappa*xt_irs_H_vec;
        [vt_irs_H,indv_irs_H] = min(0.5*pt_irs_H_vec.^2 + lambda*0.5*xt_irs_H_vec.^2 + lambdai*0.5*(i_irs__vec-ilag_irs_).^2 + bet*Ev_irs_H);
        
        ppol_irs_H(i_s,t) = pt_irs_H_vec(indv_irs_H);
        xpol_irs_H(i_s,t) = xt_irs_H_vec(indv_irs_H);
        ipol_irs_H(i_s,t) = i_irs__vec(indv_irs_H);
        val_irs_H(i_s,t) = vt_irs_H;
        
        xt_irs_L_vec = Ex_irs_L -sig*(i_irs__vec-Ep_irs_L-rL);
        pt_irs_L_vec = bet*Ep_irs_L + kappa*xt_irs_L_vec;
        [vt_irs_L,indv_irs_L] = min(0.5*pt_irs_L_vec.^2 + lambda*0.5*xt_irs_L_vec.^2 + lambdai*0.5*(i_irs__vec-ilag_irs_).^2 + bet*Ev_irs_L);
        ppol_irs_L(i_s,t) = pt_irs_L_vec(indv_irs_L);
        xpol_irs_L(i_s,t) = xt_irs_L_vec(indv_irs_L);
        ipol_irs_L(i_s,t) = i_irs__vec(indv_irs_L);
        val_irs_L(i_s,t) = vt_irs_L;
    
    end

end

% --Simulate path: r^n = rL rL rH rH--
vpath_irs = zeros(nperiod,1);
ppath_irs = zeros(nperiod,1);
xpath_irs = zeros(nperiod,1);
ipath_irs = zeros(nperiod,1);
rpath_irs = zeros(nperiod,1);
dpath_irs = zeros(nperiod,1);
i_0 = 0;

rnpath = rH*ones(1,nperiod);
rnpath(1:1) = rL;
drL    = rnpath == rL;

% t = 1; r^n = rL
    T = 1;
    
    indi = St == i_0;
    indv1 =find(indi == 1);

    v1_irs = val_irs_L(indv1,T);
    i_1    = ipol_irs_L(indv1,T);
    d_1    = sign(i_1-i_0);
    p_1    = ppol_irs_L(indv1,T);
    x_1    = xpol_irs_L(indv1,T);
    Ep_2   = (p_1-kappa*x_1)/bet;
    r_1    = i_1 - Ep_2;
    
    vpath_irs(T) = v1_irs; 
    ipath_irs(T) = i_1; 
    ppath_irs(T) = p_1; 
    xpath_irs(T) = x_1; 
    dpath_irs(T) = d_1;
    rpath_irs(T) = r_1;

% t = 1:3 r^n = rL
for T = 2:nperiod
    
    drLT = drL(T);
    
    i_lag = ipath_irs(T-1);
    d_lag = dpath_irs(T-1);
    
    indi   = St == i_lag;
    indv_2 = find(indi == 1);
    
    v_t = drLT*val_irs_L(indv_2,T) + (1-drLT)*val_irs_H(indv_2,T);
    vpath_irs(T) = v_t;
    i_t  = drLT*ipol_irs_L(indv_2,T) + (1-drLT)*ipol_irs_H(indv_2,T);
    d_t  = sign(i_t-i_lag);
    p_t  = drLT*ppol_irs_L(indv_2,T) + (1-drLT)*ppol_irs_H(indv_2,T);
    x_t  = drLT*xpol_irs_L(indv_2,T) + (1-drLT)*xpol_irs_H(indv_2,T);
    Ep_t = (p_t-kappa*x_t)/bet;
    r_t = i_t - Ep_t;
    
    ipath_irs(T)  = i_t; 
    ppath_irs(T)  = p_t; 
    xpath_irs(T)  = x_t;
    dpath_irs(T)  = d_t;
    rpath_irs(T)  = r_t;
    
end



%
% Model 3: Pure discretion
%

ppol_pd_H = zeros(1,nperiod);
xpol_pd_H = zeros(1,nperiod);
ipol_pd_H = zeros(1,nperiod);
val_pd_H  = zeros(1,nperiod);

ppol_pd_L = zeros(1,nperiod);
xpol_pd_L = zeros(1,nperiod);
ipol_pd_L = zeros(1,nperiod);
val_pd_L  = zeros(1,nperiod);

ivec = (0:inti:1.5)/4;

% --Solve model backward--
% t = 3
    t = nperiod;
    
    x_pd_4H_vec = -sig*(ivec-rH);
    p_pd_4H_vec = kappa*x_pd_4H_vec;
    [v_pd_H,indv_pd_H] = min(0.5*p_pd_4H_vec.^2 + lambda*0.5*x_pd_4H_vec.^2);
    
    ppol_pd_H(:,t) = p_pd_4H_vec(indv_pd_H);
    xpol_pd_H(:,t) = x_pd_4H_vec(indv_pd_H);
    ipol_pd_H(:,t) = ivec(indv_pd_H);
    val_pd_H(:,t) = v_pd_H;
    
    x4_pd_L_vec = -sig*(ivec-rL);
    p4_pd_L_vec = kappa*x4_pd_L_vec;
    [v_pd_L,indv_pd_L] = min(0.5*p4_pd_L_vec.^2 + lambda*0.5*x4_pd_L_vec.^2);
    ppol_pd_L(:,t) = p4_pd_L_vec(indv_pd_L);
    xpol_pd_L(:,t) = x4_pd_L_vec(indv_pd_L);
    ipol_pd_L(:,t) = ivec(indv_pd_L);
    val_pd_L(:,t) = v_pd_L;
    
% t = 0,1,2
for t = nperiod-1:-1:1
%for i_s = 1:ns
    
    Ep_pd_H = probH*ppol_pd_L(t+1) + (1-probH)*ppol_pd_H(t+1);
    Ex_pd_H = probH*xpol_pd_L(t+1) + (1-probH)*xpol_pd_H(t+1);
    Ev_pd_H = probH*val_pd_L(t+1)  + (1-probH)*val_pd_H(t+1);
    
    Ep_pd_L = probL*ppol_pd_L(t+1) + (1-probL)*ppol_pd_H(t+1);
    Ex_pd_L = probL*xpol_pd_L(t+1) + (1-probL)*xpol_pd_H(t+1);
    Ev_pd_L = probL*val_pd_L(t+1)  + (1-probL)*val_pd_H(t+1);
    
    xt_pd_H_vec = Ex_pd_H-sig*(ivec-Ep_pd_H-rH);
    pt_pd_H_vec = bet*Ep_pd_H +kappa*xt_pd_H_vec;
    [vtH,indv_pd_H] = min(0.5*pt_pd_H_vec.^2 + lambda*0.5*xt_pd_H_vec.^2 + bet*Ev_pd_H);
    
    ppol_pd_H(t) = pt_pd_H_vec(indv_pd_H);
    xpol_pd_H(t) = xt_pd_H_vec(indv_pd_H);
    ipol_pd_H(t) = ivec(indv_pd_H);
    val_pd_H(t) = vtH;
    
    xt_pd_L_vec = Ex_pd_L -sig*(ivec-Ep_pd_L-rL);
    pt_pd_L_vec = bet*Ep_pd_L + kappa*xt_pd_L_vec;
    [vt_pd_L,indv_pd_L] = min(0.5*pt_pd_L_vec.^2 + lambda*0.5*xt_pd_L_vec.^2 + bet*Ev_pd_L);
    ppol_pd_L(t) = pt_pd_L_vec(indv_pd_L);
    xpol_pd_L(t) = xt_pd_L_vec(indv_pd_L);
    ipol_pd_L(t) = ivec(indv_pd_L);
    val_pd_L(t) = vt_pd_L;
    
end

% --Simulate path: r^n = rL rL rH rH--
vpath_pd = zeros(nperiod,1);
ppath_pd = zeros(nperiod,1);
xpath_pd = zeros(nperiod,1);
ipath_pd = zeros(nperiod,1);
rnpath_pd = zeros(nperiod,1);
rpath_pd = zeros(nperiod,1);
i_0 = rH;

% t = 1; r^n = rL
    T = 1;
    v1_pd = val_pd_L(T);
    vpath_pd(T) = v1_pd; 
    i_1 = ipol_pd_L(T);
    p_1 = ppol_pd_L(T);
    x_1 = xpol_pd_L(T);
    Ep_2 = (p_1-kappa*x_1)/bet;
    r_1 = i_1 - Ep_2; 
    ppath_pd(T) = p_1; 
    xpath_pd(T) = x_1;
    rpath_pd(T) = r_1;
    ipath_pd(T) = i_1;
    %rnpath_pd(T) = rL;

% t = 2: r^n = rL
for T = 2:nperiod
    
    drLT = drL(T);
    
    v_t_pd = drLT*val_pd_L(T)+(1-drLT)*val_pd_H(T);
    vpath_pd(T) = v_t_pd;
    i_t = drLT*ipol_pd_L(T)+(1-drLT)*ipol_pd_H(T);
    p_t = drLT*ppol_pd_L(T)+(1-drLT)*ppol_pd_H(T);
    x_t = drLT*xpol_pd_L(T)+(1-drLT)*xpol_pd_H(T);
    Ep_t = (p_t-kappa*x_t)/bet;
    r_t = i_t - Ep_t; 
    ipath_pd(T) = i_t; 
    ppath_pd(T) = p_t; 
    xpath_pd(T) = x_t;
    rpath_pd(T) = r_t;

end


%
% Model 4: Commitment
%

rngrid   = [rL rH];
rn       = rngrid';
rnpath   = [rL rH rH rH];
i_choice = [0 0.5 1 1.5]/4;

[i0mat,i1mat,i2mat,i3mat] = ndgrid(i_choice,i_choice,i_choice,i_choice);

imat = [i0mat(:),i1mat(:),i2mat(:),i3mat(:)];
[ni,~] = size(imat);

nrn = size(rngrid,2);
indL = 1;
indH = 2;

vpath_com3 = zeros(ni,nperiod+1,nrn);
ppath_com3 = zeros(ni,nperiod+1,nrn);
xpath_com3 = zeros(ni,nperiod+1,nrn);

vimps_com  = zeros(ni,nperiod);
pimps_com  = zeros(ni,nperiod);
ximps_com  = zeros(ni,nperiod);
epimps_com = zeros(ni,nperiod);

% --Solve model backwardly--
for ii = 1:ni
    
    ipath_i = imat(ii,:);
    
    for t = nperiod:-1:1
        
        rn_t = rnpath(t);
        i_t  = ipath_i(t);
        dum_rH = rn_t == rH; 
        
        Ep = [probL*(ppath_com3(ii,t+1,indL)) + (1-probL)*(ppath_com3(ii,t+1,indH));probH*(ppath_com3(ii,t+1,indL)) + (1-probH)*(ppath_com3(ii,t+1,indH))];
        Ex = [probL*(xpath_com3(ii,t+1,indL)) + (1-probL)*(xpath_com3(ii,t+1,indH));probH*(xpath_com3(ii,t+1,indL)) + (1-probH)*(xpath_com3(ii,t+1,indH))];
        Ev = [probL*(vpath_com3(ii,t+1,indL)) + (1-probL)*(vpath_com3(ii,t+1,indH));probH*(vpath_com3(ii,t+1,indL)) + (1-probH)*(vpath_com3(ii,t+1,indH))];
        
        x_t = Ex - sig*(i_t - Ep - rn);
        p_t = bet*Ep + kappa*x_t;
        v_t = 0.5*(p_t.^2 + lambda*x_t.^2) + bet*Ev;
        
        ppath_com3(ii,t,indL) = p_t(indL); ppath_com3(ii,t,indH) = p_t(indH);
        xpath_com3(ii,t,indL) = x_t(indL); xpath_com3(ii,t,indH) = x_t(indH);
        vpath_com3(ii,t,indL) = v_t(indL); vpath_com3(ii,t,indH) = v_t(indH);
        
        vimps_com(ii,t) = (1-dum_rH)*v_t(indL) + dum_rH*v_t(indH);
        pimps_com(ii,t) = (1-dum_rH)*p_t(indL) + dum_rH*p_t(indH);
        ximps_com(ii,t) = (1-dum_rH)*x_t(indL) + dum_rH*x_t(indH);
        epimps_com(ii,t) = (1-dum_rH)*Ep(indL) + dum_rH*Ep(indH);
    
    end
    
    
end


[vmin,i_optpath] = min(vimps_com(:,1));
indv_double = sum(vimps_com(:,1) == vmin);
v1_com = vmin;

if indv_double ~= 1
    warning('Minimum value is not unique.')
end

xpath_com = ximps_com(i_optpath,:);
ppath_com = pimps_com(i_optpath,:);
ipath_com = imat(i_optpath,:);
rpath_com = ipath_com-epimps_com(i_optpath,:);


%
% Plot figure
%

tgrid = 0:3;
titlefontsize = 17;
    
figure(1)
subplot(1,2,1)
    h = plot(tgrid,4*ipath_com,tgrid,4*ipath_pd,tgrid,4*ipath_irs,tgrid,4*ipath_irr,tgrid,1.5*ones(size(tgrid)));
    legend('Commitment','Pure discretion','Interest rate smoothing','With irreversibility','Interpreter','Latex','fontsize',10,'EdgeColor','n')
    set(h(1),'color','black','linewidth',1);
    set(h(2),'color','black','linewidth',2,'linestyle',':','Marker','^');
    set(h(3),'color','red','linewidth',1.5,'linestyle','--','Marker','x');
    set(h(4),'color','blue','linewidth',1,'Marker','o');
    set(h(5),'color','black','linestyle','--');
    xticks(tgrid)
    ylim([-0.5 2])
    title(['$p = ',num2str(probH),'$'],'Interpreter','Latex','fontsize',titlefontsize);
    ylabel('Nominal interest rate','fontsize',titlefontsize,'FontName','Times')
    xlabel('Periods','fontsize',titlefontsize,'FontName','Times');
    
    h_axes = gca;
    h_axes.XAxis.FontSize = titlefontsize;
    h_axes.XAxis.FontName = 'Times';
    
    h_axes = gca;
    h_axes.YAxis.FontSize = titlefontsize;
    h_axes.YAxis.FontName = 'Times';
    
subplot(1,2,2)
    h = plot(tgrid,4*rnpath,tgrid,4*rpath_com,tgrid,4*rpath_pd,tgrid,4*rpath_irs,tgrid,4*rpath_irr,tgrid,1.5*ones(size(tgrid)));
    set(h(1),'color',[0.6,0.6,0.6],'linewidth',1);
    set(h(2),'color','black','linewidth',1);
    set(h(3),'color','black','linewidth',2,'linestyle',':','Marker','^');
    set(h(4),'color','red','linewidth',1.5,'linestyle','--','Marker','x');
    set(h(5),'color','blue','linewidth',1,'Marker','o');
    set(h(6),'color','black','linestyle','--');
    xticks(tgrid)
    ylim([-2 2])
    title(['$p = ',num2str(probH),'$'],'Interpreter','Latex','fontsize',titlefontsize);
    ylabel('Real interest rate','fontsize',titlefontsize,'FontName','Times');
    xlabel('Periods','fontsize',titlefontsize,'FontName','Times');
    legend(h(1),'Natural rate','Interpreter','Latex','fontsize',10,'EdgeColor','n')
    
    h_axes = gca;
    h_axes.XAxis.FontSize = titlefontsize;
    h_axes.XAxis.FontName = 'Times';
    
    h_axes = gca;
    h_axes.YAxis.FontSize = titlefontsize;
    h_axes.YAxis.FontName = 'Times';
    
    
%
% 2: i^* and probablities
%

params = [sig,kappa,bet,rL,rH,lambda];

pHgrid = (0:0.001:1)';
npH = size(pHgrid,1);
i0grid = zeros(npH,1);
i1grid = zeros(npH,1);
i2grid = zeros(npH,1);
v1grid = zeros(npH,1);

v1i0grid = zeros(npH,1);
v1in0grid = zeros(npH,1);

v105grid = zeros(npH,1);
v11grid = zeros(npH,1);
v115grid = zeros(npH,1);

indv1grid = zeros(npH,1);
indv2grid = zeros(npH,1);

for j = 1:npH
    
    probHj    = pHgrid(j);
    
    % recall function
    [i0j,i1j,i2j,v1j,v1i0j,v1in0j,v205j,v21j,v215j,indv1j,indv2j] = func_i1pH(probHj,probL,params);
    
    i0grid(j) = 4*i0j;
    i1grid(j) = 4*i1j;
    i2grid(j) = 4*i2j;
    v1grid(j) = v1j;
    
    v1i0grid(j)  = v1i0j;
    v1in0grid(j) = v1in0j;
    
    v105grid(j) = v205j;
    v11grid(j)  = v21j;
    v115grid(j) = v215j;
    
    indv1grid(j) = indv1j;
    indv2grid(j) = indv2j;
     
end

ind_n0 = find(i0grid ~= 0);

N1 = v1i0grid-v1in0grid;
N1(ind_n0) = NaN;
[~,ind_probHstar] = min(abs(N1));
pHstar = pHgrid(ind_probHstar);

N105 = v1i0grid-v105grid;

%
% Plot figure
%
figure(2)
    h = plot(pHgrid,i1grid,[pHstar pHstar],[-10 10]);
        xlabel('$p_{\rm H}$','FontSize',titlefontsize,'Interpreter','Latex')
        set(h(1),'LineWidth',1.8,'color','black')
        set(h(2),'LineWidth',0.5,'color','black','LineStyle',':')
        ylabel('$4\times i_1^*$','FontSize',titlefontsize,'Interpreter','Latex')
        ylim([0 1.5])
    
    h_axes = gca;
    h_axes.XAxis.FontSize = titlefontsize;
    h_axes.XAxis.FontName = 'Times';
    
    h_axes = gca;
    h_axes.YAxis.FontSize = titlefontsize;
    h_axes.YAxis.FontName = 'Times';