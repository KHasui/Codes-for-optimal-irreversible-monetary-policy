function [i_1,i_2,i_3,v2,v2i0,v2in0,v205,v21,v215,indv1,indv2] = func_i1pH(probH,probL,params)
% Matlab function for deriving optimal path of interest rate and value under the model with policy irreversibility
% INPUT
%  probH: Prob(rnt=rL|rn_{t-1}=rH)
%  probL: Prob(rnt=rH|rn_{t-1}=rL)
%  params: Vector of model parameters
%
% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.

sig    = params(1);
kappa  = params(2);
bet    = params(3);
rL     = params(4);
rH     = params(5);
lambda = params(6);


imin = 0;
imax = 1.5;
inti = .5;
i_grid = (imin:inti:imax)/4;
d_grid = [-1,0,1];

[imat,dmat] = ndgrid(i_grid,d_grid);
si = imat(:);
St = [dmat(:),imat(:)];

[~,amat] = ndgrid(i_grid,[0,0,1]);
[~,bmat] = ndgrid(i_grid,[1,0,0]);
[~,cmat] = ndgrid(i_grid,[0,1,1]);

a = amat(:).*St(:,2);
b = bmat(:).*St(:,2) + imax/4*cmat(:);

ns = size(St,1);
nivec = size(i_grid,2);

ppol_H = zeros(ns,4);
xpol_H = zeros(ns,4);
ipol_H = zeros(ns,4);
val_H  = zeros(ns,4);

ppol_L = zeros(ns,4);
xpol_L = zeros(ns,4);
ipol_L = zeros(ns,4);
val_L  = zeros(ns,4);


% final period (t = 3 in the paper)
t = 4;
for i_s = 1:ns
    
    
    a_is = a(i_s); b_is = b(i_s); i4_vec = a_is:inti/4:b_is;
    
    x4H_vec = -sig*(i4_vec-rH);
    p4H_vec = kappa*x4H_vec;
    [vH,indvH] = min(0.5*p4H_vec.^2 + lambda*0.5*x4H_vec.^2);
    
    ppol_H(i_s,t) = p4H_vec(indvH);
    xpol_H(i_s,t) = x4H_vec(indvH);
    ipol_H(i_s,t) = i4_vec(indvH);
    val_H(i_s,t) = vH;
    
    x4L_vec = -sig*(i4_vec-rL);
    p4L_vec = kappa*x4L_vec;
    [vL,indvL] = min(0.5*p4L_vec.^2 + lambda*0.5*x4L_vec.^2);
    ppol_L(i_s,t) = p4L_vec(indvL);
    xpol_L(i_s,t) = x4L_vec(indvL);
    ipol_L(i_s,t) = i4_vec(indvL);
    val_L(i_s,t) = vL;
    
end


% t = 2 in the paper
t = 3;
for i_s = 1:ns
    
    
    ilag = si(i_s);
    a_is = a(i_s); b_is = b(i_s); i3_vec = a_is:inti/4:b_is;
    
    ni3 = size(i3_vec,2);
    
    Ep4H = zeros(1,ni3);
    Ex4H = zeros(1,ni3);
    Ev4H = zeros(1,ni3);
    
    Ep4L = zeros(1,ni3);
    Ex4L = zeros(1,ni3);
    Ev4L = zeros(1,ni3);
    
    for ii = 1:ni3
        
        it = i3_vec(ii);  dt = sign(it-ilag);
        hit_ind_d = St(:,1)==dt;
        hit_ind_i = St(:,2)==it;
        hit_ind = hit_ind_i+hit_ind_d;
        ind4 = find(hit_ind == 2);
        
        Ep4H(ii) = probH*ppol_L(ind4,t+1) + (1-probH)*ppol_H(ind4,t+1);
        Ex4H(ii) = probH*xpol_L(ind4,t+1) + (1-probH)*xpol_H(ind4,t+1);
        Ev4H(ii) = probH*val_L(ind4,t+1)  + (1-probH)*val_H(ind4,t+1);
        
        Ep4L(ii) = probL*ppol_L(ind4,t+1) + (1-probL)*ppol_H(ind4,t+1);
        Ex4L(ii) = probL*xpol_L(ind4,t+1) + (1-probL)*xpol_H(ind4,t+1);
        Ev4L(ii) = probL*val_L(ind4,t+1)  + (1-probL)*val_H(ind4,t+1);
        
        
    end
    
    x3H_vec = Ex4H-sig*(i3_vec-Ep4H-rH);
    p3H_vec = bet*Ep4H +kappa*x3H_vec;
    [v3H,indvH] = min(0.5*p3H_vec.^2 + lambda*0.5*x3H_vec.^2 + bet*Ev4H);
    
    ppol_H(i_s,t) = p3H_vec(indvH);
    xpol_H(i_s,t) = x3H_vec(indvH);
    ipol_H(i_s,t) = i3_vec(indvH);
    val_H(i_s,t) = v3H;
    
    x3L_vec = Ex4L -sig*(i3_vec-Ep4L-rL);
    p3L_vec = bet*Ep4L + kappa*x3L_vec;
    [v3L,indvL] = min(0.5*p3L_vec.^2 + lambda*0.5*x3L_vec.^2 + bet*Ev4L);
    
    ppol_L(i_s,t) = p3L_vec(indvL);
    xpol_L(i_s,t) = x3L_vec(indvL);
    ipol_L(i_s,t) = i3_vec(indvL);
    val_L(i_s,t) = v3L;
    
end


% t = 1 in the paper
t = 2;
val2H0 = zeros(ns,nivec);
for i_s = 1:ns
    
    
    ilag = si(i_s);
    a_is = a(i_s);
    b_is = b(i_s);
    i2_vec = a_is:inti/4:b_is;
    
    ni2 = size(i2_vec,2);
    
    Ep3H = zeros(1,ni2);
    Ex3H = zeros(1,ni2);
    Ev3H = zeros(1,ni2);
    
    Ep3L = zeros(1,ni2);
    Ex3L = zeros(1,ni2);
    Ev3L = zeros(1,ni2);
    
    Ep3H0 = zeros(1,ni2);
    Ex3H0 = zeros(1,ni2);
    Ev3H0 = zeros(1,ni2);
    
    Ep3L0 = zeros(1,ni2);
    Ex3L0 = zeros(1,ni2);
    Ev3L0 = zeros(1,ni2);
    
    for ii = 1:ni2
        it = i2_vec(ii);  dt = sign(it-ilag);
        hit_ind_d = St(:,1)==dt;
        hit_ind_i = St(:,2)==it;
        hit_ind = hit_ind_i+hit_ind_d;
        ind3 = find(hit_ind == 2);
        
        Ep3H(ii) = probH*ppol_L(ind3,t+1) + (1-probH)*ppol_H(ind3,t+1);
        Ex3H(ii) = probH*xpol_L(ind3,t+1) + (1-probH)*xpol_H(ind3,t+1);
        Ev3H(ii) = probH*val_L(ind3,t+1) + (1-probH)*val_H(ind3,t+1);
        
        Ep3L(ii) = probL*ppol_L(ind3,t+1) + (1-probL)*ppol_H(ind3,t+1);
        Ex3L(ii) = probL*xpol_L(ind3,t+1) + (1-probL)*xpol_H(ind3,t+1);
        Ev3L(ii) = probL*val_L(ind3,t+1) + (1-probL)*val_H(ind3,t+1);
    end
    
    
    for ii = 1:nivec
    %for it = 0
        it = i_grid(ii);  dt = sign(it-ilag);
        hit_ind_d = St(:,1)==dt;
        hit_ind_i = St(:,2)==it;
        hit_ind = hit_ind_i+hit_ind_d;
        ind3 = find(hit_ind == 2);
        
        Ep3H0(ii) = probH*ppol_L(ind3,t+1) + (1-probH)*ppol_H(ind3,t+1);
        Ex3H0(ii) = probH*xpol_L(ind3,t+1) + (1-probH)*xpol_H(ind3,t+1);
        Ev3H0(ii) = probH*val_L(ind3,t+1) + (1-probH)*val_H(ind3,t+1);
        
        Ep3L0(ii) = probL*ppol_L(ind3,t+1) + (1-probL)*ppol_H(ind3,t+1);
        Ex3L0(ii) = probL*xpol_L(ind3,t+1) + (1-probL)*xpol_H(ind3,t+1);
        Ev3L0(ii) = probL*val_L(ind3,t+1) + (1-probL)*val_H(ind3,t+1);
    end
    
    x2H_vec = Ex3H-sig*(i2_vec-Ep3H-rH);
    p2H_vec = bet*Ep3H +kappa*x2H_vec;
    [v2H,indvH] = min(0.5*p2H_vec.^2 + lambda*0.5*x2H_vec.^2 + bet*Ev3H);    
    
    ppol_H(i_s,t) = p2H_vec(indvH);
    xpol_H(i_s,t) = x2H_vec(indvH);
    ipol_H(i_s,t) = i2_vec(indvH);
    val_H(i_s,t) = v2H;
    
    x2L_vec = Ex3L -sig*(i2_vec-Ep3L-rL);
    p2L_vec = bet*Ep3L + kappa*x2L_vec;
    [v2L,indvL] = min(0.5*p2L_vec.^2 + lambda*0.5*x2L_vec.^2 + bet*Ev3L);
    ppol_L(i_s,t) = p2L_vec(indvL);
    xpol_L(i_s,t) = x2L_vec(indvL);
    ipol_L(i_s,t) = i2_vec(indvL);
    val_L(i_s,t) = v2L;
    
    
    x2H0 = Ex3H0-sig*(i_grid-Ep3H0-rH);
    p2H0 = bet*Ep3H0 +kappa*x2H0;
    v2H0 = 0.5*p2H0.^2 + lambda*0.5*x2H0.^2 + bet*Ev3H0;
    val2H0(i_s,:) = v2H0;
    
end


% intial period: (t = 0 in the paper)
t = 1;
for i_s = 1:ns
    
    
    ilag = si(i_s);
    a_is = a(i_s); b_is = b(i_s);
    i1_vec = a_is:inti/4:b_is;
    
    ni1 = size(i1_vec,2);
    
    Ep2H = zeros(1,ni1);
    Ex2H = zeros(1,ni1);
    Ev2H = zeros(1,ni1);
    
    Ep2L = zeros(1,ni1);
    Ex2L = zeros(1,ni1);
    Ev2L = zeros(1,ni1);
    
    for ii = 1:ni1
        
        it = i1_vec(ii);  dt = sign(it-ilag);
        hit_ind_d = St(:,1)==dt;
        hit_ind_i = St(:,2)==it;
        hit_ind = hit_ind_i+hit_ind_d;
        ind2 = find(hit_ind == 2);
        
        Ep2H(ii) = probH*ppol_L(ind2,t+1) + (1-probH)*ppol_H(ind2,t+1);
        Ex2H(ii) = probH*xpol_L(ind2,t+1) + (1-probH)*xpol_H(ind2,t+1);
        Ev2H(ii) = probH*val_L(ind2,t+1) + (1-probH)*val_H(ind2,t+1);
        
        Ep2L(ii) = probL*ppol_L(ind2,t+1) + (1-probL)*ppol_H(ind2,t+1);
        Ex2L(ii) = probL*xpol_L(ind2,t+1) + (1-probL)*xpol_H(ind2,t+1);
        Ev2L(ii) = probL*val_L(ind2,t+1) + (1-probL)*val_H(ind2,t+1);
        
    end
    
    x1H_vec = Ex2L-sig*(i1_vec-Ep2H-rH);
    p1H_vec = bet*Ep2H +kappa*x1H_vec;
    [v1H,indvH] = min(0.5*p1H_vec.^2 + lambda*0.5*x1H_vec.^2 + bet*Ev2H);
    
    ppol_H(i_s,t) = p1H_vec(indvH);
    xpol_H(i_s,t) = x1H_vec(indvH);
    ipol_H(i_s,t) = i1_vec(indvH);
    val_H(i_s,t) = v1H;
    
    x1L_vec = Ex2L -sig*(i1_vec-Ep2L-rL);
    p1L_vec = bet*Ep2L + kappa*x1L_vec;
    [v1L,indvL] = min(0.5*p1L_vec.^2 + lambda*0.5*x1L_vec.^2 + bet*Ev2L);
    ppol_L(i_s,t) = p1L_vec(indvL);
    xpol_L(i_s,t) = x1L_vec(indvL);
    ipol_L(i_s,t) = i1_vec(indvL);
    val_L(i_s,t) = v1L;
    
end

% Simulate path: r^n = rL rL rH rH

i_0 = 0;
d_0 = 0;



% t = 1; r^n = rL
T = 1;
indd = St(:,1) == d_0;
indi = St(:,2) == i_0;
indv1 =find(indd+indi == 2);
i_1 = ipol_L(indv1,T);
d_1 = sign(i_1-i_0);



% t = 2: r^n = rH
T = T + 1;
indd = St(:,1) == d_1;
indi = St(:,2) == i_1;
indv2 =find(indd+indi == 2);
v2 = val_H(indv2,T);
i_2 = ipol_H(indv2,T);
d_2 = sign(i_2-i_1);

% t = 3: r^n = rH
T = T + 1;
indd = St(:,1) == d_2;
indi = St(:,2) == i_2;
indv3 = indd+indi == 2;
i_3 = ipol_H(indv3,T);


% for t = 2, and n2 = 0
v2i0  = val2H0(indv2,1);
v2in0 = min(val2H0(indv2,2:nivec));
v205 = val2H0(indv2,2);
v21 = val2H0(indv2,3);
v215 = val2H0(indv2,4);



