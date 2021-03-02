function v = objective_irrev(policyrate,g,ilag,dltlag,gnext,Fv,Fy,Fp,params,expands,prb)
% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.
%  Input:
%        policyrate: nominal interest rate (running var)
%        g: demand shock (state var)
%        ilag: lagged nominal interest rate (state var)
%        gnext: tomorrow possible demand shock 
%        prb: probability of state
%  Output:
%        v: value
%

% extract params
bet = params(2);
sig = params(3);
kpp = params(7);
lamda  = params(8);
lamdaf = params(9);

dlti   = policyrate-ilag;
% delta : -1 0 1
dlt = sign(dlti);

% expected value, inflation, and output gap
Ev = Fv(gnext,policyrate*expands,dlt*expands)*prb;
Ep = Fp(gnext,policyrate*expands,dlt*expands)*prb;
Ey = Fy(gnext,policyrate*expands,dlt*expands)*prb;

y = Ey - sig*(policyrate-Ep) + g;
p = bet*Ep + kpp*y;

% penalty function
f = 1/2*dltlag*(1+dltlag)*(min(dlti,0))^2-1/2*dltlag*(1-dltlag)*(max(dlti,0))^2;

% Bellman eq
v = p^2 + lamda*y^2 + lamdaf*f + bet*Ev;
