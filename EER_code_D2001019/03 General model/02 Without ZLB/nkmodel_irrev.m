function [y,p] = nkmodel_irrev(policyrate,g,ilag,gnext,Fy,Fp,params,expands,prb)
% Codes for "Optimal irreversible monetary policy"
% Copyright(c) Kohei Hasui, Matsuyama University (Feb 2021).
% Note: Please use this code your own responsibility. I would not be responsible for any damage or error that might occur by using this code.
% Matlab function: Evaluating inflation rate and output gap.
%  Input:
%        policyrate: nominal interest rate (running var)
%        g:          demand shock (state var)
%        ilag:       lagged nominal interest rate (state var)
%        gnext:      tomorrow possible demand shock 
%        prb:        probability of state
%  Output:
%        y: output gap
%        p: inflation rate
%

% extract params
bet = params(2);
sig = params(3);
kpp = params(7);

dlti = policyrate-ilag; % Di
% delta : -1 0 1
dlt  = sign(dlti);

% expected inflation and output gap
Ep = Fp(gnext,policyrate*expands,dlt*expands)*prb;
Ey = Fy(gnext,policyrate*expands,dlt*expands)*prb;

y = Ey - sig*(policyrate-Ep) + g;
p = bet*Ep + kpp*y;
