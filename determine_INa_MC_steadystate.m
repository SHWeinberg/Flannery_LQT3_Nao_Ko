function SS = determine_INa_MC_steadystate(V, mutant_type, mutant_flag)

y0 = [1 zeros(1,12)];
[~,y] = ode15s(@(t,y) ode_MC(t, y, V, mutant_type, mutant_flag), [0 10000], y0);
% plot(t, y);
SS = y(end,:);

function dyt = ode_MC(t, y, V, mt, mf)
Npatch = 1;
switch mt
    case 'Y1795C'
        if mf
            UL = .5e-6*ones(Npatch,1);
            LU = 6e-4*ones(Npatch,1);
        else
            UL = 1e-7*ones(Npatch,1);
            LU = 3.8e-3*ones(Npatch,1);
        end
    case 'I1768V'
        UL = 1e-7*ones(Npatch,1);
        LU = 3e-4*ones(Npatch,1);
    case '1795insD'
        if mf
            UL = 1e-7*ones(Npatch,1);
            LU = 9.5e-4*ones(Npatch,1);
        else
            UL = zeros(Npatch,1);
            LU = zeros(Npatch,1);
        end
    case 'dKPQ'
        if mf
            UL = 2e-6*ones(Npatch,1);
            LU = 1e-4*ones(Npatch,1);
        else
            UL = zeros(Npatch,1);
            LU = zeros(Npatch,1);
        end
        
end

alpha11 = a11(V,mf,mt);
beta11 = b11(V,mt);
alpha12 = a12(V,mf,mt);
beta12 = b12(V,mt);
alpha13 = a13(V,mf,mt);
beta13 = b13(V,mf,mt);
alpha2 = a2(V,mf,mt);
alpha3 = a3(V,mf,mt);
beta3 = b3(V,mf,mt);
alpha4 = a4(V,mf,mt);
beta4 = b4(V,mf,mt);
alpha5 = a5(V,mf,mt);
beta5 = b5(V,mf,mt);

% beta2 = b2(V,mf);
beta2 = alpha13.*alpha2.*alpha3./(beta13.*beta3);

states = reshape(y,13,Npatch);
dstates = zeros(13,Npatch);


dstates(1,:) =    beta3.*states(6,:)' + beta11.*states(2,:)' - states(1,:)'.*(alpha3 + alpha11);
dstates(2,:) =  alpha11.*states(1,:)' + beta3.*states(7,:)' + beta12.*states(3,:)' - states(2,:)'.*(alpha3 + alpha12 + beta11);
dstates(3,:) = alpha2.*states(9,:)' + alpha12.*states(2,:)' + beta4.*states(4,:)' + beta3.*states(8,:)' - states(3,:)'.*(alpha3 + alpha4 + beta2 + beta12);
dstates(4,:) =                       alpha4.*states(3,:)' + beta5.*states(5,:)' - states(4,:)'.*(alpha5 + beta4);
dstates(5,:) =                                        alpha5.*states(4,:)' - beta5.*states(5,:)';
dstates(6,:) =              LU.*states(10,:)' + alpha3.*states(1,:)' + beta11.*states(7,:)' - states(6,:)'.*(UL + alpha11 + beta3);
dstates(7,:) = LU.*states(11,:)' + alpha3.*states(2,:)' + alpha11.*states(6,:)' + beta12.*states(8,:)' - states(7,:)'.*(UL + alpha12 + beta3 + beta11);
dstates(8,:) = LU.*states(12,:)' + alpha3.*states(3,:)' + alpha12.*states(7,:)' + beta13.*states(9,:)' - states(8,:)'.*(UL + alpha13 + beta3 + beta12);
dstates(9,:) =             LU.*states(13,:)' + alpha13.*states(8,:)' + beta2.*states(3,:)' - states(9,:)'.*(UL + alpha2 + beta13);
dstates(10,:) =                          UL.*states(6,:)' + beta11.*states(11,:)' - states(10,:)'.*(LU + alpha11);
dstates(11,:) =              UL.*states(7,:)' + alpha11.*states(10,:)' + beta12.*states(12,:)' - states(11,:)'.*(LU + alpha12 + beta11);
dstates(12,:) =          UL.*states(8,:)' + alpha12.*states(11,:)' + beta13.*states(13,:)' - states(12,:)'.*(LU + alpha13 + beta12);
dstates(13,:) =                      UL.*states(9,:)' + alpha13.*states(12,:)' - states(13,:)'.*(LU + beta13);

dyt = reshape(dstates,13*Npatch,1);


function x = a11(V,mflag, mt)
switch mt
    case {'Y1795C','I1768V','1795insD'}
        x = 3.802./(0.1027*exp(-V/17) + 0.2*exp(-V/150));
    case 'dKPQ'
        x = (1+.25*mflag)*(3.002 + 0.80*10^(-mflag))./(0.1027*exp(-V/17) + 0.2*exp(-V/150));
end


function x = a12(V, mflag, mt)
switch mt
    case {'Y1795C','I1768V','1795insD'}
        x = 3.802./(0.1027*exp(-V/15) + 0.23*exp(-V/150));
    case 'dKPQ'
        x = (1+.25*mflag)*(3.002 + 0.80*10^(-mflag))./(0.1027*exp(-V/15) + 0.23*exp(-V/150));
        
end


function x = a13(V, mflag,mt)
switch mt
    case {'Y1795C','I1768V','1795insD'}
        x = 3.802./(0.1027*exp(-V/12) + 0.25*exp(-V/150));
    case 'dKPQ'
       x = (1+.25*mflag)*(3.002 + 0.80*10^(-mflag))./(0.1027*exp(-V/12) + 0.25*exp(-V/150));
end


function x = b11(V, mt)
switch mt
    case 'Y1795C'
        x = 0.77*exp(-V/20.3);
    case 'I1768V'
        x = 0.4*exp(-V/20.3);
    case {'1795insD','dKPQ'}
        x = 0.1917*exp(-V/20.3);
end


function x = b12(V, mt)
switch mt
    case 'Y1795C'
        x = 0.77*exp(-(V-5)/20.3);
    case 'I1768V'
        x = 0.4*exp(-(V-5)/20.3);
    case {'1795insD','dKPQ'}
        x = 0.2*exp(-(V-5)/20.3);
end


function x = b13(V,mflag, mt)
switch mt
    case 'Y1795C'
        if mflag
            x = 0.535*exp(-(V-10)/20.3);
        else
            x = 0.17*exp(-(V-10)/20.3);
        end
    case 'I1768V'
        x = 0.4*exp(-(V-10)/20.3)/4.5;
    case {'1795insD','dKPQ'}
        x = 0.22*exp(-(V-10)/20.3);
end



function x = b2(V,mflag, mt)
x = a13(V, mflag,mt).*a2(V,mflag, mt).*a3(V, mflag, mt)./(b13(V,mflag, mt).*b3(V, mflag, mt));




function x = a3(V, mflag, mt)
switch mt
    case 'Y1795C'
        x = 3.7933e-7*exp(-(V)/7.7);
    case 'I1768V'
        x = (1+mflag)*1.897e-6*exp(-(V)/7.7);
    case '1795insD'
        x = 3.7933e-7*exp(-V/7.7)/(1+1.5*mflag);
    case 'dKPQ'
        x = (1+mflag*19)*3.7933e-9*exp(-V/5.2);
end


function x = b4(V, mflag, mt)
switch mt
    case {'Y1795C','1795insD','dKPQ'}
        x = a3(V, mflag, mt);
    case 'I1768V'
        x = 5*a3(V, mflag, mt);
end


function x = b5(V,mflag, mt)
switch mt
    case 'Y1795C'
        x = a3(V,mflag,mt)/50;
    case 'I1768V'
        x = 5*a3(V,mflag,mt)/50;
    case '1795insD'
        x = a3(V,mflag,mt)/(20 + ~mflag*30);
    case 'dKPQ'
        x = 0;
end




function x = b3(V, mflag,mt)
switch mt
    case {'Y1795C','I1768V','1795insD'}
        x = 0.0084 + 0.00002*V;
    case 'dKPQ'
        x = (0.0084 + 0.00002*V)*(1+mflag);
end



function x = a2(V,mflag, mt)
switch mt
    case 'Y1795C'
        x = 2.04*exp((V)/29.68)/(1+mflag);
    case 'I1768V'
        x = 9.178*exp((V)/29.68)/4.5;
    case '1795insD'
        x = 9.178*exp(V/29.68);
    case 'dKPQ'
        x = 9.178*exp(V/(29.68+mflag*70.32));
end


function x = a4(V,mflag, mt)
switch mt
    case {'Y1795C','1795insD','dKPQ'}
        x = a2(V,mflag, mt)/100;
    case 'I1768V'
        x = 1.5*a2(V,mflag,mt)/100;
        
end


function x = a5(V,mflag, mt)
switch mt
    case 'Y1795C'
        x =  a2(V,mflag,mt)/9.5e4;
    case 'I1768V'
        x = 1.5*a2(V,mflag,mt)/9.5e4;
    case '1795insD'
        x = a2(V,mflag,mt)/(3.5e4 + ~mflag*6e4);
    case 'dKPQ'
        x = 0;
end



