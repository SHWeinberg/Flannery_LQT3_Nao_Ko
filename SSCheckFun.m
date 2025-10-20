function SSCheckBIN = SSCheckFun(APD901, APD902, CaTA1, CaTA2)

% calculate percent change beat-to-beat
dAPD90 = (APD902 - APD901)/APD901*100;
dCaTA = (CaTA2 - CaTA1)/CaTA1*100;

% beat-to-beat < 1%
if abs(dAPD90) < 1 && abs(dCaTA) < 1
    SSCheckBIN = 1;
else
    SSCheckBIN = 0; 
end

end