%This function employs the Euler method
function [auxc, auxx, ckc, ckx, pinc, pinx, mpc, mpx, pxyin, pxya] =...
    Reaction_Part_GitHub(AUXp, AUXc, AUXx, CKc, CKx, TDIF, PINc, PINx, MPc,...
    MPx, PXYin, PXYa, Fa, d_Aux, d_PIN, r1, r2, r3, r4, r5, r6, r7, r8,...
    d_PXY_a, d_PXY_in, d_CK ,d_MP, dR) 

    %We use explicit Euler to solve the reaction equations
    %For all species, we use the concentration of species X at time t 
    %to calculate the concentration of species at the next time step, t+1
    

    %%Aux
    auxc = AUXc+dR*( - d_Aux*AUXc+Fa*AUXp-r8*PINc*AUXc+0.5*r8*PINx*AUXx);
    auxx = AUXx+dR*( - d_Aux*AUXx+0.5*r8*PINc*AUXc-r8*PINx*AUXx);
    
    %%Ck
    ckc  =CKc+dR*(-d_CK*CKc-r1*AUXc*CKc); %without diffusion term
    ckx  =CKx+dR*(-d_CK*CKx-r1*AUXx*CKx); %without diffusion term
    
    %%PIN1
    pinc =PINc+dR*(r6*MPc-r7*CKc*PINc-d_PIN*PINc);
    pinx = PINx + dR*(r6*MPx - r7*CKx*PINx-d_PIN*PINx);
    
    
    %%PXYin
    pxyin = PXYin+dR*(r4*MPc-d_PXY_in*PXYin-r2*PXYin*TDIF);
    
    %%PXYa
    pxya = PXYa+dR*(r2*PXYin*TDIF-d_PXY_a*PXYa);
    
    
    %%MP
    mpc = MPc+dR*(r5*AUXc-r3*PXYa*MPc-d_MP*MPc);
    mpx = MPx+dR*(r5*AUXx-d_MP*MPx);
    
    
end