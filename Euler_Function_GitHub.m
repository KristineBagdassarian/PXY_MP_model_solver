%This function calls a function to calculate the reaction changes using
%Euler method. This is then updated with a diffusion matrix Hop

function[AUXc ,AUXx, CKc, CKx, PINc, PINx, MPc, MPx, PXYin, PXYa]=Euler_Function_GitHub(AUXp,...
                                AUXc, AUXx, CKc, CKx, TDIF, PINc, PINx,...
                                MPc, MPx, PXYin, PXYa,Fa, d_Aux, d_PIN,...
                                r1, r2, r3, r4, r5, r6, r7, r8, d_PXY_a,...
                                d_PXY_in ,d_CK, d_MP, dR, dt, Hop,CKp)

                   
           %solve reaction part using the Euler method                 
     for i = 1:dt/dR  
        [AUXc ,AUXx, CKc, CKx, PINc, PINx, MPc, MPx, PXYin, PXYa] = ...
        Reaction_Part_GitHub(AUXp, AUXc, AUXx, CKc, CKx, TDIF, PINc, PINx,...
        MPc, MPx, PXYin, PXYa,Fa, d_Aux, d_PIN, r1, r2, r3, r4, r5,...
        r6, r7, r8, d_PXY_a, d_PXY_in ,d_CK, d_MP, dR);  
     end % End looping over the reaction time steps
     
    % Implicit diffusion 
    % apply the diffusion matrix
    CK = [CKp,CKc,CKx]; %vector that contains concen of CK in diff cells 
    
    CK = (Hop*CK')'; 
    CKc = CK(2); 
    CKx = CK(3);
    