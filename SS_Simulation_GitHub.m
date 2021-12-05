%this function ensures the concentrations of AUXc ,AUXx, CKc, CKx, PINc,
%PINx, MPc, MPx, PXYin, PXYa reach a steady state


function[AUXc ,AUXx, CKc, CKx, PINc, PINx, MPc, MPx, PXYin, PXYa]=SS_Simulation_GitHub(AUXp,...
    Fa, CKp, TDIF, d_Aux, d_PIN, ...
    d_MP, d_CK, d_PXY_in, d_PXY_a, r1, r2, r3, r4, r5, r6, r7, r8)

%starting parameter to be recorded in case a NaN was found
parameter_vector_begin=[AUXp,Fa, CKp, TDIF, d_Aux, d_PIN, ...
    d_MP, d_CK, d_PXY_in, d_PXY_a, r1, r2, r3, r4, r5, r6, r7, r8];

        %% Initial stuff for all the simulations
        AUXc = 0;
        AUXx = 0;
        CKc = 0;
        CKx = 0;
        PINc=0;
        PINx=0;
        MPc=0;
        MPx=0;
        PXYin=0;
        PXYa=0;
        %cell sizes
        dx1      = 1.9;              % um: Ph+cam radius                           
        dx2      = 6;                %ca+xy radius
        dx3      = 5.4;              %xy radius
        Dconst  = 220;               % uM2.s-1: diffusion coefficient. 
        t=0;                         %starting t
        dt=0.0000005;                %time increment
        dR=0.0000005;                %reaction increment
        Diff   = Dconst * [dt/dx1^2, dt/dx2^2, dt/dx3^2]; % Diffusion coefficient 
        %scaled for cell size and diffusion time step. 
        counter=0;
        % Explicit Diffusion matrix.
        Hop = [1, 0, 0; Diff(1), (1-Diff(1)-Diff(2)), Diff(2);...
            0, Diff(2), 1-Diff(2)-Diff(3)];

        %% Ensuring a steady state is achieved
        % holding Conc_old to check for steady state
        Conc_old = AUXc + AUXx + CKc + CKx + PINc + PINx + MPc + MPx + PXYin + PXYa;
        Tno = 3;               % number of seconds difference in states below threshold
        Tcount = 0;             % counting seconds
        atSS =0;
        %%%%%
        rec_sec=0.001;
        % main program
        while atSS == 0   %while steady state is not achieved

            t=t+dt;
             counter=counter+1;
            % Euler thing here
            [AUXc ,AUXx, CKc, CKx, PINc, PINx, MPc, MPx, PXYin, PXYa]=Euler_Function_GitHub(AUXp,...
                                        AUXc, AUXx, CKc, CKx, TDIF, PINc, PINx,...
                                       MPc, MPx, PXYin, PXYa,Fa, d_Aux, d_PIN,...
                                        r1, r2, r3, r4, r5, r6, r7, r8, d_PXY_a,...
                                        d_PXY_in ,d_CK, d_MP, dR, dt, Hop,CKp);

        %check for NaNs
          if (isnan(AUXc)==1)||(isnan(AUXx)==1)||(isnan(CKc)==1)||(isnan(CKx)==1)||(isnan(PINc)==1)||(isnan(PINx)==1)||(isnan(MPc)==1)||(isnan(MPx)==1)||(isnan(PXYa)==1)||(isnan(PXYin)==1)
           savename='Parameters_giving_NaN.csv';
           csvwrite(savename,parameter_vector_begin);
               disp('nan found at parameter_vector1')
             pause()
                %isnan
          end
             %% 
            if(mod(counter,1/dt)==0) % take a snap shot every second
              % count number of consecutive sec mean distance below throshold 
              Conc_now = AUXc + AUXx + CKc + CKx + PINc + PINx + MPc + MPx + PXYin + PXYa;

               if abs((mean(Conc_old-Conc_now))) < 0.0001 
                   Tcount = Tcount + 1;
                % reset Tcount
               elseif Tcount > 0
                    Tcount = 0;
                end

                % stop simulation
                if Tcount == Tno    
                    atSS = 1; 
                end

                Conc_old = Conc_now;
    
            end
        end