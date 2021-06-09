%The goal of this programme is to take paramters AUXp,Fa, CKp, TDIF, d_Aux, d_PIN, ...
%d_MP, d_CK, d_PXY_in, d_PXY_a, r1, r2, r3, r4, r5, r6, r7, r8 and
%outputs final steady state concentrations AUXc ,AUXx, CKc, CKx, PINc, 
%PINx, MPc, MPx, PXYin, PXYa.

%The parameters are read from a file. The programme then runs these parameters
%through the solver and outputs a table of parameters, three columns of 999 
%for visual separation, and the final steady state concentrations of the relevant
%components. 

%The programme uses a function SS_Simulation_GitHub.m to ensure the concentrations reach a
%steady state. At every step dt, Euler_Function_GitHub.m calculates the
%reaction part of the equations using the Euler method. This is then
%updated with a diffusion matrix Hop. The above calclations are repeated 
%until a steady state is reached. 

%% Main programme

clear

%% Read table
Table = load('Full_network_p41.csv','w');

%% Insert parameters from a file
Table(1,:)=[];
%Table(:,19:31)=[];
[n,m]=size(Table);
num_parameters=m;
Table_width=num_parameters+3+10;
Table_for_printing=zeros(Table_width,1)';
%% Specify vectors where results will be recorded

for i=1:length(Table)
%parameter vectors
parameter_vector=Table(i,:);
[AUXc ,AUXx, CKc, CKx, PINc, PINx, MPc, MPx, PXYin, PXYa]=...
               SS_Simulation_GitHub(parameter_vector(1), parameter_vector(2),...
               parameter_vector(3), parameter_vector(4),...
               parameter_vector(5), parameter_vector(6),...
               parameter_vector(7), parameter_vector(8),...
               parameter_vector(9), parameter_vector(10),...
               parameter_vector(11), parameter_vector(12),...
               parameter_vector(13), parameter_vector(14), ...
               parameter_vector(15),parameter_vector(16),...
               parameter_vector(17),parameter_vector(18));
           
           
          %ensure there are no NaNs. If there are, just run a pause
          %algorithm
           if (isnan(AUXc)==1)||(isnan(AUXx)==1)||(isnan(CKc)==1)||(isnan(CKx)==1)||(isnan(PINc)==1)||(isnan(PINx)==1)||(isnan(MPc)==1)||(isnan(MPx)==1)||(isnan(PXYa)==1)||(isnan(PXYin)==1)
           savename='Parameters_giving_NaN.csv';
           csvwrite(savename,parameter_vector);
               disp('nan found at parameter_vector1')
             pause()
           end
          
           %make a vector 'entry' with the parameter, 999 999 999 to visually
           %separate your paramters from your output concentrations
           
             entry=[parameter_vector 999 999 999 AUXc, AUXx,...
            CKc, CKx, PINc, PINx, MPc, MPx, PXYin,PXYa];
           %combine the vectors
           Table_for_printing=vertcat(Table_for_printing,entry);

end

csvwrite('Table_for_printing_full_network_p41.csv',Table_for_printing);