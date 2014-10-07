% CDF_analyzer.m
%
% Calculates Gaussian CDF fitting at different localization errors
%
% Note:  
%      it needs the variable *msd* (after excuting the function
%      MSD_calculator_light.m) to function properly.
%
% OUTPUT:
%      rrr: j-cell at different localization error
%           each cell contains 1*5 cells, which are
%           1,1: 2-component gaussian fitting, major component (D1) fraction
%           1,2: 2-component gaussian fitting, D1
%           1,3: 2-component gaussian fitting, D2
%           1,4: 1-component gaussian fitting, D
%           1,5: error sum ratio of the two fitting methods (1-component/2-component)
%
%
% Tzu-Chi Yen
% last modified Sep 2014
%
 
%=============================================================
%                Diffusion mode analysis
%=============================================================    
%     % for simulated Brownian motion
%     msd_b =  MSD_calculator_light({traj});
%     r_sq_1ms = msd_b{1}{2}{1}.^2+msd_b{1}{3}{1}.^2;
    

for j=3:3;
   loc_error = j*0.001;    % unit: micron 
    
for i=1:1
    time_step=0.0005;                % inverse frame-rate, unit:s

    % for experimental data
    r_sq_1ms = msd{1}{2}{i}.^2+msd{1}{3}{i}.^2;
    %r_sq_1ms = sqrt(msd{1}{2}{1}.^2);
    
    
        r_sq = 0:1e-6:max(r_sq_1ms)+1e-4;
        s_l=length(r_sq_1ms);

        hist_r_sq_1ms = hist(r_sq_1ms,r_sq);
        PDF_r_sq_1ms = cumsum(hist_r_sq_1ms)./s_l; 
        options=optimset('display','iter-detailed','TolX',1e-20,'TolFun',1e-20,'MaxFunEvals',1000);

        fun = @(c,r) 1-exp(-r_sq/(4*loc_error^2+c(1)^2));
        c0 = 0.2;
        [cc_x_1 cc_x_1_resnorm cc_x_1_residual]=lsqcurvefit(fun,c0,r_sq,PDF_r_sq_1ms,0,1,options);
        PDF_fit_1=fun(cc_x_1,r_sq);
        Err_1 = sum((PDF_r_sq_1ms-PDF_fit_1).^2);

        fun = @(c,r) 1-c(1)*exp(-r_sq/(4*loc_error^2+c(2)^2))-(1-c(1))*exp(-r_sq/(4*loc_error^2+c(3)^2));
        c0=[0.5 0.2 0.05];
        [cc_x_2 cc_x_2_resnorm cc_x_2_residual]=lsqcurvefit(fun,c0,r_sq,PDF_r_sq_1ms,[0 0 0],[1,1,1],options);
        PDF_fit_2=fun(cc_x_2,r_sq);
        Err_2 = sum((PDF_r_sq_1ms-PDF_fit_2).^2);
              
   
        res_fit{i,1}=PDF_fit_1;
        res_fit{i,2}=PDF_fit_2;
        res_fit{i,3}=PDF_r_sq_1ms;     
        
        

        cc_x_1=cc_x_1^2/4/time_step/i;                 %D1	
        cc_x_2(1,2)=cc_x_2(1,2)^2/4/time_step/i;       %D1-2-component
        cc_x_2(1,3)=cc_x_2(1,3)^2/4/time_step/i;       %D2-2-component

        
        resul(i,1)=cc_x_2(1,1);
        resul(i,2)=cc_x_2(1,2);
        resul(i,3)=cc_x_2(1,3);
        resul(i,4)=cc_x_1;
        resul(i,5)=Err_1/Err_2;
        
%         
end

rrr{j}=resul;
resul
end
%         res_r_sq{i}=r_sq;
% 
%         res_resid{i,1}=cc_x_1_residual;
%         res_resid{i,2}=cc_x_2_residual;
% 
%         res_resnorm(i,1)=cc_x_1_resnorm;
%         res_resnorm(i,2)=cc_x_2_resnorm;
% 
% 
%         res_fit{i,1}=PDF_fit_1;
%         res_fit{i,2}=PDF_fit_2;
%         res_fit{i,3}=PDF_r_sq_1ms;


    %     hold on
    % plot(r_sq,PDF_r_sq_1ms);
    % plot(r_sq,PDF_fit_1,'r-');
    % plot(r_sq,PDF_fit_2,'r-','Color','g');

