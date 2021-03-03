%%%% Solution of Dust Charging %%%%%
%%% Derivation of Dispersion relation %%%%

%  Variation with Dust Temperature
   clc;clear; 

% Defined values 
  ep0 = 8.85e-12;
  ev = 1e-16;
  W = 3*ev;                       % Work function
  J = 4e22;                        % Photon flux density
  Tp = 11.62*ev;                  % Average energy of photo electrons
  Y = 1e-6;                          % yield of photon
  r = 0.1e-6;                       % dust grain radius
  Te = 0.5*ev;
  Ti = 0.33*ev;
  Td = 0.024*ev;
  ni0 = 5e9;
  ne0 = 2.5e9;
  nd0 = 1.25e6;
  zd0 = 2e3;
  md = 1.15e-16;
  mi = 6.6e-26;
  mn = 6.6e-26;
  me = 9.1e-31;
  lambda = 0.05; 
  
  e = 1.6e-19; 
  h = 6.6e-34; 
 sigma_d = [0.04,0.03,0.02];
      nu_dn = 50;
       nn = 3.14e10;                     % neutral density
   
% Calculated values %
 for i = 1:3
  Nd0 = nd0/nd0;
  alpha = (4/(ni0*(((8*Te)/(pi*me))^(1/2))));
  gamma = ((16*pi*me*Td^2)/(h^3*(((8*Te)/(pi*me))^(1/2))*ni0))*exp(-W/Td);
  Z =  (e^2*zd0)/(r*Te);
  K = ((me*Ti)/(mi*Te))^(1/2);
  beta = Tp/Te;
  sigma_i = Ti/Te;
  sigma_d_1 = Td/Te;
  delta = ne0/ni0;
  Delta = sqrt(ni0/(nd0*zd0));
  R = ((((Z/sigma_i)-1)*(K/sigma_i))-(delta*(1+Z)))/(Z*(delta+(K/sigma_i)...
      + ((2*gamma*Z)/sigma_d(i))+((Y*J*alpha)/beta)+(gamma*((Z/sigma_d(i))-1))));
  P = (Nd0*R)/(Delta)-(((Delta*Nd0)/sigma_i)*(1+(delta*sigma_i)));  
  Q = (3*sigma_d(i)*Delta^2);

    w_pd = sqrt((nd0*e^2*zd0^2)/(ep0*md));
    lambda_eB = sqrt((ep0*Te)/(ne0*e^2));
%     nu_dn = (4./3).*pi.*(a.^2).*vn.*nn.*(mn./md); 

    ks = linspace(0.001,8,80);
    t = ks.*ks;
    t1 = ks.*ks.*lambda_eB.*lambda_eB;
    for k = 1:length(ks)
        
        A = 1 - (P/(ks(k)^2*Nd0));
        B = -((Delta*nu_dn)/(pi))+((P*nu_dn)/(pi*ks(k)^2*Nd0*delta));
        C = (Nd0)- ((3*sigma_d(i)*Delta*P)/(Nd0))+(Q*ks(k)^2);
        
        x1(k) = (C/(A));
        X1(k) = sqrt(C/(A));
        
        Y1(k) = (3*B)/(2*A);
    end
    
    
    hold all;
    figure(1) 
    plot(ks,(X1),'linewidth',2);
    legend('\sigma_d - 0.04','\sigma_d - 0.03','\sigma_d - 0.02');
    grid on
%     xlabel('\kappa^2');
%     ylabel('\omega_r^2');
     xlabel('\kappa');
    ylabel('\omega_r');


    
%     hold all;
%     figure(2) 
%     plot(ks,(Y1));
% %     hleg1 = legend('Td - 1/40','Td - 3','Td -  5','Td -10');
%     xlabel('\kappa');
%     ylabel('\omega_i');


%     color = 'rmbk' 
%     markers  = '+p*d'
%     hold all;
%     figure(1) 
%     plot(ks,real(x1),[color(i),markers(i)]);
%     hleg1 = legend('Td - 1/40','Td - 3','Td -  5','Td -10');
%     xlabel('k');
%     ylabel('w');
  end