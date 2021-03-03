clc;clear 
global sigma_i Td Te nu_dn delta Delta zd

 ev = 1.6e-19; r = 1e-8;Te = 0.5*ev; mi = 6.6e-26; me = 9.1e-31; Ti = 0.3*ev; W = 3*ev;
 ne0 =  5e10;ni0 =  2.5e10;nd0 = 1.25e7; e = 1.6e-19; h = 6.6e-34; Y =0.5;J = 4e22;Tpe = 11.2*ev;          % mars
 zd0 = 2000; nu_dn = 1.65e2; Td = 0.05;

zd = (e^2*zd0)/(r*Te);delta = ne0/ni0;sigma_i = Ti/Te;Delta = sqrt(ni0/(nd0*zd));alpha = 4/(ni0*sqrt((8*Te)/(pi*me)));
K = sqrt((me*Ti)/(mi*Te)); beta = Tpe/Te;  


xspan = linspace(0.01,0.1,1000);
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[x,y] = ode45('sub_dust',xspan,[0.0 0.0 1 1e-6],options);

t = y(:,3).*y(:,4);

%%%%%%%%%% Calculation of surface potential %%%%%%%%%%%%%%%%%%%%

% W = [1*ev,2*ev,3*ev];
Td2 = linspace(0.001,1,100);            % experiment
Td1 = Td2.*ev;
%  Td = linspace(0.01,0.1,10);         % Martian temp
for j = 1:3
    td1 = length(Td1);
     for i = 1:td1
    %      zd1 = 200;
        gamma = (16*pi*me*Td1(i)^2*exp(-W/Td1(i)))/(h^3*sqrt((8*Te)/(pi*me))*ni0);
        func = @(phi_d) -delta*exp(y(:,1))*(1+ ((e*phi_d*beta)/Te))+(K*exp(-y(:,1)/sigma_i)*exp(-((e*phi_d*Tpe)/(Te*Ti))))+gamma *(1+((e*phi_d)/Te))*exp(-e*phi_d/Td1(i))...
                     + (Y*J*alpha*exp(-e*phi_d/Te));    

%         func_1 = @(z) -(delta*exp(y(:,1))*(1+ ((e^2*zd1*z)/(r*Te))))+(K*exp(-y(:,1)/sigma_i)*exp(-((e^2*zd1*z)/(r*Ti))))+gamma *(1+((e^2*zd1*z)/(r*Te)))*exp(-((e^2*zd1*z)/(r*Td1(i))))...
%                      + (Y*J*alpha*exp(((e^2*zd1*z)/(r*Tpe))));          

        phi_d0 = 1;
%         z0 = 10;
        phi_d1(i)  = fsolve(func,phi_d0);

%         z1(i)  = fsolve(func,z0);

     end
     
    figure(1)
    hold all
    plot(Td2,(phi_d1),'linewidth',2)
%     legend('W - 7eV','W - 4eV','W - 1eV')
    grid on
    ylabel('\phi_d')
%     xlim([0,0.35])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(2)
% subplot(2,2,1);plot(x,y(:,1));
% ylabel('\phi')
% subplot(2,2,2);plot(x,y(:,2));
% ylabel('E')
% subplot(2,2,3);plot(x,y(:,3));
% ylabel('N_d')
% subplot(2,2,4);plot(x,y(:,4));
% ylabel('v_d')

% 
%     p = diff(y(:,2));
% 
% figure(3)
% plot(x,t)

