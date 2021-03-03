function dy = sub_dust(~,y)
global sigma_i Td Te nu_dn delta Delta zd
     dy = ones(4,1);
     
%      phi = y(1);
%      E = y(2);
%      Nd = y(3);
%      vd = y(4);
     
    
     dy(1) = -y(2);
     dy(2) = exp(-y(1)/sigma_i)+(delta*exp(y(1)))- (zd*y(3)/Delta^2);
     dy(3) = ((y(3)*zd)/((y(4).*y(4))-((3*Td)/Te)))*y(2)-((3*nu_dn*y(3))/((1-((3*Td)/(Te*y(4).*y(4))))*pi));
     dy(4) = -(((zd/y(4))+ ((3*Td*zd)/((Te*y(4).*y(4).*y(4))-(3*Td*y(4)))))*y(2))-((9*Td*nu_dn)/(Te*y(4)*pi*(1-((3*Td)/(Te*y(4).*y(4))))))...
             - ((3*nu_dn*y(4))/pi);
%      dy(4)=34;
%      disp([Nd vd]);
end