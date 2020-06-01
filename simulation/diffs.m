function dres = diffs(t,res)
J02 = -1082.63*10^-6;
ae = 6378136; 
GM = 398600441.8e6; 
Xate=res(1);
Yate=res(2);
Zate=res(3);
r=sqrt(Xate^2 + Yate^2 + Zate^2);
GM1 = GM/r^2;
x01 = Xate/r;
y01 = Yate/r;
z01 = Zate/r;
p = ae/r;

dres = res(:);
dres(1) = res(4);
dres(2) = res(5);
dres(3) = res(6);

dres(4) = -GM1*x01 + 1.5*J02*GM1*x01*(p^2)*(1 - 5*z01^2);
dres(5) = -GM1*y01 + 1.5*J02*GM1*y01*(p^2)*(1 - 5*z01^2);
dres(6) = -GM1*z01 + 1.5*J02*GM1*z01*(p^2)*(3 - 5*z01^2);
end