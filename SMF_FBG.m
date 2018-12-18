clear;
clc;

a=4.5e-6;   %Radius
n1=1.449;   %Core refractive index
n2=1.444;   %Cladding refractive index
lamda_Brag=1550e-9; %Center wavelength
L=0.2;      %Length of the fiber grating
lamda=1e-9*linspace(1549.8,1550.2,500); %Wavelength span
delta_n_eff=2*1e-5; %Refractive rate modulation depth
v=1;    %Fringe visibility

k0=2*pi/lamda_Brag;
V=k0*a*sqrt(n1^2-n2^2);
U=0.5:0.01:V; %U^2+W^2=V^2, in order to make sure there is no imaginary part, set U<=V
W=sqrt(V.^2-U.^2);
f=besselj(0,U)./U./besselj(1,U)-besselk(0,W)./W./besselk(1,W);

for i=1:(length(U)-1)
    if f(i)*f(i+1)<=0
        U=U(i):1e-4:U(i+1);
        W=sqrt(V.^2-U.^2);
        f=besselj(0,U)./U./besselj(1,U)-besselk(0,W)./W./besselk(1,W);
        for j=1:(length(U)-1)
            if f(j)*f(j+1)<=0
                US=(U(j)+U(j+1))/2;
                break
            end
        end
        break
    end
end

WS=sqrt(V.^2-US.^2);
beta=sqrt(k0.^2.*n1.^2-US.^2./a.^2);
n_eff=beta./k0;
aa=besselj(0,US)./besselk(0,WS); %Core and cladding bessel function relative coefficient. 
%e=besselj(0,US*r/a) for core
%e=aa*besselk(0,US*r/a) for cladding

%Integral and normalization
%对FBG，模式的交叠积分就是自身的交叠积分
e21=integral(@(r)(besselj(0,US.*r./a)).^2.*r,0,a);
e22=integral(@(r)(aa.*besselk(0,WS.*r./a)).^2.*r,a,62.5e-6);
gamma=e21./(e21+e22);

sigma=zeros(1,length(lamda));
kappa=zeros(1,length(lamda));
delta=zeros(1,length(lamda));
sigma1=zeros(1,length(lamda));
a=zeros(1,length(lamda));
b=zeros(1,length(lamda));
r=zeros(1,length(lamda));
t=zeros(1,length(lamda));

for k=1:500
    sigma(k)=pi./lamda(k)*n1*delta_n_eff*gamma;  %DC coupling coefficient
    kappa(k)=v/2*sigma(k);  %AC coupling coefficient
    delta(k)=2*pi*n_eff.*(1./lamda(k)-1/lamda_Brag);  %Detuning
    sigma1(k)=delta(k)+sigma(k);  %DC self coupling coefficient
    
    a(k)=sinh(sqrt(kappa(k).^2-sigma1(k).^2)*L).^2;
    b(k)=cosh(sqrt(kappa(k).^2-sigma1(k).^2)*L).^2-sigma1(k).^2/kappa(k).^2;
    r(k)=a(k)/b(k);  %Refractive index
    t(k)=1-r(k);
end

plot(lamda*1e9,r,'b')
hold on
grid on
plot(lamda*1e9,t,'r')
title('FBG Reflection Spectrum');
xlabel('Wavelength/nm');
ylabel('Reflection & Transmission');
legend('Reflection','Transmission');