clear;
clc;

n_eff=1.46; %Effective refractive index
lamda_B=1550*1e-9;  %Center wavelength
L=0.1;  %Length of the fiber grating
delta_n_eff=2*1e-5; %Refractive rate modulation depth

phi=input('Phase shift is: ');
T=input('The length ratio is: ');
lamda=1e-9*linspace(1549.85,1550.15,500);  %Wavelength span

delta=zeros(1,length(lamda));
s=zeros(1,length(lamda));
s11=zeros(1,length(lamda));
s12=zeros(1,length(lamda));
s21=zeros(1,length(lamda));
s22=zeros(1,length(lamda));
ss11=zeros(1,length(lamda));
ss12=zeros(1,length(lamda));
ss21=zeros(1,length(lamda));
ss22=zeros(1,length(lamda));
R=zeros(1,length(lamda));
r=zeros(1,length(lamda));

z_i_1= L/(1+T); %Uniform fiber grating length before phase shift
z_i_2=L- z_i_1; %Uniform fiber grating length after phase shift
fpi=[exp(-1i*phi/2),0;0, exp(1i*phi/2)];  %Phase shift matrix
F=[1 0;0 1];    %Initialize transmission matrix

for num=1:500
    delta(num)=2*pi*n_eff*(1./lamda(num)-1./lamda_B);%Inter-mode detuning
    kappa=pi./lamda(num)*delta_n_eff;
    s(num)=sqrt(kappa.^2-delta(num).^2);
    
    s11(num)=(cosh(s(num)*z_i_1)-1i*(delta(num)/s(num))*sinh(s(num)*z_i_1));
    s12(num)=-1i*(kappa./s(num)).*sinh(s(num)*z_i_1);
    s21(num)=1i*(kappa./s(num)).*sinh(s(num)*z_i_1);   
    s22(num)=cosh(s(num)*z_i_1)+1i*(delta(num)/s(num))*sinh(s(num)*z_i_1);
    f1=[s11(num) s12(num);s21(num) s22(num)];

    ss11(num)=(cosh(s(num)*z_i_2)-1i*(delta(num)/s(num))*sinh(s(num)*z_i_2));
    ss12(num)=-1i*(kappa./s(num)).*sinh(s(num)*z_i_2);
    ss21(num)=1i*(kappa./s(num)).*sinh(s(num)*z_i_2);   
    ss22(num)=cosh(s(num)*z_i_2)+1i*(delta(num)/s(num))*sinh(s(num)*z_i_2);
    f2=[ss11(num) ss12(num);ss21(num) ss22(num)];

    f_out=f2*fpi*f1*F;
    r(num)=f_out(3)/f_out(1);
    R(num)=(abs(r(num)))^2; %Amplitude of the refractive index
end

plot(lamda*1e9,R);
grid on;
xlabel('Wavelength/nm');
ylabel('Reflection');