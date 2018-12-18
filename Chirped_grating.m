clear;
clc;

L=0.04;            %Length of the fiber grating
n_eff=1.45;        %Effective refractive index
C=10*1e-9;         %Chirp coefficient
N=50;              %Segmentation
M=1501;            %Total points of computation

lamda1=1549;
lamda2=1551;
lamda=linspace(lamda1,lamda2,M)*1e-9; %Wavelength span
delta_lamda=(lamda2-lamda1)/M*1e-9;

r=zeros(1,length(lamda));
phi=zeros(1,length(lamda));
tau=zeros(1,length(lamda));

for k=1:M
    F=[1,0;0,1];  %Initialize transmission matrix
    for i=1:N
        %£¨1£©Uniform
        delta_n_eff=0.00005;
        %£¨2£©Gaussian apodization
        %delta_n_eff=0.00005*exp((-64*(-L/2+i*L/N)^4)/L^4);
        lamda_D=(1550-C*L/2+C*i*L/N)*1e-9;  %The center wavelength of each fiber grating segment
        sigma=2*pi*n_eff*(1/lamda(k)-1/lamda_D)+2*pi*delta_n_eff/lamda(k)+(4*pi*n_eff)*C*(-L/2+i*L/N)/lamda_D^2;  %Self coupling coefficient
        kappa=pi*delta_n_eff/lamda(k);  %AC coupling coefficient
        omega=sqrt(kappa^2-sigma^2);
        
        f11=cosh(omega*L/N)-1i*(sigma/omega)*sinh(omega*L/N);
        f12=1i*(kappa/omega)*sinh(omega*L/N);
        f21=-1i*(kappa/omega)*sinh(omega*L/N);
        f22=cosh(omega*L/N)+1i*(sigma/omega)*sinh(omega*L/N);
        F=F*[f11,f12;f21,f22];  %Transmission matrix
    end
    r(k)=(abs(F(3)/F(1)))^2;  %Magnitude of refractive rate
    phi(k)=phase(F(3)/F(1));  %Phase of refractive rate
end

tau(1)=phi(1);  %Delay
tau(2)=phi(2);
tau(3)=phi(3);

for i=4:M
    if (abs(phi(i-1)-phi(i))<=1)   %Derivative at the phase discontinuities
        tau(i)=((lamda1+i*delta_lamda)^2*1e-18/(2*pi*3e-4)*(phi(i-1)-phi(i))/delta_lamda);
    else
        tau(i)=((lamda1+i*delta_lamda)^2*1e-18/(2*pi*3e-4)*(phi(i-3)-phi(i-2))/delta_lamda);
    end
end

subplot(2,1,1)
plot(lamda*1e9,r,'b');
grid on
xlabel('Wavelength/nm');
ylabel('Reflection');
title('Reflection Spectrum of Chirped Grating');

subplot(2,1,2)
plot(lamda*1e9,tau,'r');
grid on
xlabel('Wavelength/nm');
ylabel('Time Delay');
title('Time Delay of Chirped Grating');