function varargout = Chirped_grating_GUI(varargin)
% CHIRPED_GRATING_GUI MATLAB code for Chirped_grating_GUI.fig
%      CHIRPED_GRATING_GUI, by itself, creates a new CHIRPED_GRATING_GUI or raises the existing
%      singleton*.
%
%      H = CHIRPED_GRATING_GUI returns the handle to a new CHIRPED_GRATING_GUI or the handle to
%      the existing singleton*.
%
%      CHIRPED_GRATING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHIRPED_GRATING_GUI.M with the given input arguments.
%
%      CHIRPED_GRATING_GUI('Property','Value',...) creates a new CHIRPED_GRATING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Chirped_grating_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Chirped_grating_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Chirped_grating_GUI

% Last Modified by GUIDE v2.5 18-Dec-2018 12:14:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Chirped_grating_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Chirped_grating_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Chirped_grating_GUI is made visible.
function Chirped_grating_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Chirped_grating_GUI (see VARARGIN)

% Choose default command line output for Chirped_grating_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Chirped_grating_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Chirped_grating_GUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(~, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(~, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla reset;

C=str2num(get(handles.edit1,'string'));
L=str2num(get(handles.edit2,'string'));

n_eff=1.45;        %Effective refractive index
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

axes(handles.axes1)
plot(handles.axes1,lamda*1e9,r,'b');
grid on
xlabel('Wavelength/nm');
ylabel('Reflection');
title('Reflection Spectrum of Chirped Grating');

axes(handles.axes2)
plot(handles.axes2,lamda*1e9,tau,'r');
grid on
xlabel('Wavelength/nm');
ylabel('Time Delay');
title('Time Delay of Chirped Grating');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(~, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla reset;

C=str2num(get(handles.edit1,'string'));
L=str2num(get(handles.edit2,'string'));

n_eff=1.45;        %Effective refractive index
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
        %£¨2£©Gaussian apodization
        delta_n_eff=0.00005*exp((-64*(-L/2+i*L/N)^4)/L^4);
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

axes(handles.axes1)
plot(handles.axes1,lamda*1e9,r);
grid on
xlabel('Wavelength/nm');
ylabel('Reflection');
title('Reflection Spectrum of Chirped Grating');

axes(handles.axes2)
plot(handles.axes2,lamda*1e9,tau,'r');
grid on
xlabel('Wavelength/nm');
ylabel('Time Delay');
title('Time Delay of Chirped Grating');
