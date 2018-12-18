function varargout = FBG_GUI(varargin)
% FBG_GUI MATLAB code for FBG_GUI.fig
%      FBG_GUI, by itself, creates a new FBG_GUI or raises the existing
%      singleton*.
%
%      H = FBG_GUI returns the handle to a new FBG_GUI or the handle to
%      the existing singleton*.
%
%      FBG_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FBG_GUI.M with the given input arguments.
%
%      FBG_GUI('Property','Value',...) creates a new FBG_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FBG_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FBG_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FBG_GUI

% Last Modified by GUIDE v2.5 18-Dec-2018 11:33:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FBG_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FBG_GUI_OutputFcn, ...
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


% --- Executes just before FBG_GUI is made visible.
function FBG_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FBG_GUI (see VARARGIN)

% Choose default command line output for FBG_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FBG_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FBG_GUI_OutputFcn(~, ~, handles) 
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



function edit3_Callback(~, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(~, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, ~, ~)
% hObject    handle to edit4 (see GCBO)
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

delta_n_eff=str2num(get(handles.edit1,'string'));
L=str2num(get(handles.edit2,'string'));

a=4.5e-6;   %Radius
n1=1.449;   %Core refractive index
n2=1.444;   %Cladding refractive index
lamda_Brag=1550e-9; %Center wavelength
%L=0.2;      %Length of the fiber grating
lamda=1e-9*linspace(1549.8,1550.2,500); %Wavelength span
%delta_n_eff=2*1e-5; %Refractive rate modulation depth
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
set(handles.edit3,'string',n_eff);

aa=besselj(0,US)./besselk(0,WS); %Core and cladding bessel function relative coefficient. 
%e=besselj(0,US*r/a) for core
%e=aa*besselk(0,US*r/a) for cladding

%Integral and normalization
%对FBG，模式的交叠积分就是自身的交叠积分
e21=integral(@(r)(besselj(0,US.*r./a)).^2.*r,0,a);
e22=integral(@(r)(aa.*besselk(0,WS.*r./a)).^2.*r,a,62.5e-6);
gamma=e21./(e21+e22);
set(handles.edit4,'string',gamma);

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
    
    a(k)=sinh(sqrt(kappa(k)^2-sigma1(k)^2)*L)^2;
    b(k)=cosh(sqrt(kappa(k)^2-sigma1(k)^2)*L)^2-sigma1(k)^2/kappa(k)^2;
    r(k)=a(k)/b(k);  %Refractive index
    t(k)=1-r(k);
end

plot(handles.axes1,lamda*1e9,r,'b');
hold on
grid on
plot(handles.axes1,lamda*1e9,t,'r');
title('FBG Reflection Spectrum');
xlabel('Wavelength/nm');
ylabel('Reflection & Transmission');
legend('Reflection','Transmission');
