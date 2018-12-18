function varargout = Phase_shift_grating_GUI(varargin)
% PHASE_SHIFT_GRATING_GUI MATLAB code for Phase_shift_grating_GUI.fig
%      PHASE_SHIFT_GRATING_GUI, by itself, creates a new PHASE_SHIFT_GRATING_GUI or raises the existing
%      singleton*.
%
%      H = PHASE_SHIFT_GRATING_GUI returns the handle to a new PHASE_SHIFT_GRATING_GUI or the handle to
%      the existing singleton*.
%
%      PHASE_SHIFT_GRATING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHASE_SHIFT_GRATING_GUI.M with the given input arguments.
%
%      PHASE_SHIFT_GRATING_GUI('Property','Value',...) creates a new PHASE_SHIFT_GRATING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Phase_shift_grating_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Phase_shift_grating_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Phase_shift_grating_GUI

% Last Modified by GUIDE v2.5 18-Dec-2018 12:37:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Phase_shift_grating_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Phase_shift_grating_GUI_OutputFcn, ...
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


% --- Executes just before Phase_shift_grating_GUI is made visible.
function Phase_shift_grating_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Phase_shift_grating_GUI (see VARARGIN)

% Choose default command line output for Phase_shift_grating_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Phase_shift_grating_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Phase_shift_grating_GUI_OutputFcn(~, ~, handles) 
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

phi=str2num(get(handles.edit1,'string'));
T=str2num(get(handles.edit2,'string'));

n_eff=1.46; %Effective refractive index
lamda_B=1550*1e-9;  %Center wavelength
L=0.1;  %Length of the fiber grating
delta_n_eff=2*1e-5; %Refractive rate modulation depth
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

plot(handles.axes1,lamda*1e9,R);
grid on;
xlabel('Wavelength/nm');
ylabel('Reflection');
title('Reflection Spectrum of Phase-Shift Grating');
