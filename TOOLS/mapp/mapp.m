function varargout = mapp(varargin)
% MAPP MATLAB code for mapp.fig
%      MAPP, by itself, creates a new MAPP or raises the existing
%      singleton*.
%
%      H = MAPP returns the handle to a new MAPP or the handle to
%      the existing singleton*.
%
%      MAPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPP.M with the given input arguments.
%
%      MAPP('Property','Value',...) creates a new MAPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mapp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mapp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mapp

% Last Modified by GUIDE v2.5 12-Dec-2014 09:59:34

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mapp_OpeningFcn, ...
    'gui_OutputFcn',  @mapp_OutputFcn, ...
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
end

% --- Executes just before mapp is made visible.
function mapp_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mapp (see VARARGIN)
global strct_config
strct_config.workpath= mfilename('fullpath');
strct_config.workpath=strct_config.workpath(1:end-4);

load_config(handles);

% Choose default command line output for mapp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mapp wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = mapp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line outpendut from handles structure

varargout{1} = handles.output;
end


function ptspa_s_Callback(hObject, eventdata, handles)
% hObject    handle to ptspa_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptspa_s as text
%        str2double(get(hObject,'String')) returns contents of ptspa_s as a double
end

% --- Executes during object creation, after setting all properties.
function ptspa_s_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptspa_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptspa_w_Callback(hObject, ~, handles)
% hObject    handle to ptspa_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptspa_w as text
%        str2double(get(hObject,'String')) returns contents of ptspa_w as a double
end

% --- Executes during object creation, after setting all properties.
function ptspa_w_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptspa_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptspa_redest_Callback(hObject, eventdata, handles)
% hObject    handle to ptspa_redest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptspa_redest as text
%        str2double(get(hObject,'String')) returns contents of ptspa_redest as a double
end

% --- Executes during object creation, after setting all properties.
function ptspa_redest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptspa_redest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptstr_s_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptstr_s as text
%        str2double(get(hObject,'String')) returns contents of ptstr_s as a double
end

% --- Executes during object creation, after setting all properties.
function ptstr_s_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptstr_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptstr_w_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptstr_w as text
%        str2double(get(hObject,'String')) returns contents of ptstr_w as a double
end

% --- Executes during object creation, after setting all properties.
function ptstr_w_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptstr_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptstr_redest_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_redest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptstr_redest as text
%        str2double(get(hObject,'String')) returns contents of ptstr_redest as a double
end

% --- Executes during object creation, after setting all properties.
function ptstr_redest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptstr_redest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptstr_redeen_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_redeen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptstr_redeen as text
%        str2double(get(hObject,'String')) returns contents of ptstr_redeen as a double
end

% --- Executes during object creation, after setting all properties.
function ptstr_redeen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptstr_redeen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptspa_zonename_Callback(hObject, eventdata, handles)
% hObject    handle to ptspa_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptspa_zonename as text
%        str2double(get(hObject,'String')) returns contents of ptspa_zonename as a double
end

% --- Executes during object creation, after setting all properties.
function ptspa_zonename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptspa_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptspa_ma_Callback(hObject, eventdata, handles)
% hObject    handle to ptspa_ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptspa_ma as text
%        str2double(get(hObject,'String')) returns contents of ptspa_ma as a double
end

% --- Executes during object creation, after setting all properties.
function ptspa_ma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptspa_ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptspa_re_Callback(hObject, eventdata, handles)
% hObject    handle to ptspa_re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptspa_re as text
%        str2double(get(hObject,'String')) returns contents of ptspa_re as a double
end

% --- Executes during object creation, after setting all properties.
function ptspa_re_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptspa_re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptstr_zonename_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptstr_zonename as text
%        str2double(get(hObject,'String')) returns contents of ptstr_zonename as a double
end

% --- Executes during object creation, after setting all properties.
function ptstr_zonename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptstr_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptstr_ma_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptstr_ma as text
%        str2double(get(hObject,'String')) returns contents of ptstr_ma as a double
end

% --- Executes during object creation, after setting all properties.
function ptstr_ma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptstr_ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptstr_re_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptstr_re as text
%        str2double(get(hObject,'String')) returns contents of ptstr_re as a double
end

% --- Executes during object creation, after setting all properties.
function ptstr_re_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptstr_re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptd_zonename_Callback(hObject, eventdata, handles)
% hObject    handle to ptd_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptd_zonename as text
%        str2double(get(hObject,'String')) returns contents of ptd_zonename as a double
end

% --- Executes during object creation, after setting all properties.
function ptd_zonename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptd_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptd_ma_Callback(hObject, eventdata, handles)
% hObject    handle to ptd_ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptd_ma as text
%        str2double(get(hObject,'String')) returns contents of ptd_ma as a double
end

% --- Executes during object creation, after setting all properties.
function ptd_ma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptd_ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptd_re_Callback(hObject, eventdata, handles)
% hObject    handle to ptd_re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptd_re as text
%        str2double(get(hObject,'String')) returns contents of ptd_re as a double
end


% --- Executes during object creation, after setting all properties.
function ptd_re_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptd_re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in ptd_pe.
function ptd_pe_Callback(hObject, eventdata, handles)
% hObject    handle to ptd_pe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ptd_pe
end


function pta_zonename_Callback(hObject, eventdata, handles)
% hObject    handle to pta_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pta_zonename as text
%        str2double(get(hObject,'String')) returns contents of pta_zonename as a double
end

% --- Executes during object creation, after setting all properties.
function pta_zonename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pta_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function pta_ma_Callback(hObject, eventdata, handles)
% hObject    handle to pta_ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pta_ma as text
%        str2double(get(hObject,'String')) returns contents of pta_ma as a double
end

% --- Executes during object creation, after setting all properties.
function pta_ma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pta_ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function pta_re_Callback(hObject, eventdata, handles)
% hObject    handle to pta_re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pta_re as text
%        str2double(get(hObject,'String')) returns contents of pta_re as a double
end

% --- Executes during object creation, after setting all properties.
function pta_re_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pta_re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pta_uplus_flag.
function pta_uplus_flag_Callback(hObject, eventdata, handles)
% hObject    handle to pta_uplus_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if get(hObject, 'Value')==1
%     set(handles.pta_uplus, 'Visible', 'On');
%     set(handles.text122, 'Visible', 'On');
% else
%     set(handles.pta_uplus, 'Visible', 'Off');
%     set(handles.text122, 'Visible', 'Off');
% end
% Hint: get(hObject,'Value') returns toggle state of pta_uplus_flag
end

% --- Executes on button press in pta_vort.
function pta_vort_Callback(hObject, eventdata, handles)
% hObject    handle to pta_vort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pta_vort
end

% --- Executes on button press in pta_gren.
function pta_gren_Callback(hObject, eventdata, handles)
% hObject    handle to pta_gren (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pta_gren
end

% --- Executes on button press in pta_tke.ptspa_ma
function pta_tke_Callback(hObject, eventdata, handles)
% hObject    handle to pta_tke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pta_tke
end

% --- Executes on button press in pta_ac.
function pta_ac_Callback(hObject, eventdata, handles)

end

% --- Executes on button press in pta_sc.
function pta_sc_Callback(hObject, eventdata, handles)
% hObject    handle to pta_sc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pta_sc
end


function mph_zonename_Callback(hObject, eventdata, handles)
% hObject    handle to mph_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mph_zonename as text
%        str2double(get(hObject,'String')) returns contents of mph_zonename as a double
end

% --- Executes during object creation, after setting all properties.
function mph_zonename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mph_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function mvh_zonename_Callback(hObject, eventdata, handles)
% hObject    handle to mvh_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mvh_zonename as text
%        str2double(get(hObject,'String')) returns contents of mvh_zonename as a double
end

% --- Executes during object creation, after setting all properties.
function mvh_zonename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mvh_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in run_mph.
function run_mph_Callback(hObject, eventdata, handles)

global strct_config

logbar(handles.status_text,'Starte MergePressureHistories...');
s=get(handles.mph_samples_opt, 'String');
n=strct_config.PHprobes;
j=get(handles.mph_red_opt, 'String');

markiert=get(handles.cgns_path_text, 'Value');

if get(handles.overall_zonename_flag, 'Value')==1
    text=get(handles.overall_zonename, 'String');
    z=text{markiert(1)};
else
    z=get(handles.mph_zonename, 'String');
end

command=sprintf('cd %s && %sC-TOOLS/bin/MergePressureHistories -n %d -z %s',strct_config.folder,strct_config.workpath,n,z);

if get(handles.mph_red, 'Value') == 1
    command=sprintf('%s -j %s', command, j);
end

if get(handles.mph_samples, 'Value') == 1
    command=sprintf('%s -s %s', command, s);
end

logbar(handles.status_text,command);
[status, cmdout] = system(command);
logbar(handles.status_text,cmdout);

%set(handles.analyze_mph, 'Visible', 'On');

end

% --- Executes on button press in run_mvh.
function run_mvh_Callback(hObject, eventdata, handles)

global strct_config
logbar(handles.status_text,'Starte MergeVelocityHistories...');
s=get(handles.mvh_samples_opt, 'String');
n=strct_config.TKEprobes;
j=get(handles.mvh_red_opt, 'String');
markiert=get(handles.cgns_path_text, 'Value');

if get(handles.overall_zonename_flag, 'Value')==1
    text=get(handles.overall_zonename, 'String');
    z=text{markiert(1)};
else
    z=get(handles.mvh_zonename, 'String');
end

command=sprintf('cd %s && %sC-TOOLS/bin/MergeVelocityHistories -n %d -z %s',strct_config.folder,strct_config.workpath,n,z);
if get(handles.mvh_red, 'Value') == 1
    command=sprintf('%s -j %s', command, j);
end
if get(handles.mvh_samples, 'Value') == 1
    command=sprintf('%s -s %s', command, s);
end

logbar(handles.status_text,command);
[status, cmdout] = system(command);
logbar(handles.status_text,cmdout);

% set(handles.analyze_mvh, 'Visible', 'On');

end


function pta_uplus_Callback(hObject, eventdata, handles)
% hObject    handle to pta_uplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pta_uplus as text
%        str2double(get(hObject,'String')) returns contents of pta_uplus as a double
end

% --- Executes during object creation, after setting all properties.
function pta_uplus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pta_uplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in run_ptspa.
function run_ptspa_Callback(hObject, eventdata, handles)

global strct_config
logbar(handles.status_text,'Starte PostTKESpanwise...');
markiert=get(handles.cgns_path_text, 'Value');

if get(handles.overall_zonename_flag, 'Value')==1
    text=get(handles.overall_zonename, 'String');
    z=text{markiert(1)};
else
    z=get(handles.ptspa_zonename, 'String');
end

w=get(handles.ptspa_w, 'String');
s=get(handles.ptspa_s, 'String');
redest=get(handles.ptspa_redest, 'String');

ma=strct_config.Mach;
re=sprintf('%f',strct_config.Re);

command=sprintf('cd %s && %sC-TOOLS/bin/postTKESpanwiseCascade -m %s -r %s -f %s -g %s -x %s -s %s -w %s -z %s',strct_config.folder,strct_config.workpath,ma,re,strct_config.cgns_path,strct_config.cgns_avrg,redest,s,w,z);
logbar(handles.status_text,command)
[status, cmdout] = system(command);


logbar(handles.status_text,cmdout);

set(handles.analyze_ptspa, 'Visible', 'On');


end

% --- Executes on button press in run_ptstr.
function run_ptstr_Callback(hObject, eventdata, handles)

global strct_config
logbar(handles.status_text,'Starte PostTKEStreamwise...');
markiert=get(handles.cgns_path_text, 'Value');

if get(handles.overall_zonename_flag, 'Value')==1
    text=get(handles.overall_zonename, 'String');
    z=text{markiert(1)};
else
    z=get(handles.ptstr_zonename, 'String');
end

w=get(handles.ptstr_w, 'String');
s=get(handles.ptstr_s, 'String');
redest=get(handles.ptstr_redest, 'String');
redeen=get(handles.ptstr_redeen, 'String');

ma=strct_config.Mach;
re=sprintf('%f',strct_config.Re);

command=sprintf('cd %s && %sC-TOOLS/bin/postTKEStreamwiseCascade -m %s -r %s -f %s -g %s -x %s -s %s -w %s -z %s -y %s',strct_config.folder,strct_config.workpath,ma,re,strct_config.cgns_path,strct_config.cgns_avrg,redest,s,w,z,redeen);
logbar(handles.status_text,command);
[status, cmdout] = system(command);


logbar(handles.status_text,cmdout);

set(handles.analyze_ptstr, 'Visible', 'On');


end

% --- Executes on button press in run_ptd.
function run_ptd_Callback(hObject, eventdata, handles)    

global strct_config
logbar(handles.status_text,'Starte PostTransientData...');
markiert=get(handles.cgns_path_text, 'Value');

if get(handles.overall_zonename_flag, 'Value')==1
    text=get(handles.overall_zonename, 'String');
    z=text{markiert(1)};
else
    z=get(handles.pta_zonename, 'String');
end

pe_flag=get(handles.ptd_pe, 'Value');
conv_flag=get(handles.ptd_conv, 'Value');

ma=strct_config.Mach;
re=sprintf('%f',strct_config.Re);

command=sprintf('cd %s && %sC-TOOLS/bin/postTransientData -f %s -m %s -r %s -z %s',strct_config.folder,strct_config.workpath,strct_config.cgns_path,ma,re,z);
if pe_flag==1
    command=sprintf('%s -p',command);
end
if conv_flag==1
    command=sprintf('%s -c',command);
end

logbar(handles.status_text,command);
[status, cmdout] = system(command);

logbar(handles.status_text,cmdout);

end

% --- Executes on button press in run_pta.
function run_pta_Callback(hObject, eventdata, handles)

global strct_config
logbar(handles.status_text,'Starte PostTimeAverage...');

markiert=get(handles.cgns_path_text, 'Value');

if get(handles.overall_zonename_flag, 'Value')==1
    text=get(handles.overall_zonename, 'String');
    z=text{markiert(1)};
else
    z=get(handles.ptd_zonename, 'String');
end

uplus_flag=get(handles.pta_uplus_flag, 'Value');
gren_flag=get(handles.pta_gren, 'Value');
ac_flag=get(handles.pta_ac, 'Value');
wu_flag=get(handles.pta_wu, 'Value');
fluc_flag=get(handles.pta_fluc, 'Value');


ma=strct_config.Mach;
re=sprintf('%f',strct_config.Re);
aoa=sprintf('%.4f',strct_config.AoA/180*pi());

pta_pos=abs(str2double(get(handles.lst_pos, 'String')));
if get(handles.side_selection, 'Value')==2
    pta_pos=(-1)*pta_pos;
end
command=sprintf('cd %s && %sC-TOOLS/bin/postTimeAverage -f %s -m %s -r %s -z %s',strct_config.folder,strct_config.workpath,strct_config.cgns_avrg,ma,re,z);
if uplus_flag==1
    command=sprintf('%s -u',command);
end
if gren_flag==1
    command=sprintf('%s -g',command);
end
if ac_flag==1
    command=sprintf('%s -a %s',command,aoa);
end
if wu_flag==1
    command=sprintf('%s -w',command);
end
if fluc_flag==1
    command=sprintf('%s -q',command);
end

%Profile exportieren
if get(handles.pta_export_profiles, 'Value') == 1
    command=sprintf('%s -b %s',command,pta_pos);
end

logbar(handles.status_text,command);
[status, cmdout] = system(command);

end


% --- Executes on button press in add_file.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to add_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


function cgns_path_text_Callback(hObject, eventdata, handles)

selected=get(handles.cgns_path_text, 'Value');
path=get(handles.cgns_path_text, 'String');
set_file_names(path{selected(1)}, handles);

end

% --- Executes during object creation, after setting all properties.
function cgns_path_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cgns_path_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function status_text_Callback(hObject, eventdata, handles)
% hObject    handle to status_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of status_text as text
%        str2double(get(hObject,'String')) returns contents of status_text as a double
end


% --- Executes during object creation, after setting all properties.
function status_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

global status_text_c_tools;
status_text_c_tools=hObject;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function add_file_Callback(hObject, eventdata, handles)

global strct_config

strct_config.analyze_index=1;

[FileName, Folder] = uigetfile(sprintf('%s/*.cgns', strct_config.folder), 'Select Data');
filepath = [Folder, FileName];

if  FileName~=0
    text=['Datei hinzugefügt: ', filepath] ;
    logbar(handles.status_text,text);
    
    cgns_text=get(handles.cgns_path_text, 'String');
    cgns_text{numel(cgns_text)+1}=sprintf('%s%s',Folder,FileName);
    
    zonename=get(handles.overall_zonename, 'String');
    zonename{numel(zonename)+1}=strct_config.standart_zonename;
    set(handles.overall_zonename, 'String', zonename);
    
    if numel(cgns_text) > 1
        set(handles.overall_zonename_flag, 'Value', 1);
        set(handles.overall_zonename_flag, 'Enable', 'Off');
        set(handles.overall_zonename, 'Enable', 'On');
        set(handles.mph_zonename, 'Enable', 'Off')
        set(handles.mvh_zonename, 'Enable', 'Off');
        set(handles.ptspa_zonename, 'Enable', 'Off');
        set(handles.ptstr_zonename, 'Enable', 'Off');
        set(handles.pta_zonename, 'Enable', 'Off');
        set(handles.ptd_zonename, 'Enable', 'Off');
        set(handles.overall_zonename, 'Enable', 'On');
    end
    
    
    
    set(handles.cgns_path_text, 'String', cgns_text);
    
    set_file_names(filepath, handles);
end

end

function set_file_names(filepath, handles)

global strct_config

[strct_config.folder,FileName_cgns,ext] = fileparts(filepath);
strct_config.cgns_path=filepath;

FileName_cgns=sprintf('%s%s',FileName_cgns,ext);

avrg=strrep(FileName_cgns, '.cgns', '_VarianzTimeAverage.cgns');

strct_config.cgns_avrg= sprintf('%s/%s',strct_config.folder,avrg);

if exist(sprintf('%s/PressureHistory_Points_All.dat',strct_config.folder),'file')==2
    set(handles.text166, 'Visible', 'On');
else
    set(handles.text166, 'Visible', 'Off');
end
if exist(sprintf('%s/TKEHistory_Points_All.dat',strct_config.folder),'file')==2
    set(handles.text165, 'Visible', 'On');
else
    set(handles.text165, 'Visible', 'Off');
end
if exist(sprintf('%s/TKEStreamwiseCascade.dat',strct_config.folder),'file')==2
    set(handles.text168, 'Visible', 'On');
else
    set(handles.text168, 'Visible', 'Off');
end
if exist(sprintf('%s/TKESpanwiseCascade.dat',strct_config.folder),'file')==2
    set(handles.text167, 'Visible', 'On');
else
    set(handles.text167, 'Visible', 'Off');
end

end

% --- Executes on button press in analyze_mph.
function analyze_mph_Callback(hObject, eventdata, handles)

global strct_config
global strct_probes

strct_config.FlagMappCallFA=1;
strct_config.methode=get(handles.methode_mph, 'Value') ;
strct_config.unit_x = 1;
strct_config.unit_y = 1;

strct_config.RawDataPath= sprintf('%s/PressureHistory_Points_All.dat',strct_config.folder);
strct_config.RawDataFileName= sprintf('PressureHistory_Points_All.dat');

%strct_config.yPlusSoll=get(handles.yplus_ptstr, 'String');
strct_config.hanning = get(handles.hanning_flag_mph, 'Value'); %hanning window flag
strct_config.avrg = get(handles.average_flag_mph, 'Value'); %mittelwert flag

strct_config.firstSample= 1; %erstes Sample
strct_config.ar= 30; %AR-Ordnung f�r pBurg

if strct_config.methode==3
    strct_config.pwelch_windowsize= str2double(get(handles.mph_windowsize, 'String')) ;
end
if get(handles.mph_highpass, 'Value') == 1
    strct_config.co_frequency= str2double(get(handles.mph_co_frequency, 'String')) ;
else
    strct_config.co_frequency=0;
end


logbar(handles.status_text,'Starte Analyse und Export von PressureHistory_Points_all.dat ...');

importData('');
strct_config.lastSample= strct_probes(1).samples; %letztes Sample setzen
start_selection();
strct_config.updated=1;
frequenzanalyse('export_tecplot_Callback',0,0,'');

end

% --- Executes on button press in analyze_mvh.
function analyze_mvh_Callback(hObject, eventdata, handles)

global strct_config
global strct_probes

strct_config.FlagMappCallFA=1;
strct_config.unit_x = 1;
strct_config.unit_y = 2;
strct_config.RawDataPath= sprintf('%s/TKEHistory_Points_All.dat',strct_config.folder);
strct_config.RawDataFileName= sprintf('TKEHistory_Points_All.dat');

strct_config.methode=get(handles.methode_mvh, 'Value') ;
%strct_config.yPlusSoll=get(handles.yplus_ptspa, 'String');
strct_config.hanning = get(handles.hanning_flag_mvh, 'Value'); %hanning window flag
strct_config.avrg = get(handles.average_flag_mvh, 'Value'); %mittelwert flag

strct_config.firstSample= 1; %erstes Sample
strct_config.ar= 30; %AR-Ordnung f�r pBurg

if strct_config.methode==3
    strct_config.pwelch_windowsize= str2double(get(handles.mvh_windowsize, 'String')) ;
end
if get(handles.mvh_highpass, 'Value') == 1
    strct_config.co_frequency= str2double(get(handles.mvh_co_frequency, 'String')) ;
else
    strct_config.co_frequency=0;
end

logbar(handles.status_text,'Starte Analyse und Export von TKEHistory_Points_all.dat ...');

importData('');
strct_config.lastSample= strct_probes(1).samples; %letztes Sample setzen
start_selection();
strct_config.updated=1;
frequenzanalyse('export_tecplot_Callback',0,0,'');

end

% --- Executes on button press in analyze_ptspa.
function analyze_ptspa_Callback(hObject, eventdata, handles)

global strct_config strct_probes

strct_config.FlagMappCallFA=1;
strct_config.unit_x = 3;
strct_config.unit_y = 2;
strct_config.methode=get(handles.methode_ptspa, 'Value') ;
strct_config.yPlusSoll=get(handles.yplus_ptspa, 'String');
strct_config.hanning = get(handles.hanning_flag_ptspa, 'Value'); %hanning window flag
strct_config.avrg = get(handles.average_flag_ptspa, 'Value'); %mittelwert flag
strct_config.firstSample= 1; %erstes Sample
strct_config.ar= 30; %AR-Ordnung f�r pBurg
s_end=str2double(get(handles.ptspa_s, 'String'));
if strct_config.methode==3
    strct_config.pwelch_windowsize= str2double(get(handles.ptspa_windowsize, 'String')) ;
end
if get(handles.ptspa_highpass, 'Value') == 1
    strct_config.co_frequency= str2double(get(handles.ptspa_co_frequency, 'String')) ;
else
    strct_config.co_frequency=0;
end


%Schleife über 'Anzahl Samples'
logbar(handles.status_text, sprintf('Starte Analyse und Export von TKESpanwiseCascade_1-%d.dat...', s_end));
for i=1:s_end
    strct_config.RawDataPath= sprintf('%s/TKESpanwiseCascade_%d.dat',strct_config.folder, i);
    strct_config.RawDataFileName= sprintf('TKESpanwiseCascade_%d.dat',i);
    
    importData('');
    strct_config.lastSample= strct_probes(1).samples; %letztes Sample setzen
    start_selection();
    
    %Schleife über die ausgewählten Proben
    for probe=1:length(strct_config.markiert)
        if probe == 1 %erste ausgewählte Probe wird immer addiert
            process_signal(strct_config.markiert(probe));
            statistics('',strct_config.markiert(probe));
            if i==1 %erste Datei
                strct_probes_avrg(strct_config.markiert(probe)).psdY=(1/s_end).*strct_probes(strct_config.markiert(probe)).psdY;
                strct_probes_avrg(strct_config.markiert(probe)).psdX=(1/s_end).*strct_probes(strct_config.markiert(probe)).psdX;
            else    %werte durch s_end teilen und addieren
                strct_probes_avrg(strct_config.markiert(probe)).psdY=strct_probes_avrg(strct_config.markiert(probe)).psdY+(1/s_end).*strct_probes(strct_config.markiert(probe)).psdY;
                strct_probes_avrg(strct_config.markiert(probe)).psdX=strct_probes_avrg(strct_config.markiert(probe)).psdX+(1/s_end).*strct_probes(strct_config.markiert(probe)).psdX;
            end
        else if strct_config.markiert(probe) ~= strct_config.markiert(probe-1) % Prüfen, ob die Datei bereits addiert wurde;
            process_signal(strct_config.markiert(probe));
            statistics('',strct_config.markiert(probe));
            if i==1 %erste Datei
                strct_probes_avrg(strct_config.markiert(probe)).psdY=(1/s_end).*strct_probes(strct_config.markiert(probe)).psdY;
                strct_probes_avrg(strct_config.markiert(probe)).psdX=(1/s_end).*strct_probes(strct_config.markiert(probe)).psdX;
            else    %werte durch s_end teilen und addieren
                strct_probes_avrg(strct_config.markiert(probe)).psdY=strct_probes_avrg(strct_config.markiert(probe)).psdY+(1/s_end).*strct_probes(strct_config.markiert(probe)).psdY;
                strct_probes_avrg(strct_config.markiert(probe)).psdX=strct_probes_avrg(strct_config.markiert(probe)).psdX+(1/s_end).*strct_probes(strct_config.markiert(probe)).psdX;
            end
            end
        end
    end
    % Restliche Einträge mit Nullen auffüllen, damit die Dimension stimmt
    for probe=strct_config.markiert(end)+1:strct_config.probes+1
        strct_probes_avrg(probe).psdY=[];
        strct_probes_avrg(probe).psdX=[];
    end
    
end
%Werte übertragen und Export starten
for i=1:strct_config.probes+1
    strct_probes(i).psdX=strct_probes_avrg(i).psdX;
    strct_probes(i).psdY=strct_probes_avrg(i).psdY;
end

strct_config.RawDataPath= sprintf('%s/TKESpanwiseCascade.dat',strct_config.folder);
strct_config.RawDataFileName= sprintf('TKESpanwiseCascade.dat');

start_selection();
strct_config.updated=0;
frequenzanalyse('export_tecplot_Callback',0,0,'');

end

% --- Executes on button press in analyze_ptstr.
function analyze_ptstr_Callback(hObject, eventdata, handles)

global strct_config strct_probes

strct_config.FlagMappCallFA=1;
strct_config.unit_x = 2;
strct_config.unit_y = 2;
strct_config.methode=get(handles.methode_ptstr, 'Value') ;
strct_config.yPlusSoll=get(handles.yplus_ptstr, 'String');
strct_config.hanning = get(handles.hanning_flag_ptstr, 'Value'); %hanning window flag
strct_config.avrg = get(handles.average_flag_ptstr, 'Value'); %mittelwert flag
strct_config.firstSample= 1; %erstes Sample
strct_config.ar= 30; %AR-Ordnung f�r pBurg
s_end=str2double(get(handles.ptstr_s, 'String'));
if strct_config.methode==3
    strct_config.pwelch_windowsize= str2double(get(handles.ptstr_windowsize, 'String')) ;
end
if get(handles.ptstr_highpass, 'Value') == 1
    strct_config.co_frequency= str2double(get(handles.ptstr_co_frequency, 'String')) ;
else
    strct_config.co_frequency=0;
end


logbar(handles.status_text, sprintf('Starte Analyse und Export von TKEStreamwiseCascade_1-%d.dat...', s_end));
for i=1:s_end
    strct_config.RawDataPath= sprintf('%s/TKEStreamwiseCascade_%d.dat',strct_config.folder, i);
    strct_config.RawDataFileName= sprintf('TKEStreamwiseCascade_%d.dat',i);
    
    importData('');
    strct_config.lastSample= strct_probes(1).samples; %letztes Sample setzen
    start_selection();
    
    %Schleife über die ausgewählten Proben
    for probe=1:length(strct_config.markiert)
        if probe==1 % erste wird immer addiert
            process_signal(strct_config.markiert(probe));
            statistics('',strct_config.markiert(probe));
            if i==1 %erste Datei
                strct_probes_avrg(strct_config.markiert(probe)).psdY=(1/s_end).*strct_probes(strct_config.markiert(probe)).psdY;
                strct_probes_avrg(strct_config.markiert(probe)).psdX=(1/s_end).*strct_probes(strct_config.markiert(probe)).psdX;
            else %werte durch s_end teilen und addieren
                strct_probes_avrg(strct_config.markiert(probe)).psdY=strct_probes_avrg(strct_config.markiert(probe)).psdY+(1/s_end).*strct_probes(strct_config.markiert(probe)).psdY;
                strct_probes_avrg(strct_config.markiert(probe)).psdX=strct_probes_avrg(strct_config.markiert(probe)).psdX+(1/s_end).*strct_probes(strct_config.markiert(probe)).psdX;
            end
        else if strct_config.markiert(probe) ~= strct_config.markiert(probe-1) %Prüfen, ob die Datei bereits addiert wurde; 
            process_signal(strct_config.markiert(probe));
            statistics('',strct_config.markiert(probe));
            if i==1 %erste Datei
                strct_probes_avrg(strct_config.markiert(probe)).psdY=(1/s_end).*strct_probes(strct_config.markiert(probe)).psdY;
                strct_probes_avrg(strct_config.markiert(probe)).psdX=(1/s_end).*strct_probes(strct_config.markiert(probe)).psdX;
            else %werte durch s_end teilen und addieren
                strct_probes_avrg(strct_config.markiert(probe)).psdY=strct_probes_avrg(strct_config.markiert(probe)).psdY+(1/s_end).*strct_probes(strct_config.markiert(probe)).psdY;
                strct_probes_avrg(strct_config.markiert(probe)).psdX=strct_probes_avrg(strct_config.markiert(probe)).psdX+(1/s_end).*strct_probes(strct_config.markiert(probe)).psdX;
            end                
            end
        end
    end
    % Restliche Einträge mit Nullen auffüllen, damit die Dimension stimmt
    for probe=strct_config.markiert(end)+1:strct_config.probes+1
        strct_probes_avrg(probe).psdY=[];
        strct_probes_avrg(probe).psdX=[];
    end
    
end
%Werte übertragen und Export starten
for i=1:strct_config.probes+1
    strct_probes(i).psdX=strct_probes_avrg(i).psdX;
    strct_probes(i).psdY=strct_probes_avrg(i).psdY;
end

strct_config.RawDataPath= sprintf('%s/TKEStreamwiseCascade.dat',strct_config.folder);
strct_config.RawDataFileName= sprintf('TKEStreamwiseCascade.dat');

start_selection();
strct_config.updated=0;
frequenzanalyse('export_tecplot_Callback',0,0,'');

end


function analyze_lst(handles)

global strct_config strct_probes
logbar(handles.status_text,'Starte LST-Tool...');

strct_config.updated=0;
strct_config.markiert=1;
strct_config.methode=7;
strct_config.unit_x=1;

position = abs(str2double(get(handles.lst_pos, 'String')));
if get(handles.side_selection, 'Value')==2
    position=(-1)*position;
end
file=sprintf('%s/BoundaryLayer_%.2f.dat',strct_config.folder,position);

[strct_probes(1).psdX,strct_probes(1).psdY]=LST_NACA(file);

strct_config.RawDataFileName=sprintf('BoundaryLayer_%.2f.dat',position);

frequenzanalyse('export_tecplot_Callback',0,0,'');

end



% --- Executes on button press in auto_flag_mph.
function auto_flag_mph_Callback(hObject, eventdata, handles)
% hObject    handle to auto_flag_mph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_flag_mph
end

% --- Executes on button press in auto_flag_mvh.
function auto_flag_mvh_Callback(hObject, eventdata, handles)
% hObject    handle to auto_flag_mvh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_flag_mvh
end

% --- Executes on button press in auto_flag_ptspa.
function auto_flag_ptspa_Callback(hObject, eventdata, handles)
% hObject    handle to auto_flag_ptspa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_flag_ptspa
end

% --- Executes on button press in auto_flag_ptstr.
function auto_flag_ptstr_Callback(hObject, eventdata, handles)
% hObject    handle to auto_flag_ptstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_flag_ptstr
end

% --- Executes on button press in auto_flag_ptd.
function auto_flag_ptd_Callback(hObject, eventdata, handles)
% hObject    handle to auto_flag_ptd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDTKESpanwiseCascade.datATA)

% Hint: get(hObject,'Value') returns toggle state of auto_flag_ptd
end

% --- Executes on button press in rawdata_flag_pta.
function rawdata_flag_pta_Callback(hObject, eventdata, handles)
% hObject    handle to rawdata_flag_pta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rawdata_flag_pta
end

% --- Executes on button press in auto_run.
function auto_run_Callback(hObject, eventdata, handles)
global strct_config
%führt die ausgewählten C-Tools aus; anschließend wird die erstellte Datei
%mit Frequenzanaylse.m analysiert und exportiert

logbar(handles.status_text,'Starte automatische Auswertung...');


cgns_feld=get(handles.cgns_path_text, 'String');
markiert=get(handles.cgns_path_text, 'Value');
n=1;

for i=1:numel(cgns_feld)
    
    if length(markiert) >= n
        
        if markiert(n)==i
            
            n=n+1;
            set(handles.cgns_path_text, 'Value', i);
            set_file_names(cgns_feld{i}, handles);
            
            if exist(sprintf('%s/config.ini',strct_config.folder),'file')==2
                
                strct_config.ConfigPath= sprintf('%s/config.ini',strct_config.folder);
                fprintf('Lade Datei: %s\n',char(strct_config.ConfigPath));
                
                ini = IniConfig();
                ini.ReadFile(strct_config.ConfigPath);
                strct_config.Re = (ini.GetValues('[FluidProperties]','Reynolds'));
                strct_config.Mach = (ini.GetValues('[FluidProperties]','Mach'));
                strct_config.AoA = (ini.GetValues('[BoundaryConditions]','AoA'));
                strct_config.L_ref_config = (ini.GetValues('[FluidProperties]','L'));
                strct_config.T = (ini.GetValues('[FluidProperties]','T'));
                strct_config.U_ref_config=strct_config.Mach*sqrt(1.4*287*strct_config.T);
                strct_config.TKEprobes=int32((ini.GetValues('[VelocityHistory]','NumberLocations')));
                strct_config.PHprobes=int32((ini.GetValues('[PressureHistory]','NumberLocations')));
                
                fprintf('Importierte Konfiguration: Re=%g, Mach=%g, AoA=%g, L=%g, T=%g, U_ref=%g, NumberLocationsPressure: %g, NumberLocationsTKE: %g \n',...
                    strct_config.Re,strct_config.Mach,strct_config.AoA,...
                    strct_config.L_ref_config,strct_config.T,strct_config.U_ref_config,strct_config.PHprobes,strct_config.TKEprobes );
            else
                logbar(handles.status_text, sprintf('Warnung: config.ini in %s nicht gefunden',strct_config.folder));
            end
            
            if get(handles.auto_flag_ptd, 'Value') == 1
                run_ptd_Callback(0,0,handles);
            end
            if get(handles.rawdata_flag_pta, 'Value') == 1 
                run_pta_Callback(0,0,handles);
            end
            if get(handles.lst_analyze_flag, 'Value') ==1
                analyze_lst(handles);
            end
            if get(handles.raw_data_mph, 'Value') == 1
                run_mph_Callback(0,0,handles);
            end
            if get(handles.auto_flag_mph, 'Value') == 1
                analyze_mph_Callback(0,0,handles);
            end
            if get(handles.raw_data_mvh, 'Value') == 1
                run_mvh_Callback(0,0,handles);
            end
            if get(handles.auto_flag_mvh, 'Value') == 1
                analyze_mvh_Callback(0,0,handles);
            end
            if get(handles.raw_data_ptspa, 'Value') == 1
                run_ptspa_Callback(0,0,handles);
            end
            if get(handles.auto_flag_ptspa, 'Value') == 1
                analyze_ptspa_Callback(0,0,handles);
            end
            if get(handles.raw_data_ptstr, 'Value') == 1
                run_ptstr_Callback(0,0,handles);
            end
            if get(handles.auto_flag_ptstr, 'Value') == 1
                analyze_ptstr_Callback(0,0,handles);
            end
        end
    end
end

set(handles.cgns_path_text, 'Value', markiert);

logbar(handles.status_text,'automatische Auswertung beendet');
end




function write_config(handles)
%Schreibt die Einstellungen aus der GUI in die Datei 'config'

global strct_config

%%%%%%%%%%%
ini = IniConfig();
ini.ReadFile(sprintf('%smapp.ini',strct_config.workpath));

%Standartordner und -zonename
cgns_path= get(handles.cgns_path_text, 'String');
zonename = get(handles.overall_zonename, 'String');

cgns_string=cgns_path{1};
zonename_string=zonename{1};

for i = 2:length(cgns_path)
    cgns_string=sprintf('%s,%s',cgns_string,cgns_path{i});
    zonename_string=sprintf('%s,%s',zonename_string,zonename{i});
end

ini.SetValues('[Allgemein]','folder', cgns_string);
ini.SetValues('[Allgemein]','zonename', zonename_string);

ini.SetValues('[Allgemein]','standart_folder' ,strct_config.folder);
ini.SetValues('[Allgemein]','standart_zonename', zonename{1});

ini.SetValues('[Allgemein]','global_L_ref' ,get(handles.L_ref_text, 'String'));
ini.SetValues('[Allgemein]','L_ref_flag', get(handles.L_ref_flag, 'Value'));

%Einstellungen C-Tools
ini.SetValues('[Einstellungen C-Tools]','MPH_zonename', get(handles.mph_zonename, 'String'));
ini.SetValues('[Einstellungen C-Tools]','MPH_red_flag', get(handles.mph_red, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','MPH_red_opt', get(handles.mph_red_opt, 'String'));
ini.SetValues('[Einstellungen C-Tools]','MPH_samples_flag', get(handles.mph_samples, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','MPH_samples_opt', get(handles.mph_samples_opt, 'String'));
ini.SetValues('[Einstellungen C-Tools]','MVH_zonename', get(handles.mvh_zonename, 'String'));
ini.SetValues('[Einstellungen C-Tools]','MVH_red_flag', get(handles.mvh_red, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','MVH_red_opt', get(handles.mvh_red_opt, 'String'));
ini.SetValues('[Einstellungen C-Tools]','MVH_samples_flag', get(handles.mvh_samples, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','MVH_samples_opt', get(handles.mvh_samples_opt, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTSPA_zonename', get(handles.ptspa_zonename, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTSPA_Schrittweite', get(handles.ptspa_w, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTSPA_samples', get(handles.ptspa_s, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTSPA_RedeStart', get(handles.ptspa_redest, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTSTR_zonename', get(handles.ptstr_zonename, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTSTR_Schrittweite', get(handles.ptstr_w, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTSTR_samples', get(handles.ptstr_s, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTSTR_RedeStart', get(handles.ptstr_redest, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTSTR_RedeEnde', get(handles.ptstr_redeen, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTD_zonename', get(handles.ptd_zonename, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTD_Pressure_Extr_flag', get(handles.ptd_pe, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','PTD_Convergence_flag', get(handles.ptd_conv, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','PTA_zonename', get(handles.pta_zonename, 'String'));
ini.SetValues('[Einstellungen C-Tools]','PTA_uplus_flag', get(handles.pta_uplus_flag, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','PTA_grenzschicht_flag', get(handles.pta_gren, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','aerodynamic_coefficient_flag', get(handles.pta_ac, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','PTA_wall_units_flag', get(handles.pta_wu, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','PTA_fluctuation_flag', get(handles.pta_fluc, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','raw_data_mvh_flag', get(handles.raw_data_mvh, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','raw_data_mph_flag', get(handles.raw_data_mph, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','raw_data_ptspa_flag', get(handles.raw_data_ptspa, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','raw_data_ptstr_flag', get(handles.raw_data_ptstr, 'Value'));
ini.SetValues('[Einstellungen C-Tools]','LST_Analyse_Flag', get(handles.lst_analyze_flag,'Value'  ));
ini.SetValues('[Einstellungen C-Tools]','PTA_extract_profile_flag',get(handles.pta_export_profiles,'Value'));
ini.SetValues('[Einstellungen C-Tools]','PTA_profile_pos',get(handles.lst_pos,'String' ));
ini.SetValues('[Einstellungen C-Tools]','PTA_profile_side',get(handles.side_selection,'Value' ));

%Einstellungen Frequenzanalyse
ini.SetValues('[Einstellungen Frequenzanalyse]','MPH_Methode_flag', get(handles.methode_mph, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MPH_Hanning_flag', get(handles.hanning_flag_mph, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MPH_Average_flag', get(handles.average_flag_mph, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MPH_pwelch_windowsize', get(handles.mph_windowsize, 'String'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MPH_highpass_flag', get(handles.mph_highpass, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MPH_co_frequency', get(handles.mph_co_frequency, 'String'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MVH_Methode_flag', get(handles.methode_mvh, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MVH_Hanning_flag', get(handles.hanning_flag_mvh, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MVH_Average_flag', get(handles.average_flag_mvh, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MVH_pwelch_windowsize', get(handles.mvh_windowsize, 'String'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MVH_highpass_flag', get(handles.mvh_highpass, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','MVH_co_frequency', get(handles.mvh_co_frequency, 'String'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSPA_Methode_flag', get(handles.methode_ptspa, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSPA_Hanning_flag', get(handles.hanning_flag_ptspa, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSPA_Average_flag', get(handles.average_flag_ptspa, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSPA_yplus', get(handles.yplus_ptspa, 'String'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSPA_pwelch_windowsize', get(handles.ptspa_windowsize, 'String'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSPA_highpass_flag', get(handles.ptspa_highpass, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSPA_co_frequency', get(handles.ptspa_co_frequency, 'String'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSTR_Methode_flag', get(handles.methode_ptstr, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSTR_Hanning_flag', get(handles.hanning_flag_ptstr, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSTR_Average_flag', get(handles.average_flag_ptstr, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSTR_yplus', get(handles.yplus_ptstr, 'String'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSTR_pwelch_windowsize', get(handles.ptstr_windowsize, 'String'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSTR_highpass_flag', get(handles.ptstr_highpass, 'Value'));
ini.SetValues('[Einstellungen Frequenzanalyse]','PTSTR_co_frequency', get(handles.ptstr_co_frequency, 'String'));

ini.SetValues('[Einstellungen LST]','alpha_start', get(handles.alpha_start, 'String'));
ini.SetValues('[Einstellungen LST]','alpha_end', get(handles.alpha_end, 'String'));
ini.SetValues('[Einstellungen LST]','alpha_steps', get(handles.alpha_steps, 'String'));

ini.WriteFile(sprintf('%smapp.ini',strct_config.workpath));
logbar(handles.status_text,'Einstellungen gespeichert');

end

% --- Executes on button press in save_config.
function save_config_Callback(hObject, eventdata, handles)

write_config(handles);

end

function load_config(handles)

global strct_config


if exist(sprintf('%smapp.ini',strct_config.workpath),'file')==2
   
    ini = IniConfig();
    ini.ReadFile(sprintf('%smapp.ini',strct_config.workpath));
        
    cgns_path=ini.GetValues('[Allgemein]','folder');
    zonename=ini.GetValues('[Allgemein]','zonename');
    
    cgns_path=strsplit(cgns_path, ',');
    zonename=strsplit(zonename, ',');
    
    %Standartordner, Standartzonename, Referenzlänge
    
    set(handles.cgns_path_text, 'String', cgns_path);
    set(handles.overall_zonename, 'String',zonename); 
    
    strct_config.folder=ini.GetValues('[Allgemein]','standart_folder');
    strct_config.standart_zonename = ini.GetValues('[Allgemein]','standart_zonename');
    
    strct_config.global_L_ref = ini.GetValues('[Allgemein]','global_L_ref');
    strct_config.L_ref_flag = ini.GetValues('[Allgemein]','L_ref_flag');
    
    set(handles.L_ref_flag, 'Value', strct_config.L_ref_flag) ;
    set(handles.L_ref_text, 'String', strct_config.global_L_ref) ;
    
    %Werte für C-Tools
    
    set(handles.mph_zonename, 'String',ini.GetValues('[Einstellungen C-Tools]','MPH_zonename') );
    set(handles.mph_red, 'Value',ini.GetValues('[Einstellungen C-Tools]','MPH_red_flag') );
    set(handles.mph_red_opt,'String', ini.GetValues('[Einstellungen C-Tools]','MPH_red_opt') );
    set(handles.mph_samples, 'Value',ini.GetValues('[Einstellungen C-Tools]','MPH_samples_flag') );
    set(handles.mph_samples_opt,'String', ini.GetValues('[Einstellungen C-Tools]','MPH_samples_opt') );    
    set(handles.mvh_zonename,'String', ini.GetValues('[Einstellungen C-Tools]','MVH_zonename') );
    set(handles.mvh_red, 'Value',ini.GetValues('[Einstellungen C-Tools]','MVH_red_flag') );
    set(handles.mvh_red_opt,'String', ini.GetValues('[Einstellungen C-Tools]','MVH_red_opt') );
    set(handles.mvh_samples, 'Value',ini.GetValues('[Einstellungen C-Tools]','MVH_samples_flag') );
    set(handles.mvh_samples_opt,'String', ini.GetValues('[Einstellungen C-Tools]','MVH_samples_opt') );    
    set(handles.ptspa_zonename,'String', ini.GetValues('[Einstellungen C-Tools]','PTSPA_zonename' ));
    set(handles.ptspa_w, 'String',ini.GetValues('[Einstellungen C-Tools]','PTSPA_Schrittweite') );
    set(handles.ptspa_s,'String', ini.GetValues('[Einstellungen C-Tools]','PTSPA_samples') );
    set(handles.ptspa_redest,'String', ini.GetValues('[Einstellungen C-Tools]','PTSPA_RedeStart') );
    set(handles.ptstr_zonename,'String', ini.GetValues('[Einstellungen C-Tools]','PTSTR_zonename') );
    set(handles.ptstr_w, 'String',ini.GetValues('[Einstellungen C-Tools]','PTSTR_Schrittweite') );
    set(handles.ptstr_s, 'String',ini.GetValues('[Einstellungen C-Tools]','PTSTR_samples') );
    set(handles.ptstr_redest, 'String',ini.GetValues('[Einstellungen C-Tools]','PTSTR_RedeStart', '%*s %s') );
    set(handles.ptstr_redeen, 'String',ini.GetValues('[Einstellungen C-Tools]','PTSTR_RedeEnde' ));
    set(handles.ptd_zonename, 'String',ini.GetValues('[Einstellungen C-Tools]','PTD_zonename') );
    set(handles.ptd_pe,'Value', ini.GetValues('[Einstellungen C-Tools]','PTD_Pressure_Extr_flag') );
    set(handles.ptd_conv,'Value', ini.GetValues('[Einstellungen C-Tools]','PTD_Convergence_flag') );
    set(handles.pta_zonename,'String', ini.GetValues('[Einstellungen C-Tools]','PTA_zonename') );
    set(handles.pta_uplus_flag, 'Value',ini.GetValues('[Einstellungen C-Tools]','PTA_uplus_flag') );
    set(handles.pta_gren,'Value', ini.GetValues('[Einstellungen C-Tools]','PTA_grenzschicht_flag') );
    set(handles.pta_ac,'Value', ini.GetValues('[Einstellungen C-Tools]','aerodynamic_coefficient_flag') );
    set(handles.pta_wu, 'Value',ini.GetValues('[Einstellungen C-Tools]','PTA_wall_units_flag') );
    set(handles.lst_analyze_flag,'Value', ini.GetValues('[Einstellungen C-Tools]','LST_Analyse_Flag') );
    set(handles.pta_export_profiles,'Value', ini.GetValues('[Einstellungen C-Tools]','PTA_extract_profile_flag') );
    set(handles.lst_pos,'String', ini.GetValues('[Einstellungen C-Tools]','PTA_profile_pos') );
    set(handles.side_selection,'Value', ini.GetValues('[Einstellungen C-Tools]','PTA_profile_side') );

    set(handles.raw_data_mvh, 'Value',ini.GetValues('[Einstellungen C-Tools]','raw_data_mvh_flag') );
    set(handles.raw_data_mph, 'Value',ini.GetValues('[Einstellungen C-Tools]','raw_data_mph_flag') );
    set(handles.raw_data_ptspa, 'Value',ini.GetValues('[Einstellungen C-Tools]','raw_data_ptspa_flag') );
    set(handles.raw_data_ptstr, 'Value',ini.GetValues('[Einstellungen C-Tools]','raw_data_ptstr_flag') );
    
    %Werte für Frequenzanalyse
    set(handles.methode_mph, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','MPH_Methode_flag') );
    set(handles.hanning_flag_mph, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','MPH_Hanning_flag') );
    set(handles.average_flag_mph, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','MPH_Average_flag') );
    set(handles.mph_windowsize, 'String',ini.GetValues('[Einstellungen Frequenzanalyse]','MPH_pwelch_windowsize') );
    set(handles.mph_highpass, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','MPH_highpass_flag') );
    set(handles.mph_co_frequency, 'String',ini.GetValues('[Einstellungen Frequenzanalyse]','MPH_co_frequency') );
    set(handles.methode_mvh, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','MVH_Methode_flag') );
    set(handles.hanning_flag_mvh, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','MVH_Hanning_flag') );
    set(handles.average_flag_mvh, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','MVH_Average_flag') );
    set(handles.mvh_windowsize, 'String',ini.GetValues('[Einstellungen Frequenzanalyse]','MVH_pwelch_windowsize') );
    set(handles.mvh_highpass, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','MVH_highpass_flag') );
    set(handles.mvh_co_frequency, 'String',ini.GetValues('[Einstellungen Frequenzanalyse]','MVH_co_frequency') );
    set(handles.methode_ptspa, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSPA_Methode_flag') );
    set(handles.hanning_flag_ptspa, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSPA_Hanning_flag') );
    set(handles.average_flag_ptspa, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSPA_Average_flag') );
    PTSPA_yplus=ini.GetValues('[Einstellungen Frequenzanalyse]','PTSPA_yplus' );
    set(handles.ptspa_windowsize, 'String',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSPA_pwelch_windowsize') );
    set(handles.ptspa_highpass, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSPA_highpass_flag') );
    set(handles.ptspa_co_frequency, 'String',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSPA_co_frequency') );
    set(handles.methode_ptstr, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSTR_Methode_flag') );
    set(handles.hanning_flag_ptstr,'Value',  ini.GetValues('[Einstellungen Frequenzanalyse]','PTSTR_Hanning_flag') );
    set(handles.average_flag_ptstr, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSTR_Average_flag') );
    PTSTR_yplus=ini.GetValues('[Einstellungen Frequenzanalyse]','PTSTR_yplus') ;
    set(handles.ptstr_windowsize, 'String',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSTR_pwelch_windowsize') );
    set(handles.ptstr_highpass, 'Value',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSTR_highpass_flag') );
    set(handles.ptstr_co_frequency, 'String',ini.GetValues('[Einstellungen Frequenzanalyse]','PTSTR_co_frequency') );
    
    strct_config.alpha_start = ini.GetValues('[Einstellungen LST]','alpha_start') ;
    strct_config.alpha_end = ini.GetValues('[Einstellungen LST]','alpha_end');
    strct_config.alpha_steps = ini.GetValues('[Einstellungen LST]','alpha_steps') ;
    
    set(handles.alpha_start, 'String', strct_config.alpha_start );
    set(handles.alpha_end, 'String', strct_config.alpha_end );
    set(handles.alpha_steps, 'String', strct_config.alpha_steps );
    
    
    PTSTR_yplus_string = sprintf('%d',PTSTR_yplus(1));
    PTSPA_yplus_string = sprintf('%d',PTSPA_yplus(1));
    
    for i=2:length(PTSPA_yplus)
        PTSPA_yplus_string=sprintf('%s %d',PTSPA_yplus_string,PTSPA_yplus(i));
    end
    for i=2:length(PTSTR_yplus)
        PTSTR_yplus_string=sprintf('%s %d',PTSTR_yplus_string,PTSTR_yplus(i));
    end
    
     set(handles.yplus_ptspa, 'String',PTSPA_yplus_string);
     set(handles.yplus_ptstr, 'String',PTSTR_yplus_string);
        
        
    if get(handles.mph_red, 'Value')==1
        set(handles.mph_red_opt, 'Visible', 'On');
    else
        set(handles.mph_red_opt, 'Visible', 'Off');
    end
    
    if get(handles.mvh_red, 'Value')==1
        set(handles.mvh_red_opt, 'Visible', 'On');
    else
        set(handles.mvh_red_opt, 'Visible', 'Off');
    end
    
    if get(handles.mph_samples, 'Value')==1
        set(handles.mph_samples_opt, 'Visible', 'On');
    else
        set(handles.mph_samples_opt, 'Visible', 'Off');
    end
    
    if get(handles.mvh_samples, 'Value')==1
        set(handles.mvh_samples_opt, 'Visible', 'On');
    else
        set(handles.mvh_samples_opt, 'Visible', 'Off');
    end    
    
    if get(handles.methode_mph, 'Value')==3
        set(handles.mph_windowsize, 'Visible', 'On');
        set(handles.text135, 'Visible', 'On');
    else
        set(handles.mph_windowsize, 'Visible', 'Off');
        set(handles.text135, 'Visible', 'Off');
    end
    if get(handles.methode_mvh, 'Value')==3
        set(handles.mvh_windowsize, 'Visible', 'On');
        set(handles.text136, 'Visible', 'On');
    else
        set(handles.mvh_windowsize, 'Visible', 'Off');
        set(handles.text136, 'Visible', 'Off');
    end
    if get(handles.methode_ptspa, 'Value')==3
        set(handles.ptspa_windowsize, 'Visible', 'On');
        set(handles.text137, 'Visible', 'On');
    else
        set(handles.ptspa_windowsize, 'Visible', 'Off');
        set(handles.text137, 'Visible', 'Off');
    end
    if get(handles.methode_ptstr, 'Value')==3
        set(handles.ptstr_windowsize, 'Visible', 'On');
        set(handles.text138, 'Visible', 'On');
    else
        set(handles.ptstr_windowsize, 'Visible', 'Off');
        set(handles.text138, 'Visible', 'Off');
    end
    if get(handles.mph_highpass, 'Value')==1
        set(handles.mph_co_frequency, 'Visible', 'On');
        set(handles.text131, 'Visible', 'On');
    else
        set(handles.mph_co_frequency, 'Visible', 'Off');
        set(handles.text131, 'Visible', 'Off');
    end
    if get(handles.mvh_highpass, 'Value')==1
        set(handles.mvh_co_frequency, 'Visible', 'On');
        set(handles.text132, 'Visible', 'On');
    else
        set(handles.mvh_co_frequency, 'Visible', 'Off');
        set(handles.text132, 'Visible', 'Off');
    end
    if get(handles.ptspa_highpass, 'Value')==1
        set(handles.ptspa_co_frequency, 'Visible', 'On');
        set(handles.text133, 'Visible', 'On');
    else
        set(handles.ptspa_co_frequency, 'Visible', 'Off');
        set(handles.text133, 'Visible', 'Off');
    end
    if get(handles.ptstr_highpass, 'Value')==1
        set(handles.mph_co_frequency, 'Visible', 'On');
        set(handles.text134, 'Visible', 'On');
    else
        set(handles.ptstr_co_frequency, 'Visible', 'Off');
        set(handles.text134, 'Visible', 'Off');
    end
    
    if length(cgns_path) > 1
        set(handles.overall_zonename_flag, 'Value', 1);
        set(handles.overall_zonename_flag, 'Enable', 'Off');
        set(handles.overall_zonename, 'Enable', 'On');
        set(handles.mph_zonename, 'Enable', 'Off')
        set(handles.mvh_zonename, 'Enable', 'Off');
        set(handles.ptspa_zonename, 'Enable', 'Off');
        set(handles.ptstr_zonename, 'Enable', 'Off');
        set(handles.pta_zonename, 'Enable', 'Off');
        set(handles.ptd_zonename, 'Enable', 'Off');
        set(handles.overall_zonename, 'Enable', 'On');
    end
    
    if get(handles.L_ref_flag, 'Value')==1
        set(handles.L_ref_text, 'Visible', 'On');
        set(handles.text160, 'Visible', 'On');             
    else
        set(handles.L_ref_text, 'Visible', 'Off');
        set(handles.text160, 'Visible', 'Off');  
    end
    
    if get(handles.lst_analyze_flag, 'Value')==1
        set(handles.alpha_start, 'Visible', 'On');
        set(handles.alpha_end, 'Visible', 'On'); 
        set(handles.alpha_steps, 'Visible', 'On');
        set(handles.text161, 'Visible', 'On');
        set(handles.text162, 'Visible', 'On'); 
        set(handles.text163, 'Visible', 'On');
    else
        set(handles.alpha_start, 'Visible', 'Off');
        set(handles.alpha_end, 'Visible', 'Off'); 
        set(handles.alpha_steps, 'Visible', 'Off');
        set(handles.text161, 'Visible', 'Off');
        set(handles.text162, 'Visible', 'Off'); 
        set(handles.text163, 'Visible', 'Off');
    end   
    
    
    set_file_names(cgns_path{1}, handles);
 
else
    logbar(handles.status_text,'Warnung: mapp.ini nicht gefunden');
end

end


% --- Executes on button press in hanning_flag_mph.
function hanning_flag_mph_Callback(hObject, eventdata, handles)
% hObject    handle to hanning_flag_mph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hanning_flag_mph
end

% --- Executes on button press in average_flag_mph.
function average_flag_mph_Callback(hObject, eventdata, handles)
% hObject    handle to average_flag_mph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of average_flag_mph
end


% --- Executes on selection change in methode_mph.
function methode_mph_Callback(hObject, eventdata, handles)
% hObject    handle to methode_mph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methode_mph contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methode_mph

if get(hObject, 'Value') == 3
    set(handles.mph_windowsize, 'Visible', 'On');
    set(handles.text135, 'Visible', 'On');
else
    set(handles.mph_windowsize, 'Visible', 'Off');
    set(handles.text135, 'Visible', 'Off');
end

end

% --- Executes during object creation, after setting all properties.
function methode_mph_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methode_mph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in hanning_flag_mvh.
function hanning_flag_mvh_Callback(hObject, eventdata, handles)
% hObject    handle to hanning_flag_mvh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hanning_flag_mvh
end

% --- Executes on button press in average_flag_mvh.
function average_flag_mvh_Callback(hObject, eventdata, handles)
% hObject    handle to average_flag_mvh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of average_flag_mvh
end

% --- Executes on selection change in methode_mvh.
function methode_mvh_Callback(hObject, eventdata, handles)
% hObject    handle to methode_mvh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methode_mvh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methode_mvh

if get(hObject, 'Value') == 3
    set(handles.mvh_windowsize, 'Visible', 'On');
    set(handles.text136, 'Visible', 'On');
else
    set(handles.mvh_windowsize, 'Visible', 'Off');
    set(handles.text136, 'Visible', 'Off');
end

end

% --- Executes during object creation, after setting all properties.
function methode_mvh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methode_mvh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in hanning_flag_ptspa.
function hanning_flag_ptspa_Callback(hObject, eventdata, handles)
% hObject    handle to hanning_flag_ptspa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hanning_flag_ptspa
end

% --- Executes on button press in average_flag_ptspa.
function average_flag_ptspa_Callback(hObject, eventdata, handles)
% hObject    handle to average_flag_ptspa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of average_flag_ptspa
end

% --- Executes on selection change in methode_ptspa.
function methode_ptspa_Callback(hObject, eventdata, handles)
% hObject    handle to methode_ptspa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methode_ptspa contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methode_ptspa
if get(hObject, 'Value') == 3
    set(handles.ptspa_windowsize, 'Visible', 'On');
    set(handles.text137, 'Visible', 'On');
else
    set(handles.ptspa_windowsize, 'Visible', 'Off');
    set(handles.text137, 'Visible', 'Off');
end

end

% --- Executes during object creation, after setting all properties.
function methode_ptspa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methode_ptspa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in hanning_flag_ptstr.
function hanning_flag_ptstr_Callback(hObject, eventdata, handles)
% hObject    handle to hanning_flag_ptstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hanning_flag_ptstr
end

% --- Executes on button press in average_flag_ptstr.
function average_flag_ptstr_Callback(hObject, eventdata, handles)
% hObject    handle to average_flag_ptstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of average_flag_ptstr
end

% --- Executes on selection change in methode_ptstr.
function methode_ptstr_Callback(hObject, eventdata, handles)
% hObject    handle to methode_ptstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methode_ptstr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methode_ptstr

if get(hObject, 'Value') == 3
    set(handles.ptstr_windowsize, 'Visible', 'On');
    set(handles.text138, 'Visible', 'On');
else
    set(handles.ptstr_windowsize, 'Visible', 'Off');
    set(handles.text138, 'Visible', 'Off');
end

end

% --- Executes during object creation, after setting all properties.
function methode_ptstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methode_ptstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function yplus_ptspa_Callback(hObject, eventdata, handles)
% hObject    handle to yplus_ptspa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yplus_ptspa as text
%        str2double(get(hObject,'String')) returns contents of yplus_ptspa as a double
end

% --- Executes during object creation, after setting all properties.
function yplus_ptspa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yplus_ptspa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function lst_pos_Callback(hObject, eventdata, handles)
% hObject    handle to lst_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lst_pos as text
%        str2double(get(hObject,'String')) returns contents of lst_pos as a double
end

% --- Executes during object creation, after setting all properties.
function lst_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lst_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in overall_zonename.
function overall_zonename_Callback(hObject, eventdata, handles)
% hObject    handle to overall_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)    set(handles.mph_zonename, 'Enable', 'Off');
set(handles.mvh_zonename, 'Enable', 'Off');
set(handles.ptspa_zonename, 'Enable', 'Off');
set(handles.ptstr_zonename, 'Enable', 'Off');
set(handles.pta_zonename, 'Enable', 'Off');
set(handles.ptd_zonename, 'Enable', 'Off');
set(handles.overall_zonename, 'Enable', 'On');

% Hints: contents = cellstr(get(hObject,'String')) returns overall_zonename contents as cell array
%        contents{get(hObject,'Value')} returns selected item from overall_zonename
end

% --- Executes during object creation, after setting all properties.
function overall_zonename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overall_zonename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in overall_zonename_flag.
function overall_zonename_flag_Callback(hObject, eventdata, handles)

flag=get(hObject, 'Value');
if flag==1
    set(handles.mph_zonename, 'Enable', 'Off');
    set(handles.mvh_zonename, 'Enable', 'Off');
    set(handles.ptspa_zonename, 'Enable', 'Off');
    set(handles.ptstr_zonename, 'Enable', 'Off');
    set(handles.pta_zonename, 'Enable', 'Off');
    set(handles.ptd_zonename, 'Enable', 'Off');
    set(handles.overall_zonename, 'Enable', 'On');
else
    set(handles.mph_zonename, 'Enable', 'On');
    set(handles.mvh_zonename, 'Enable', 'On');
    set(handles.ptspa_zonename, 'Enable', 'On');
    set(handles.ptstr_zonename, 'Enable', 'On');
    set(handles.pta_zonename, 'Enable', 'On');
    set(handles.ptd_zonename, 'Enable', 'On');
    set(handles.overall_zonename, 'Enable', 'Off');
end
end

% --- Executes on button press in pta_wu.
function pta_wu_Callback(hObject, eventdata, handles)
% hObject    handle to pta_wu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pta_wu
end


% --- Executes on button press in delete_file.
function delete_file_Callback(hObject, eventdata, handles)

text_cgns = get(handles.cgns_path_text, 'String');
text_zonename = get(handles.overall_zonename, 'String');

markiert = get(handles.cgns_path_text, 'Value');

for i=1:numel(text_cgns)-markiert
    text_cgns{i+markiert-1}=text_cgns{i+markiert};
    text_zonename{i+markiert-1}=text_zonename{i+markiert};
end

text_cgns(numel(text_cgns))=[];
text_zonename(numel(text_zonename))=[];

set(handles.cgns_path_text, 'Value', 1);
set(handles.overall_zonename, 'Value', 1);
set(handles.cgns_path_text, 'String', text_cgns);
set(handles.overall_zonename, 'String', text_zonename);

if numel(text_cgns)>=1
    set_file_names(text_cgns{1}, handles);
end

if numel(text_cgns)<2
    set(handles.overall_zonename_flag, 'Enable', 'On');
    set(handles.overall_zonename_flag, 'Value', 1);
end

end


% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, eventdata, handles)
global processBuilder 
% 
% command=sprintf('pkill -SIGTERM %s',pid);
% 
% [status, cmdout] = system(command);
% logbar(handles.status_text,cmdout);
processBuilder.destroy();
logbar(handles.status_text,'Prozess beendet');

end



function pta_aoa_Callback(hObject, eventdata, handles)
% hObject    handle to pta_aoa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pta_aoa as text
%        str2double(get(hObject,'String')) returns contents of pta_aoa as a double
end

% --- Executes during object creation, after setting all properties.
function pta_aoa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pta_aoa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in skip_ctools_flag.
function skip_ctools_flag_Callback(hObject, eventdata, handles)
% hObject    handle to skip_ctools_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of skip_ctools_flag
end


% --- Executes on button press in raw_data_mph.
function raw_data_mph_Callback(hObject, eventdata, handles)
% hObject    handle to raw_data_mph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw_data_mph
end

% --- Executes on button press in raw_data_mvh.
function raw_data_mvh_Callback(hObject, eventdata, handles)
% hObject    handle to raw_data_mvh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw_data_mvh
end

% --- Executes on button press in raw_data_ptspa.
function raw_data_ptspa_Callback(hObject, eventdata, handles)
% hObject    handle to raw_data_ptspa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw_data_ptspa
end

% --- Executes on button press in raw_data_ptstr.
function raw_data_ptstr_Callback(hObject, eventdata, handles)
% hObject    handle to raw_data_ptstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw_data_ptstr
end


% --- Executes on button press in pta_fluc.
function pta_fluc_Callback(hObject, eventdata, handles)
% hObject    handle to pta_fluc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pta_fluc
end

% --- Executes on button press in ptd_conv.
function ptd_conv_Callback(hObject, eventdata, handles)
% hObject    handle to ptd_conv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ptd_conv
end

function start_selection()
global strct_config 

if strct_config.unit_x==2 || strct_config.unit_x==3
    frequenzanalyse('select_yplus_Callback',0,0,'');
elseif strct_config.unit_x==1 || strct_config.unit_x==4
    for i=1:strct_config.probes
        markiert(i)=i;
    end
    strct_config.markiert=markiert;
end

 
end

function run_command(command, dir)

global processBuilder 

processBuilder = java.lang.ProcessBuilder(command);
processBuilder.directory(java.io.File(dir));
cmdProcess = processBuilder.start();

% Set up a reader to read the output from the command prompt
reader = java.io.BufferedReader(java.io.InputStreamReader(...
    cmdProcess.getInputStream() ));

% Loop until there is some output
nextLine = char( reader.readLine );
while isempty(nextLine)
    nextLine = char( reader.readLine );
end

% Then loop until there is no more output
while ~isempty(nextLine);
    fprintf('bash: %s\n', nextLine);
    nextLine = char( reader.readLine );
end

% Get the exit value of the process
exitValue = cmdProcess.exitValue ;

end


% --- Executes on button press in mph_highpass.
function mph_highpass_Callback(hObject, eventdata, handles)
% hObject    handle to mph_highpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mph_highpass

if get(hObject, 'Value')==1
    set(handles.mph_co_frequency, 'Visible', 'On');
    set(handles.text131, 'Visible', 'On');
else
    set(handles.mph_co_frequency, 'Visible', 'Off');
    set(handles.text131, 'Visible', 'Off');
end

end


function mph_co_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to mph_co_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mph_co_frequency as text
%        str2double(get(hObject,'String')) returns contents of mph_co_frequency as a double
end

% --- Executes during object creation, after setting all properties.
function mph_co_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mph_co_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in mvh_highpass.
function mvh_highpass_Callback(hObject, eventdata, handles)
% hObject    handle to mvh_highpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mvh_highpass
if get(hObject, 'Value')==1
    set(handles.mvh_co_frequency, 'Visible', 'On');
    set(handles.text132, 'Visible', 'On');
else
    set(handles.mvh_co_frequency, 'Visible', 'Off');
    set(handles.text132, 'Visible', 'Off');
end

end


function mvh_co_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to mvh_co_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mvh_co_frequency as text
%        str2double(get(hObject,'String')) returns contents of mvh_co_frequency as a double
end

% --- Executes during object creation, after setting all properties.
function mvh_co_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mvh_co_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in ptspa_highpass.
function ptspa_highpass_Callback(hObject, eventdata, handles)
% hObject    handle to ptspa_highpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ptspa_highpass
if get(hObject, 'Value')==1
    set(handles.ptspa_co_frequency, 'Visible', 'On');
    set(handles.text133, 'Visible', 'On');
else
    set(handles.ptspa_co_frequency, 'Visible', 'Off');
    set(handles.text133, 'Visible', 'Off');
end

end


function ptspa_co_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to ptspa_co_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptspa_co_frequency as text
%        str2double(get(hObject,'String')) returns contents of ptspa_co_frequency as a double
end

% --- Executes during object creation, after setting all properties.
function ptspa_co_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptspa_co_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in ptstr_highpass.
function ptstr_highpass_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_highpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ptstr_highpass
if get(hObject, 'Value')==1
    set(handles.ptstr_co_frequency, 'Visible', 'On');
    set(handles.text134, 'Visible', 'On');
else
    set(handles.ptstr_co_frequency, 'Visible', 'Off');
    set(handles.text134, 'Visible', 'Off');
end

end


function ptstr_co_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_co_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptstr_co_frequency as text
%        str2double(get(hObject,'String')) returns contents of ptstr_co_frequency as a double
end

% --- Executes during object creation, after setting all properties.
function ptstr_co_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptstr_co_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function mph_windowsize_Callback(hObject, eventdata, handles)
% hObject    handle to mph_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mph_windowsize as text
%        str2double(get(hObject,'String')) returns contents of mph_windowsize as a double
end

% --- Executes during object creation, after setting all properties.
function mph_windowsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mph_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function mvh_windowsize_Callback(hObject, eventdata, handles)
% hObject    handle to mvh_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mvh_windowsize as text
%        str2double(get(hObject,'String')) returns contents of mvh_windowsize as a double
end

% --- Executes during object creation, after setting all properties.
function mvh_windowsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mvh_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptspa_windowsize_Callback(hObject, eventdata, handles)
% hObject    handle to ptspa_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptspa_windowsize as text
%        str2double(get(hObject,'String')) returns contents of ptspa_windowsize as a double
end

% --- Executes during object creation, after setting all properties.
function ptspa_windowsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptspa_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ptstr_windowsize_Callback(hObject, eventdata, handles)
% hObject    handle to ptstr_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptstr_windowsize as text
%        str2double(get(hObject,'String')) returns contents of ptstr_windowsize as a double
end

% --- Executes during object creation, after setting all properties.
function ptstr_windowsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptstr_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in mph_red.
function mph_red_Callback(hObject, eventdata, handles)
% hObject    handle to mph_red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.mph_red, 'Value')==1
    set(handles.mph_red_opt, 'Visible', 'On');
else
    set(handles.mph_red_opt, 'Visible', 'Off');
end

% Hint: get(hObject,'Value') returns toggle state of mph_red
end


function mph_red_opt_Callback(hObject, eventdata, handles)
% hObject    handle to mph_red_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mph_red_opt as text
%        str2double(get(hObject,'String')) returns contents of mph_red_opt as a double
end

% --- Executes during object creation, after setting all properties.
function mph_red_opt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mph_red_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in mvh_red.
function mvh_red_Callback(hObject, eventdata, handles)
% hObject    handle to mvh_red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mvh_red

if get(handles.mvh_red, 'Value')==1
    set(handles.mvh_red_opt, 'Visible', 'On');
else
    set(handles.mvh_red_opt, 'Visible', 'Off');
end

end

function mvh_red_opt_Callback(hObject, eventdata, handles)
% hObject    handle to mvh_red_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mvh_red_opt as text
%        str2double(get(hObject,'String')) returns contents of mvh_red_opt as a double
end

% --- Executes during object creation, after setting all properties.
function mvh_red_opt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mvh_red_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in lst_analyze_flag.
function lst_analyze_flag_Callback(hObject, eventdata, handles)
global strct_config

if get(handles.lst_analyze_flag, 'Value')==1
    set(handles.alpha_start, 'Visible', 'On');
    set(handles.alpha_end, 'Visible', 'On');
    set(handles.alpha_steps, 'Visible', 'On');
    set(handles.text161, 'Visible', 'On');
    set(handles.text162, 'Visible', 'On');
    set(handles.text163, 'Visible', 'On');
else
    set(handles.alpha_start, 'Visible', 'Off');
    set(handles.alpha_end, 'Visible', 'Off');
    set(handles.alpha_steps, 'Visible', 'Off');
    set(handles.text161, 'Visible', 'Off');
    set(handles.text162, 'Visible', 'Off');
    set(handles.text163, 'Visible', 'Off');
end

strct_config.alpha_start=str2double(get(handles.alpha_start, 'String'));
strct_config.alpha_end=str2double(get(handles.alpha_end, 'String'));
strct_config.alpha_steps=str2double(get(handles.alpha_steps, 'String'));

end

% --- Executes on button press in pta_export_profiles.
function pta_export_profiles_Callback(hObject, eventdata, handles)
   
end

function yplus_ptstr_Callback(hObject, eventdata, handles)
% hObject    handle to yplus_ptstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yplus_ptstr as text
%        str2double(get(hObject,'String')) returns contents of yplus_ptstr as a double
end

% --- Executes during object creation, after setting all properties.
function yplus_ptstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yplus_ptstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in side_selection.
function side_selection_Callback(hObject, eventdata, handles)
% hObject    handle to side_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns side_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from side_selection
end

% --- Executes during object creation, after setting all properties.
function side_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to side_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in L_ref_flag.
function L_ref_flag_Callback(hObject, eventdata, handles)
global strct_config

if get(handles.L_ref_flag, 'Value')==1
    set(handles.L_ref_text, 'Visible', 'On');
    set(handles.text160, 'Visible', 'On');
else
    set(handles.L_ref_text, 'Visible', 'Off');
    set(handles.text160, 'Visible', 'Off');
end

strct_config.global_L_ref = str2double(get(handles.L_ref_text, 'String'));
strct_config.L_ref_flag = (get(handles.L_ref_flag, 'Value'));

end


function L_ref_text_Callback(hObject, eventdata, handles)

global strct_config

strct_config.global_L_ref=str2double(get(handles.L_ref_text, 'String'));

end

% --- Executes during object creation, after setting all properties.
function L_ref_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_ref_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function alpha_end_Callback(hObject, eventdata, handles)
global strct_config

strct_config.alpha_end=str2double(get(handles.alpha_end, 'String'));
end

% --- Executes during object creation, after setting all properties.
function alpha_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function alpha_steps_Callback(hObject, eventdata, handles)
global strct_config

strct_config.alpha_steps=str2double(get(handles.alpha_steps, 'String'));
end

% --- Executes during object creation, after setting all properties.
function alpha_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function alpha_start_Callback(hObject, eventdata, handles)

global strct_config

strct_config.alpha_start=str2double(get(handles.alpha_start, 'String'));
end

% --- Executes during object creation, after setting all properties.
function alpha_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in mph_samples.
function mph_samples_Callback(hObject, eventdata, handles)
% hObject    handle to mph_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mph_samples

if get(handles.mph_samples, 'Value')==1
    set(handles.mph_samples_opt, 'Visible', 'On');
else
    set(handles.mph_samples_opt, 'Visible', 'Off');
end

end


function mph_samples_opt_Callback(hObject, eventdata, handles)
% hObject    handle to mph_samples_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mph_samples_opt as text
%        str2double(get(hObject,'String')) returns contents of mph_samples_opt as a double
end

% --- Executes during object creation, after setting all properties.
function mph_samples_opt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mph_samples_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in mvh_samples.
function mvh_samples_Callback(hObject, eventdata, handles)
% hObject    handle to mvh_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mvh_samples
if get(handles.mvh_samples, 'Value')==1
    set(handles.mvh_samples_opt, 'Visible', 'On');
else
    set(handles.mvh_samples_opt, 'Visible', 'Off');
end
end


function mvh_samples_opt_Callback(hObject, eventdata, handles)
% hObject    handle to mvh_samples_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mvh_samples_opt as text
%        str2double(get(hObject,'String')) returns contents of mvh_samples_opt as a double
end

% --- Executes during object creation, after setting all properties.
function mvh_samples_opt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mvh_samples_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
