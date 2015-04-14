function varargout = frequenzanalyse(varargin)
% FREQUENZANALYSE MATLAB code for frequenzanalyse.fig
%      FREQUENZANALYSE, by itself, creates a new FREQUENZANALYSE or raises the existing
%      singleton*.
%
%      H = FREQUENZANALYSE returns the handle to a new FREQUENZANALYSE or the handle to
%      the existing singleton*.
%
%      FREQUENZANALYSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FREQUENZANALYSE.M with the given input arguments.
%
%      FREQUENZANALYSE('Property','Value',...) creates a new FREQUENZANALYSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before frequenzanalyse_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to frequenzanalyse_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help frequenzanalyse

% Last Modified by GUIDE v2.5 10-Jul-2014 11:24:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @frequenzanalyse_OpeningFcn, ...
                   'gui_OutputFcn',  @frequenzanalyse_OutputFcn, ...
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

% --- Executes just before frequenzanalyse is made visible.
function frequenzanalyse_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to frequenzanalyse (see VARARGIN)

% Choose default command line output for frequenzanalyse
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes frequenzanalyse wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global strct_config 

if isempty(strct_config)    
    strct_config.unit_y = 1;
    strct_config.unit_x = 1;
    strct_config.FlagMappCallFA=0;
else
    strct_config.FlagMappCallFA=1;
    set(handles.add_file, 'Visible', 'Off');
    set(handles.select_unit_x, 'Visible', 'Off');
    set(handles.text_path, 'Visible', 'Off');
    set(handles.select_unit_y, 'Visible', 'Off');
    set(handles.text16, 'Visible', 'Off');
    set(handles.text15, 'Visible', 'Off');
    set(handles.start_import, 'Visible', 'Off');
    set(handles.import_information, 'Visible', 'On');
    start_import_Callback(handles.start_import, eventdata, handles);    
    return;

end

end

% --- Outputs from this function are returned to the command line.
function varargout = frequenzanalyse_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%clear 
end

% --- Executes on button press in add_file.
function add_file_Callback(hObject, eventdata, handles)
global strct_config
global text_path_handle


[strct_config.RawDataFileName, strct_config.folder] = uigetfile(sprintf('%s/*.dat', 'home'), 'Select Data');

strct_config.RawDataPath = [strct_config.folder, strct_config.RawDataFileName];
text=['Datei ausgewählt: ', strct_config.RawDataPath ] ;
logbar(handles.status_bar,text);

set(text_path_handle, 'string', strct_config.RawDataPath)

end

function text_path_Callback(hObject, eventdata, handles)
% hObject    handle to text_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_path as text
%        str2double(get(hObject,'String')) returns contents of text_path as a double


% --- Executes during object creation, after setting all properties.

path = get(hObject, 'string') ;
text=['Pfad geändert: ', path ] ;
logbar(handles.status_bar,text);

end

function text_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global index_data
index_data=1;

global text_path_handle ;
text_path_handle= hObject;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function select_unit_x_Callback(hObject, eventdata, handles)

global strct_config

contents = cellstr(get(hObject,'String')) ;
choice = contents{get(hObject,'Value')} ;

switch choice
    case 'Zeit [s]'
        strct_config.unit_x = 1 ;
    case 'Position [x]'
        strct_config.unit_x = 2 ;
        set(handles.select_unit_y, 'Value', 2)
    case 'Position [z]'
        strct_config.unit_x = 3 ;
        set(handles.select_unit_y, 'Value', 2)
end

logbar(handles.status_bar,['Einheit X-Achse: ', choice]);

end

function select_unit_y_Callback(hObject, eventdata, handles)

global strct_config

contents = cellstr(get(hObject,'String')) ;
choice = contents{get(hObject,'Value')} ;

switch choice
    case 'p'
        strct_config.unit_y = 1 ;
    case 'TKE'
        strct_config.unit_y = 2 ;
end

logbar(handles.status_bar,['Einheit Y-Achse: ', choice]);

end

function select_unit_y_CreateFcn(hObject, eventdata, handles)


end

function select_unit_x_CreateFcn(hObject, eventdata, handles)

end


% --- Executes on button press in start_import.
function start_import_Callback(hObject, eventdata, handles)

logbar(handles.status_bar,'Starte Import...');
global strct_probes
global strct_config

importData(handles);
strct_config.updated=1;

inf=sprintf('Importierte Datei: %s \nEinheit X: %d - Einheit Y: %d - Anzahl probes: %d - Anzahl samples: %d',...
    strct_config.RawDataPath, strct_config.unit_x, strct_config.unit_y,strct_config.probes,strct_probes(1).samples);
set(handles.import_information, 'String', inf);

%%Probenfenster aktualisieren
for probe=1:strct_config.probes
    contents{probe}=strct_probes(probe).name;  
end
set(handles.probes_menu, 'String', contents) ;  

%Paramtermatrix initialisieren
logbar(handles.status_bar,'Parameter konfigurieren...');

%%% Werte setzten, falls nicht durch Mapp gestartet
strct_config.lastSample= strct_probes(1).samples; %letztes Sample

if strct_config.FlagMappCallFA==0
    strct_config.yPlusSoll=sprintf('1 5 10 25 35 50 70 100 120 150 170 200');
    strct_config.hanning = 1 ; %hanning window flag
    strct_config.avrg = 1 ; %mittelwert flag
    strct_config.firstSample= 1; %erstes Sample
    strct_config.ar= 30; %AR-Ordnung f�r pBurg 
else % Auswahl der probes setzung und Exportieren
    if strct_config.unit_x==2 || strct_config.unit_x==3
        select_yplus_Callback(0, 0, handles);
    elseif strct_config.unit_x==1 || strct_config.unit_x==4
        for i=1:strct_config.probes
            markiert(i)=i;
        end
        strct_config.markiert=markiert;
        set(handles.probes_menu, 'Value', markiert);
    end
    
    export_tecplot_Callback(0,0,handles);
    %clear all;
    %delete(handles.figure1);   
    return ;
end

 
%Slider aktualisieren
set(handles.sliderFirstSample, 'Min', 1, 'Max', strct_config.lastSample);
set(handles.sliderLastSample, 'Min', 1, 'Max', strct_config.lastSample);

%Markierung aktualisieren
set(handles.probes_menu, 'Value', 1) ;
strct_config.markiert = get(handles.probes_menu, 'Value');

%Methode aktualisieren
set(handles.methode, 'Value', 1) ;
strct_config.methode=get(handles.methode, 'Value') ;

%Parameter setzen
update_parameter(handles);

%Plot aktualisieren
update_plot(handles) ;

%yPLus Eingabefeld einblenden
if strct_config.unit_x==2 || strct_config.unit_x==3
    set(handles.text11, 'Visible', 'On')
    set(handles.yPlus, 'Visible', 'On')
    set(handles.select_yplus, 'Visible', 'On');
end

end

% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)

global index_data imported_data parameter

contents = cellstr(get(handles.dataset_menu,'String')) ;
selected = get(handles.dataset_menu,'Value') ;

if selected < index_data-1
    for i=selected:index_data-2
        imported_data{i}=imported_data{i+1} ;
        parameter(:,:,1,i)=parameter(:,:,1,i+1);
    end
end    

contents{index_data-1} = '';
if index_data==2;
   inhalt=[];
else
for i=1:index_data-2
   inhalt{i}=contents{i};
end  
end

imported_data{index_data-1}='';
%parameter(:,:,1, index_data-1)="";
index_data=index_data-1;

%%Datensätze aktualisieren
set(handles.dataset_menu, 'String', inhalt) ;
set(handles.dataset_menu, 'Value', 1);

%%Proben Menu aktualisieren
selected_dataset = get(handles.dataset_menu, 'Value') ;
if selected_dataset>1
    for i=1:parameter(1,1,1,selected_dataset)
        text=['Probe ',num2str(i)]; 
        contents{i}=text;  
    end
%Markierung aktualisieren
    l=1; 
    for i=1:parameter(1,1,1,selected_dataset)
        if parameter(i+1,1,1,selected_dataset)==1
          markiert(l)=i;
          l=l+1;
        end
    end
    else
    markiert=[];
    contents=[];
    
end
set(handles.probes_menu, 'String', contents) ; 
set(handles.probes_menu, 'Value', markiert) ;
if index_data > 1
    update_parameter(handles) ;
end
update_plot(handles) ;

end

% --- Executes on selection change in probes_menu.
function probes_menu_Callback(hObject, eventdata, handles)

global strct_config 

strct_config.markiert = get(handles.probes_menu, 'Value');

update_plot(handles) ;
update_parameter(handles) ;

end

% --- Executes during object creation, after setting all properties.
function probes_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to probes_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function sample_red_Callback(hObject, ~ , handles)
global parameter imported_data

selected_dataset = get(handles.dataset_menu, 'Value');
markiert = get(handles.probes_menu, 'Value');

for i=1:parameter(1,1,1,selected_dataset)
    parameter(i+1,6,strct_config.methode+1,selected_dataset)=1;
end

parameter(1,5,strct_config.methode+1,selected_dataset)=str2double(get(hObject, 'String'));

set(handles.sliderFirstSample, 'Value', parameter(markiert(1)+1,2,strct_config.methode+1,selected_dataset));
set(handles.sliderLastSample, 'Value', parameter(markiert(1)+1,3,strct_config.methode+1,selected_dataset));

update_plot(handles) ;
parameter(markiert(1)+1,2,strct_config.methode+1,selected_dataset)=1;
parameter(markiert(1)+1,3,strct_config.methode+1,selected_dataset)=parameter(1,2,1,selected_dataset);
update_parameter(handles);
end

% --- Executes during object creation, after setting all properties.
function sample_red_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in select_all.
function select_all_Callback(hObject, eventdata, handles)
global strct_config

for probe=1:strct_config.probes
    strct_config.markiert(probe)=probe;
end
set(handles.probes_menu, 'Value', strct_config.markiert) ;

update_plot(handles) ;
update_parameter(handles) ;
end


% Slider 
function sliderFirstSample_CreateFcn(hObject, eventdata, handles)
end

function sliderLastSample_CreateFcn(hObject, eventdata, handles)
end


% --- Executes on slider movement.
function sliderFirstSample_Callback(hObject, eventdata, handles)
global strct_config
position=get(hObject, 'Value');
position=round(position);
strct_config.firstSample=position;
update_plot(handles) ;
update_parameter(handles);
end


% --- Executes on slider movement.
function sliderLastSample_Callback(hObject, eventdata, handles)
global strct_config
position=get(hObject, 'Value');
position=round(position);
strct_config.lastSample=position;
update_parameter(handles);
update_plot(handles) ;
end

% --- Executes during object creation, after setting all properties.
function sliderFirstSample_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderFirstSample_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end
% --- Executes during object creation, after setting all properties.
function sliderLastSample_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderLastSample_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end


function ar_order_Callback(hObject, eventdata, handles)
global strct_config
value=str2double(get(hObject, 'String'));
strct_config.ar=value;
update_plot(handles) ;

end


% --- Executes during object creation, after setting all properties.
function ar_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ar_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in calc_ar_order.
function calc_ar_order_Callback(hObject, eventdata, handles)

ordnung=str2double(get(handles.ar_order, 'String'));

update_plot(handles);
end


% --- Executes on button press in hanning_flag.
function hanning_flag_Callback(hObject, eventdata, handles)

global strct_config
value=get(hObject, 'Value');

strct_config.hanning=value;

update_plot(handles) ;
end

% --- Executes on button press in average_flag.
function average_flag_Callback(hObject, eventdata, handles)
global strct_config
value=get(hObject, 'Value');

strct_config.avrg=value;

update_plot(handles) ;

end


function status_bar_Callback(hObject, eventdata, handles)
% hObject    handle to status_bar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of status_bar as text
%        str2double(get(hObject,'String')) returns contents of status_bar as a double
end


% --- Executes during object creation, after setting all properties.
function status_bar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status_bar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global text1 ;
text1 = '';
global status_text ;
status_text = hObject;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function number_s_Callback(hObject, eventdata, handles)
% hObject    handle to number_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_s as text
%        str2double(get(hObject,'String')) returns contents of number_s as a double
end

% --- Executes during object creation, after setting all properties.
function number_s_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function number_w_Callback(hObject, eventdata, handles)
% hObject    handle to number_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_w as text
%        str2double(get(hObject,'String')) returns contents of number_w as a double
end

% --- Executes during object creation, after setting all properties.
function number_w_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function redelta_start_Callback(hObject, eventdata, handles)
% hObject    handle to redelta_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of redelta_start as text
%        str2double(get(hObject,'String')) returns contents of redelta_start as a double
end

% --- Executes during object creation, after setting all properties.
function redelta_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to redelta_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function redelta_ende_Callback(hObject, eventdata, handles)
% hObject    handle to redelta_ende (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of redelta_ende as text
%        str2double(get(hObject,'String')) returns contents of redelta_ende as a double
end

% --- Executes during object creation, after setting all properties.
function redelta_ende_CreateFcn(hObject, eventdata, handles)
% hObject    handle to redelta_ende (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in methode.
function methode_Callback(hObject, eventdata, handles)
global strct_config

strct_config.methode=get(handles.methode, 'Value') ;
update_plot(handles);

% hObject    handle to methode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methode
end

% --- Executes during object creation, after setting all properties.
function methode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in export_tecplot.
function export_tecplot_Callback(hObject, eventdata, handles)
global strct_config
global strct_probes

exportType=0; %StandardExport

if strct_config.methode==1
    methodename='FFT';
    Yname='Amplitude';
end
if strct_config.methode==2
    methodename='Spektogram';
    Yname='Amplitude';
end
if strct_config.methode==3
    methodename='PWelch';
    Yname='Amplitude';
end
if strct_config.methode==4
    methodename='1_3Oktave';
    Yname='Sound Pressure Level';
end
if strct_config.methode==5
    methodename='PBurg';
    Yname='Power Spectral Density';
end
if strct_config.methode==6
    methodename='STD';
    Yname='Standard deviation';
    exportType=1;
end
if strct_config.methode==7
    methodename='LST';
    Yname='alpha_i [1/m]';
end


    parameter.path=sprintf('%s/%s_%s.dat',strct_config.folder,strrep(strct_config.RawDataFileName,'.dat',''),methodename);
    parameter.path_rohdaten=sprintf('%s/%s_%s_rohdaten.dat',strct_config.folder,strrep(strct_config.RawDataFileName,'.dat',''),methodename);
    parameter.title=sprintf('%s_%s.dat',strrep(strct_config.RawDataFileName,'.dat',''),methodename);
    parameter.title_rohdaten=sprintf('%s_%s_rohdaten.dat',strrep(strct_config.RawDataFileName,'.dat',''),methodename);
    parameter.methodename=methodename;
    varnames{2}=Yname;    
    
    if strct_config.unit_x==1 %Temporal
        parameter.all=0; %export one zone containing all probes
        varnames{1}='Frequenz [Hz]';
        varnames_rohdaten{1}='Zeit [s]';
        if strct_config.unit_y==1
            varnames_rohdaten{2}='p/p<sub>0</sub> [-]';
        elseif strct_config.unit_y==2
            varnames_rohdaten{2}='TKE [-]';
        end
    elseif strct_config.unit_x==2 %Streamwise
        parameter.all=0;
        varnames{1}='Wellenzahl k<sub>x</sub>';    
        varnames_rohdaten{1}='x/c [-]';
        varnames_rohdaten{2}='TKE [-]';
    elseif strct_config.unit_x==3 %Spanwise
        parameter.all=0;
        varnames{1}='Wellenzahl k<sub>z</sub>';
        varnames_rohdaten{1}='x/c [-]';
        varnames_rohdaten{2}='TKE [-]';
    end
    
    if strct_config.methode==6 %Standard deviation
        parameter.all=1;
        varnames{1}='x/c [-]';
    end
    if strct_config.methode==7 %LST
        parameter.all=0;
        varnames{1}='Frequenz [Hz]';
    end
    
    %Daten schreiben; FA neu berechnen
    for probe=1:length(strct_config.markiert)
        parameter.zonename{probe}=strct_probes(strct_config.markiert(probe)).name;
        if strct_config.updated == 1
            process_signal(strct_config.markiert(probe));
            statistics(handles,strct_config.markiert(probe));
        end
        
        daten{probe}(1,:)=strct_probes(strct_config.markiert(probe)).psdX;
        daten{probe}(2,:)=strct_probes(strct_config.markiert(probe)).psdY;
        
        rohdaten{probe}(1,:)=strct_probes(strct_config.markiert(probe)).processedDataX;
        rohdaten{probe}(2,:)=strct_probes(strct_config.markiert(probe)).processedDataY;
        
    end
    
    export(handles,varnames,varnames_rohdaten,parameter,daten,rohdaten,exportType);
    fprintf('Konvertiere dat->plt und Lösche im Anschluss dat-Datei\n');
    command=sprintf('cd %s && /usr/local/Tecplot2015R1/bin/preplot %s && rm -f %s',strct_config.folder,parameter.path,parameter.path);
    [status, cmdout] = system(command);
    
    command=sprintf('cd %s && /usr/local/Tecplot2015R1/bin/preplot %s && rm -f %s',strct_config.folder,parameter.path_rohdaten,parameter.path_rohdaten);
    [status, cmdout] = system(command);    

    strct_config.updated=0;
end

% --- Executes on button press in export_plot.
function export_plot_Callback(hObject, eventdata, handles)
% hObject    handle to export_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


function yPlus_Callback(hObject, ~, handles)

end

% --- Executes during object creation, after setting all properties.
function yPlus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yPlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function ma_Callback(hObject, eventdata, handles)
% hObject    handle to ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ma as text
%        str2double(get(hObject,'String')) returns contents of ma as a double
end

% --- Executes during object creation, after setting all properties.
function ma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function re_Callback(hObject, eventdata, handles)
% hObject    handle to re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of re as text
%        str2double(get(hObject,'String')) daten{i}(1,:)=results{selected_dataset,1,1};returns contents of re as a double

end

% --- Executes during object creation, after setting all properties.
function re_CreateFcn(hObject, eventdata, handles)
% hObject    handle to re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function update_plot(handles)
%aktualisiert beide Plots passend zu den markierten Proben; 
%falls sich seit dem letzten Update Parameter geändert haben, wird vorher die Frequenzanalyse aktualisiert;
%sind mehrere Proben ausgewählt, werden die Plots Übereinandergelegt
%%%
%%% !!! results{datensatz,probe+1,methode}; probe=1: frequenz
%%% !!! signal{datensatz,pobe,methode,Achse}; Achse=1:X, Achse=2:Y

global strct_probes
global strct_config

logbar(handles.status_bar,'Plot aktualisieren...');

cla(handles.axes1,'reset')
cla(handles.axes3,'reset')

for i=1:length(strct_config.markiert)
    process_signal(strct_config.markiert(i));
    statistics(handles,strct_config.markiert(i));

    hold on;
    axes(handles.axes1);
    plot(strct_probes(strct_config.markiert(i)).processedDataX,strct_probes(strct_config.markiert(i)).processedDataY);
    axis 'auto xy';
    hold on;
    axes(handles.axes3);
    plot(strct_probes(strct_config.markiert(i)).psdX,strct_probes(strct_config.markiert(i)).psdY);
    axis([0 100 0 1]);
    axis 'auto y';
end


end

function update_parameter(handles)
% setzt die angezeigten Einstellungen in der GUI passend zur markierten Probe


global strct_config
global strct_probes

%Einstellung unter Signalanalyse
set(handles.hanning_flag, 'Value', strct_config.hanning);
set(handles.average_flag, 'Value', strct_config.avrg);
set(handles.yPlus, 'String', strct_config.yPlusSoll);

%Slider aktualisieren
set(handles.sliderFirstSample, 'Value', round(strct_config.firstSample));
set(handles.sliderLastSample, 'Value', round(strct_config.lastSample));
set(handles.sliderFirstSample_text, 'String', round(strct_config.firstSample));
set(handles.sliderLastSample_text, 'String', round(strct_config.lastSample));
set(handles.ar_order, 'String', strct_config.ar);

if strct_config.unit_x == 2 || strct_config.unit_x == 3 %prüfe Einheit X
    set(handles.text11, 'Visible', 'On');
    set(handles.yPlus, 'Visible', 'On');
else 
    set(handles.text11, 'Visible', 'Off');
    set(handles.yPlus, 'Visible', 'Off');
end


%Koordinaten und Name setzen
set(handles.coordinateX_text, 'String', strct_probes(strct_config.markiert(1)).CoordinateX);
set(handles.coordinateY_text, 'String', strct_probes(strct_config.markiert(1)).CoordinateY);
set(handles.names_text, 'String', strct_probes(strct_config.markiert(1)).name);
end



function coordinateX_text_Callback(hObject, eventdata, handles)
% hObject    handle to coordinateX_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%selected_dataset = get(handles.dataset_menu, 'Value'); Hints: get(hObject,'String') returns contents of coordinateX_text as text
%        str2double(get(hObject,'String')) returns contents of coordinateX_text as a double
end



function names_text_Callback(hObject, eventdata, handles)
% hObject    handle to names_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of names_text as text
%        str2double(get(hObject,'String')) returns contents of names_text as a double
end



function coordinateY_text_Callback(hObject, eventdata, handles)
% hObject    handle to coordinateY_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coordinateY_text as text
%        str2double(get(hObject,'String')) returns contents of
%        coordinateY_text as a double
end


% --- Executes during object creation, after setting all properties.
function coordinateY_text_CreateFcn(hObject, eventdata, handles)


% hObject    handle to coordinateY_text (see GCBO)
%selected_dataset = get(handles.dataset_menu, 'Value'); eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes during object creation, after setting all properties.
function coordinateX_text_CreateFcn(hObject, eventdata, handles)

% hObject    handle to coordinateX_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes during object creation, after setting all properties.
function names_text_CreateFcn(hObject, eventdata, handles)

% hObject    handle to names_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in select_c_tools.
function select_c_tools_Callback(hObject, eventdata, handles)



end

% --- Executes on button press in select_yplus.
function select_yplus_Callback(hObject, eventdata, handles)
global strct_config

if isempty(strct_config.FlagMappCallFA)
defAns{1}=strct_config.yPlusSoll;
x = inputdlg('Y+ Werte setzen:','Y+', [1 50],defAns);

strct_config.yPlusSoll = x{:} ;
end

wert=str2num(strct_config.yPlusSoll);
strct_config.markiert=zeros(length(wert),1);
for i=1:length(wert)
    strct_config.markiert(i)=find_probe_yplus(handles,wert(i)) ;
end

%Markierung aktualisieren
if isempty(handles) == 1
else
set(handles.probes_menu, 'Value', strct_config.markiert) ;

update_plot(handles);
update_parameter(handles);
end

end
