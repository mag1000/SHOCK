%Importiert den Datensatz und speichert alle relevaten Informationen in den
%entsprechenden structuren

function importData(handles)
    global strct_probes
    global strct_config
    ProbesCountFlag=0;
    switch strct_config.unit_x
        case 1
            if strct_config.unit_y == 2
            specificCharacter='TKE';
            strct_config.probes = strct_config.TKEprobes;
            else 
            specificCharacter='PH';
            strct_config.probes = strct_config.PHprobes;
            end
            flag_yPlus=0;
        case 2
            specificCharacter='TKE';
            ProbesCountFlag=1;
            flag_yPlus=1;
        case 3
            specificCharacter='TKE';
            flag_yPlus=1;
            ProbesCountFlag=1;
%         case 4
%             specificCharacter='PH'; 
%             flag_yPlus=0;
    end
    fprintf('Lade Datei: %s\n',char(strct_config.RawDataPath));
    file=fopen(strct_config.RawDataPath,'r');
    header = fgetl(file);
    header = fgetl(file);
    if ProbesCountFlag == 1
        strct_config.probes=length(findstr(header,'TKE'));
    end
    if isempty (handles) == 1
        fprintf('Lese %d probes ein...\n',strct_config.probes);
    else
        logbar(handles.status_bar,sprintf('Lese %d probes ein...',strct_config.probes));
    end
    
    format_start='%*s %*s %*s %*s';

    strct_probes(strct_config.probes+1).name=cellstr('averaged_sum');
    yPlus=zeros(strct_config.probes,1);
    for index=1:strct_config.probes
        format=sprintf('%s %%s,1',format_start);
        strct_probes(index).name=strrep(sscanf(header, format),'"','');
        format_start=sprintf('%s %%*s',format_start);
        string1=char(strct_probes(index).name);
        if flag_yPlus==1
            string2=sprintf('%s_%%*d_x%%e_y%%*e_yPlus%%*e',specificCharacter);
            strct_probes(index).CoordinateX=sscanf(string1,string2);
            string2=sprintf('%s_%%*d_x%%*e_y%%e_yPlus%%*e',specificCharacter);
            strct_probes(index).CoordinateY=sscanf(string1,string2);
            string2=sprintf('%s_%%*d_x%%*e_y%%*e_yPlus%%e',specificCharacter);
            strct_probes(index).yPlus=round(sscanf(string1,string2));
            strct_probes(index).name=sprintf('%s_%d_x%.2f_y%.4f_yPlus%d',specificCharacter,index,strct_probes(index).CoordinateX,strct_probes(index).CoordinateY,strct_probes(index).yPlus);
        else
            string2=sprintf('%s_%%*d_x%%e_y%%*e',specificCharacter);
            strct_probes(index).CoordinateX=sscanf(string1,string2);
            string2=sprintf('%s_%%*d_x%%*e_y%%e',specificCharacter);
            strct_probes(index).CoordinateY=sscanf(string1,string2);
        end
    end
    
    header = fgetl(file);
    samples=sscanf(header, '%*s %*s %*s I=%d,');
    
    %IMPORT-PROZESS
    tmp=fscanf(file,'%le\n',(strct_config.probes+1)*samples);

    for index=1:strct_config.probes
        strct_probes(index).rawDataX=tmp(1:samples);
        strct_probes(index).rawDataY=tmp((index)*samples+1:(index+1)*samples);

        strct_probes(index).samples=samples;
        strct_probes(index).dx=(strct_probes(index).rawDataX(samples)-strct_probes(index).rawDataX(1))/samples;
    end
    strct_probes(index).rawDataY(:)=strct_probes(index).rawDataY(:)-strct_probes(index).rawDataY(1);
end