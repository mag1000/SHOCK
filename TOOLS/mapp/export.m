function export(handles,varnames,varnames_rohdaten,parameter,daten,rohdaten,exportType)
    if isempty(handles)==1
        fprintf('Beginne Export in Datei: %s...\n',parameter.path);
    else
        logbar(handles.status_bar,'Beginne Export...');
    end
    
    %Export Daten
    if exportType==0
    
        file1=fopen(parameter.path,'w');
        fprintf(file1,'TITLE = "%s"\n',parameter.title);
        fprintf(file1,'VARIABLES =');
        fprintf(file1,' "%s"',varnames{1});
        for probe=1:1:length(daten)
            fprintf(file1,' "%s_probe%d"',varnames{2},probe);
        end
        fprintf(file1,'\n');

        fprintf(file1,'ZONE T=\"%s\", F=BLOCK, I=%d, DT=(SINGLE) \n',parameter.methodename,numel(daten{1}(1,:)));
        fprintf(file1,'%e\n',daten{1}(1,:));
        for probe=1:1:length(daten)
            fprintf(file1,'%e\n',daten{probe}(2,:));
        end   
        fclose(file1);
    
    elseif exportType==1 %ExportTyp z.B. für STD
        
        file1=fopen(parameter.path,'w');
        fprintf(file1,'TITLE = "%s"\n',parameter.title);
        fprintf(file1,'VARIABLES =');
        for var=1:1:length(varnames)
            fprintf(file1,' "%s"',varnames{var});
        end
        fprintf(file1,'\n');

        fprintf(file1,'ZONE T=\"%s\", F=BLOCK, I=%d, DT=(SINGLE) \n',parameter.methodename,numel(daten));
        for probe=1:1:length(daten)
            fprintf(file1,'%e\n',daten{probe}(1,:));
        end   
        for probe=1:1:length(daten)
            fprintf(file1,'%e\n',daten{probe}(2,:));
        end   
        fclose(file1);
        
    end
    
    %Export Rohdaten
    file2=fopen(parameter.path_rohdaten,'w');
    fprintf(file2,'TITLE = "%s"\n',parameter.title_rohdaten);
    fprintf(file2,'VARIABLES =');
    fprintf(file2,' "%s"',varnames_rohdaten{1});
    for probe=1:1:length(daten)
        fprintf(file2,' "%s_probe%d"',varnames_rohdaten{2},probe);
    end
    fprintf(file2,'\n');    
    fprintf(file2,'ZONE T=\"%s_Rohdaten\", F=BLOCK, I=%d, DT=(SINGLE) \n',parameter.methodename,numel(rohdaten{1}(1,:)));
    fprintf(file2,'%e\n',rohdaten{1}(1,:));
    for probe=1:1:length(rohdaten)
        fprintf(file2,'%e\n',rohdaten{probe}(2,:));
    end       
    fclose(file2);	
        
        
    if isempty(handles)==1
        fprintf('fertig\n');
    else
        logbar(handles.status_bar,'fertig');
    end
end

