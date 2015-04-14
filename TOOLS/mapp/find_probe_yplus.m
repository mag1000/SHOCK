function yPlus_probe=find_probe_yplus(handles,yPlusSoll)
global strct_config
global strct_probes

if isempty(handles) == 1
    fprintf('Suche probe für yPlus=%d...\n',yPlusSoll);
else
    logbar(handles.status_bar,sprintf('Suche probe für yPlus=%d... ',yPlusSoll));
end

distance=999999999.;
for probe=1:strct_config.probes
    distance_tmp=abs(strct_probes(probe).yPlus-yPlusSoll);
    if distance_tmp<distance
        distance=distance_tmp;
        yPlus_probe=probe;
        yPlus_value=strct_probes(probe).yPlus;
    end
end
if isempty(handles) == 1
    fprintf(' ...gefundene Übereinstimmung für probe %d bei yPlus:%d\n',yPlus_probe,yPlus_value);
else
    logbar(handles.status_bar,sprintf(' ...gefundene Übereinstimmung für probe %d bei yPlus:%d',yPlus_probe,yPlus_value));
end

end