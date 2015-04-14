function logbar(bar,add_text)
%     alter_text=cellstr(get(bar,'String'));
%     set(bar,'String',[alter_text;cellstr(add_text)] );
%     drawnow();    

fprintf('%s\n',add_text);

%     jhEdit = findjobj(bar);    
%     jVScroll = jhEdit.getVerticalScrollBar;
%     jVScroll.setValue(jVScroll.getMaximum);
%     jhEdit.repaint;
end