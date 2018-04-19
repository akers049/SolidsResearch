function slider2_Callback(hObject, eventdata)

step_length = 30;

    if(isempty(hObject.UserData))
        index = ceil(get(hObject,'Value'));
    elseif (hObject.UserData < get(hObject, 'Value'))
        index = ceil(get(hObject,'Value'));
    else
        index = floor(get(hObject, 'Value'));
    end
    set(hObject, 'UserData', index);
    
    set(hObject, 'Value', index);
    f = figure(1);
    bgcolor = f.Color;
    uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String',string(index),'BackgroundColor',bgcolor);
    
    color_black = 'k.';
    color_red = 'r.';
            
    load('eigenvalues_data.mat');
    

    wavelength_black = [];
    eigs_red = [];
    eigs_black = [];
    wavelength_red = [];
    for i = 1:length(data_cell{index})
        for k = 2:length(data_cell{index}{i}(2:end))
            if (data_cell{index}{i}(k) < -5e-6)
                wavelength_red = [wavelength_red; data_cell{index}{i}(1)]; 
                eigs_red = [eigs_red;data_cell{index}{i}(k)];
                
                % plot(data_cell{index}{i}(1),data_cell{index}{i}(k), color_red, 'MarkerSize', 7);
            else
                wavelength_black = [wavelength_black; data_cell{index}{i}(1)];
                eigs_black = [eigs_black;data_cell{index}{i}(k)];
               % plot(data_cell{index}{i}(1),data_cell{index}{i}(k), color_black, 'MarkerSize', 7);
            end
        end
    end

    figure(hObject.Parent);
    axes(hObject.Parent.UserData{1});
    plot(wavelength_black, eigs_black, color_black, 'MarkerSize', 7);
    hold on;
    plot(wavelength_red, eigs_red, color_red, 'MarkerSize', 7);
    hold off;
    axis([0 0.51 -0.05 0.31]);
    
    load('load_data.mat');
    axes(hObject.Parent.UserData{2});
    plot(load_data(:, 1), load_data(:, 3), 'k', 'LineWidth', 1.3);
    hold on;
    plot(load_data(index*step_length, 1), load_data(index*step_length, 3), 'r.', 'MarkerSize', 10);
    hold off;
    
    
    
    

end