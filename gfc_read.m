%% Aliakbar Zarkoob, AKA "XIV"
%  Gmail: XIV.Aliakbar.Zarkoob@gmail.com
%  Telegram: @XIVAliakbar

function data = gfc_read(filename)

    fid = fopen(filename);
    line_num = 0;
    while true
        line_num = line_num+1;
        line = fgetl(fid);
        content = strsplit(line);
    
        if strcmp(content{1}, 'end_of_head')
            break
        end
    
        if strcmp(content{1}, 'product_type')
            data.product_type = content{2};
        end
    
        if strcmp(content{1}, 'modelname')
            data.modelname = content{2};
        end
    
        if strcmp(content{1}, 'earth_gravity_constant')
            data.GM = str2double(content{2});
        end
    
        if strcmp(content{1}, 'radius')
            data.R = str2double(content{2});
        end
    
        if strcmp(content{1}, 'max_degree')
            data.max_degree = str2double(content{2});
        end
    
        if strcmp(content{1}, 'errors')
            data.errors = content{2};
        end
    
        if strcmp(content{1}, 'norm')
            data.norm = content{2};
        end
    
        if strcmp(content{1}, 'tide_system')
            data.tide_system = content{2};
        end
    
        if strcmp(content{1}, 'key')
            content = strsplit(line, '  ');
            data.key = strrep(strtrim(content(2:end)), ' ', '_');
        end
    
    end
    
    main = readtable(filename, 'NumHeaderLines', line_num, 'ReadVariableNames', false, 'FileType', 'text');
    main = main(:, 2:end);
    main.Properties.VariableNames = data.key;
    data.main = main;
    fclose(fid);

    C_lm = zeros(data.max_degree+1, data.max_degree+1);
    S_lm = zeros(data.max_degree+1, data.max_degree+1);
    for i = 0:max(data.max_degree)
        idx = data.main.L == i;
        C_lm(i+1, 1:sum(idx)) = data.main.C(idx);
        S_lm(i+1, 1:sum(idx)) = data.main.S(idx);
    end
    data.C_lm = C_lm;
    data.S_lm = S_lm;

end