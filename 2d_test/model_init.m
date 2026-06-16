function model_init()
    %clear
    setting_file_name = "SAC_est_2d_test_setting.json";
    mdlDir = fileparts(get_param(bdroot,"FileName"));
    p = jsondecode(fileread(fullfile(mdlDir, setting_file_name)));

    set_json_params(p)
    % 必要であればjson設定ファイルから str2func 関数を使って推定用の無名関数を作るようにする
    % ex. `f = str2func("@(x) x.^2 + 1");`

    %modelparamset()

end

function set_json_params(p)
    fn = fieldnames(p);
    for i = 1:numel(fn)
        k = fn{i};
        if ~isvarname(k)
            error("Invalid variable name in JSON key: %s", k);
        end

        v = p.(k);
        if isnumeric(v), v = double(v); end
        assignin("base", k, v);
    end

end
