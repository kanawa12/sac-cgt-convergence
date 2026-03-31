function model_init(setting_file_name)
    mdlDir = fileparts(get_param(bdroot,"FileName"));
    p = jsondecode(fileread(fullfile(mdlDir, "params", setting_file_name)));

    set_json_params(p)

    % 出力先（モデル直下 params）
    genDir = fullfile(mdlDir, "params");

    % 生成（下で定義する関数を呼ぶ）
funclist = {
    "est_expr", "am,bm,kx,ku";
    "cgt_expr", "am,bm,ap,bp"
};

for i = 1:size(funclist, 1)
    key  = funclist{i, 1};
    args = funclist{i, 2};
    generate_function(genDir, string(p.(key)), key, args); % チェックはしない
end

    % パス追加＆反映
    if ~contains(path, genDir)
    addpath(genDir);
    rehash;
    end

    % 既にロードされている関数を捨てる（上書き反映のため）
    clear calc_y_generated;

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

function generate_function(genDir, expr, funcname, paramstr)
    % ---- 生成する関数（固定名）----
    fpath = fullfile(genDir, funcname + "_generated.m");
    fid = fopen(fpath, "w");
    assert(fid>0, "Cannot write: %s", fpath);

    fprintf(fid, "function y = calc_y_generated(" + paramstr + ")\n");
    fprintf(fid, "%% Do not edit by hand.\n\n");
    fprintf(fid, "y = %s;\n", expr);
    fprintf(fid, "end\n");
    fclose(fid);
end

