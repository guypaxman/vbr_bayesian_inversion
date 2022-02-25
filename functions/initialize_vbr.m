function initialize_vbr()
    % ensure that the VBRc is in the working path
    vbr_path_file = './vbr_path.txt';

    if exist(vbr_path_file) == 2
        vbr_path = fileread(vbr_path_file);
    else
        % try the environment var
        vbr_path = getenv("VBRpath")
        if length(vbr_path) == 0
            display_no_path_error()
        end
    end

    if exist(vbr_path) ~= 7
        display_no_dir_error(vbr_path)
    end

    addpath(vbr_path)
    vbr_init();

end


function display_no_path_error()
    disp("Could not find vbr path. You can set the vbr path by one of the following:")
    disp("    + in matlab, call:")
    disp("        addpath('./functions); set_vbr_path('path/to/vbr')")
    disp("    + save the path in the top level of this repo in ./vbr_path.txt ")
    disp("    + use the VBRpath environment variable. For bash: ")
    disp("        $ export VBRpath=/path/to/vbr/")
    error("no route to the VBRc :( ")
end


function display_no_dir_error(vbr_path)
    disp(["The supplied vbr path, ", vbr_path, " does not exist\n"])
    display_no_path_error()
end
