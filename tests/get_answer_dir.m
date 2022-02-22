function adir = get_answer_dir()
    adir = './tests/testdata/';
    if exist(adir) ~= 7
        mkdir(adir)
    end
end 
