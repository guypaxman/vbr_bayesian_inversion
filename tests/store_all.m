function store_all()

    list_of_tests = get_full_test_list();
    failed_tests = {};
    ifailures = 0;
    test_results = struct();
    store_answer = 1;

    for itest = 1:numel(list_of_tests)
        test_name = list_of_tests{itest};
        test_results.(test_name) = call_a_test(test_name, store_answer);
        if test_results.(test_name).passed == 0
            % if something went wrong on storing, probably will error before this
            ifailures = ifailures + 1;
            failed_tests{ifailures} = test_name;
        end
    end

    if ifailures == 0
        disp("all tests stored succesfully")
    else
        disp("the following tests did not store answers")
        for ifailure = 1:ifailures
            disp(strcat('    ', failed_tests{ifailure}))
        end
    end
end
