function test_all()

    list_of_tests = get_full_test_list();
    failed_tests = {};
    ifailures = 0;
    test_results = struct();
    store_answer = 0;

    for itest = 1:numel(list_of_tests)
        test_name = list_of_tests{itest};
        disp(' ')
        disp(['Running ', test_name, ' ...'])
        test_results.(test_name) = call_a_test(test_name, store_answer);
        msg = test_results.(test_name).message;
        if test_results.(test_name).passed == 0
            disp('!!!! FAILURE !!!!')
            ifailures = ifailures + 1;
            failed_tests{ifailures} = test_name;
        end
        disp(['    ', msg])
    end

    disp(' ')
    disp('-------------')
    disp('Test Summary')
    disp('-------------')

    if ifailures == 0
        disp(' ')
        disp("All tests passed :D ")
    else
        disp("Detected at least one failure :'( ")
        disp("    the following tests failed: ")

        for ifailure = 1:ifailures
            disp(['        ', failed_tests{ifailure}])
        end
    end
end
