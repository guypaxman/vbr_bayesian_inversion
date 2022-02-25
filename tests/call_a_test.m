function test_result = call_a_test(test_name, store_answer)
    %%%
    % call the test function named in test_name
    %%
    
    answer_file = [test_name, '.mat'];  % the expected results to load or store
    bad_file = [test_name, '_bad.mat'];  % only saved if store_answer == 0 and test fails
    test_result_file = [test_name, '_comparison.mat'];
    test_result = struct();
    [current_vals, diff_tolerances] = feval(test_name);  % call the test function
      
    if store_answer==1
        store_answer_file(answer_file, current_vals)
        test_result.passed = 1;
        test_result.message = ["succesfully stored answer for ", test_name];
    else
        % reload the stored result
        expected = load_answer(answer_file);

        % calculate differences
        result = compare_answers(expected, current_vals, diff_tolerances);

        if result.passed
            msg = [test_name, ' passed '];
            [time_diff, msg] = get_time_diff(expected, current_vals, msg);
            test_result.message = msg;
            test_result.passed = 1;
        else
            test_result.message = ["test failed, check ", test_result_file, " and ", bad_file];
            test_result.passed = 0;
            store_answer_file(bad_file, current_vals)
            store_answer_file(test_result_file, result)
        end
    end
end

function [time_diff, msg] = get_time_diff(expected_vals, current_vals, msg)
    % only called if tests pass
    time_diff = struct();
    time_diff.absolute = expected_vals.elapsed_time - current_vals.elapsed_time;
    frac = time_diff.absolute / expected_vals.elapsed_time;
    time_diff.percent = frac * 100;

    if time_diff.absolute < 0
        slow_fast = 'slower';
    else
        slow_fast = 'faster';
    end
    msg = [msg, ' in ', num2str(current_vals.elapsed_time), ' s '];
    msg = [msg, '(which is ',num2str(abs(time_diff.absolute)), ' s  or '];
    msg = [msg, num2str(abs(time_diff.percent)), ' pct ', slow_fast, ')'];
end

function result = compare_answers(expected, current_vals, diff_tolerances)

    % first check both structs have the same fieldnames
    result = struct();
    result.passed = 1;
    failed = 0;
    expected_fields = fieldnames(expected);
    for ifield = 1:numel(expected_fields)
        fld = expected_fields{ifield};
        if ~isfield(current_vals, fld)
            failed = 1;
            result.passed = 0;
            result.message = "expected and current structures have different fields.";
        end
    end


    % compare each of the arrays noted in the arrays_to_compare cell array
    failed_array_check = 0;
    if failed == 0
        result.array_diffs = struct();

        if isfield(expected, "arrays_to_compare")
            for i_array = 1:numel(expected.arrays_to_compare)
                [expected_val, current_val, final_name] = get_arrays(expected, current_vals, i_array);
                % calculate differences betweeen arrays
                diff = abs(expected_val(:) - current_val(:));
                diff_frac = diff ./ expected_val(:);
                max_abs = max(diff);
                max_frac = max(diff_frac(expected_val~=0));

                % store and flag those results
                result.array_diffs.(final_name) = struct();
                result.array_diffs.(final_name).max_abs_diff = max_abs;
                result.array_diffs.(final_name).max_frac_diff = max_frac;

                failed_test = 0;
                if max_abs > diff_tolerances.max_abs_diff
                    fail_type = "absolute difference";
                    failed_test = 1;
                elseif max_frac > diff_tolerances.max_frac_diff
                    fail_type = "fractional difference";
                    failed_test = 1;
                end
                if failed_test == 1
                    failed_array_check = 1;  % at least one test has failed!
                    result.array_diffs.(final_name).failed = 1;
                    result.array_diffs.(final_name).fail_type = fail_type;
                end

            end
        end
    end

    if failed_array_check
        result.passed = 0;
        result.message = "failed array difference check";
    else
        result.passed = 1;
        result.message = "passed all array difference checks";
    end
end

function [expected_val, current_val, final_name] = get_arrays(expected, current_vals, i_array)

    % access the nested fields and pull out arrays
    array_struct_fields = expected.arrays_to_compare{i_array};

    expected_val = expected;
    current_val = current_vals;

    final_name = '';
    for i_struct = 1:numel(array_struct_fields)
        current_field = array_struct_fields{i_struct};
        expected_val = expected_val.(current_field);
        current_val = current_val.(current_field);
        if i_struct == 1
            final_name = current_field;
        else
            final_name=[final_name, '_', current_field];
        end
    end
end
