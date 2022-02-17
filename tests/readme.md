The testing here is a "simple" answer-testing framework. The principle is that you first run a test and some result (the answer) is saved. When you edit any code, you run the same test and compare the output to the existing answer. If you do not get any errors then your changes have not broken anything!


## storing answers 

Before running tests, you first have to store all the answers based on code you trust. To do that, first check out the `main` branch of the repository. Then in matlab or octave from the top-level of this repo, run:

```
addpath('./tests')
store_all()ave
```

This should run all the tests and store answers as `.mat` files in `tests/testdata/`. 

## running tests

once answers are stored, you can run all tests with

```
addpath('./tests')
test_all()
```

If you run `test_all()` immediately after `store_all()`, all tests should pass because no code has changed!

Once you have answers stored, you can re-run the tests to make sure the code still runs as expected. 

To run (or store!) an individual test, you can use
```
addpath('./tests')
store_test = 0;  
call_a_test('test_name', store_test)
```

where `test_name` is the test name, and corresponds to a function file in `./tests`. The `store_test` variable is an integer flag -- if it is `0`, the test will compare to current output to the existing answer. If it is `1`, a new answer will be stored.  

## adding tests 

To add a test:

First create a new file starting with `test_*.m` in `tests/`. The test function should not accept any arguments and it needs to return two structure containing the values from running the code and the allowed tolerances for the test. For a new test function called `test_a_new_thing.m`, the first line would be
```
function [current_vals, diff_tolerances] = test_a_new_thing()
```
The form of these structures is  important and there are required fields -- see `tests/test_make_vm.m` for an example. 

After the new test is written, add the new test name to the `test_list` cell array in `get_full_test_list.m`. 

You can then use

```
addpath('./tests')
store_test = 1;  
call_a_test('test_a_new_thing', store_test)
``` 

to store the answer, after which

```
call_a_test('test_a_new_thing', store_test)
```

should run the comparison to the answer you just stored. 

At this point, the new test function will automatically be picked up by `test_all` and `store_all` as well.

## trouble shooting failing tests

A test can fail in at least two ways: the code may error or the comparison may fail. 

The testing framework does not attempt to catch errors, so if code changes cause an error, the test attempt also errors and you can debug as normal. If you are making changes to the arguments of a function, the test function will also need to be updated.

If a test fails a comparison, it means the modified code ran successfully but it generates numerical differences above the test's tolerance. This is not **necessarily** bad, but it should be investigated. If you are making changes for which you expect output to change, you can note that in a pull request and once the changes are merged, other developers can generate a new answer for the change.  

When a test fails a comparison, some information will be saved. For a failed test, two `.mat` files will be saved in `./tests/testdata`: one with a suffix `_bad.mat` and one with `_comparison.mat`. For example, if `test_make_vm` fails, you could load up the two files with:

```
current_vals = load('./tests/testdata/test_make_vm_bad.mat');
results = load('./tests/testdata/test_make_vm_comparison.mat');
```

The `current_vals` structure will contain the data output for the current version of the code, the `results` structure will contain information on the comparison, including which arrays failed their checks. 