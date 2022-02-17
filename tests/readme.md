The testing here is a crude answer-testing framework. The principle is that you first run a test and some result (the answer) is saved. When you change code, you run the same test and compare the output to the existing answer. If you do not get any errors then your changes have not broken anything!


## storing answers 

So the run tests, you first have to store all the answers. To do that, first check out the `main` branch of the repository. Then in matlab or octave from the top-level of this repo, run:

```
addpath('./tests')
store_all_tests()
```

This should run all the tests and store answers in `tests/testdata/`. 

## running tests

once answers are stored, you can run all tests with

```
addpath('./tests')
test_all()
```

If you run `test_all()` immediately after `store_all_tests()`, you should not see any output. 

Once you have answers stored, you can re-run the tests to make sure the code still runs as expected. 

## adding tests 

To add a test:

1. create a new file starting with `test_*.m` in `tests/`. The test function should accept a single integer flag for storing/comparing and then should run code and  either save or compare output. See `tests/test_make_vm.m` for an example. 
2. add a call to the new function in `test_all` and `store_all`, switching the storage flag between 0 (for comparing) and 1 (for storing).
3. To store the answer of the new test without storing all the tests again, just run your new function with the storage flag set to 1. 
