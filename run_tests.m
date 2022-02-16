% to test changes

addpath('./tests')
addpath('./functions')
addpath('./inv_functions')
addpath('./GLAD25')
addpath(getenv('VBRpath'))
vbr_init;
test_all()
