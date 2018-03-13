%%
clear all
clc

% run tests in all the subdirectories

[passed,failed,incomplete] = find_and_run_tests(pwd);

display(['Total passed tests = ', num2str(passed)]);
display(['Total failed tests = ', num2str(failed)]);
display(['Total incomplete tests = ', num2str(incomplete)]);