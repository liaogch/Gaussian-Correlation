%set_env
%set path environment so that MATLAB can find the needed .m files
%without changing current working location

parent='D:\liao\MATLAB\Gaussian Correlation\';
path(path, parent);
path(path, [parent 'figure']);
path(path, [parent 'data']);
path(path, [parent 'dataresult']);
path(path, [parent 'figure\fig1']);
path(path, [parent 'figure\fig2']);