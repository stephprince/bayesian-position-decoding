function f = setfilterfunction(f, funcname, loadvariables, varargin)
%f = setfilteriterator(f, funcname, loadvariables, options)
%Sets the run function for the filter.  This function is called by the
%filter's iterator.

for i = 1:length(f)
    f(i).function.name = funcname;
    if (nargin > 2)
        f(i).function.loadvariables = loadvariables;
    else
        f(i).function.loadvariables = {};
    end
    if (nargin > 3)
        f(i).function.options = varargin;
    else
        f(i).function.options = {};
    end
    f(i).output = [];
end