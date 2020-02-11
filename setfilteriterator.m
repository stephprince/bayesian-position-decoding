function f = setfilteriterator(f, funchandle)
%f = setfilteriterator(f, funchandle)
%Sets the iterator function for the filter

for i = 1:length(f)
    f(i).iterator = funchandle;
end