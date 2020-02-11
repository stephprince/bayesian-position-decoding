function f = runfilter(f)

for an = 1:length(f)
    iterator = f(an).iterator;
    f(an) = feval(iterator,f(an));
end
    
