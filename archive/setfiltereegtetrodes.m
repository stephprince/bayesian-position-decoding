function f = setfiltereegtetrodes(f, filterString)
% f = setfiltereegtetrodes(f, filterString)
% For each epoch in the filter F, this function finds the indices to the
% tetrodesthat satisfy the given filter condtions in filterString.  
% filterstring can either be a single string, or a cell array of strings.  
%
% If it is a single string or one element cell array, the filter operates on
% the tetinfo structure and the syntax for filterString is defined in 
% EVALUATEFILTER.m. The animal and desired epochs for the filter need to be 
% predefined by createfilter calls.
% 
% If filterstring is three element cell array, the first element of 
% the cell array is taken as a the name of a function that returns a single
% number corresponding to an eegtetrode (e.g. 'geteegtet').  The second element
% is the type of eeg data desired (e.g. 'theta'), and the third
% element is the list of arguments to the function.  In this case the function
% is called once for each cell in f.an(#).data, producing an eegtetrode
% associated with each cell in the data element of f.


for an = 1:length(f)
    if isempty(f(an).animal)
        error(['You must define an animal for the filter before filtering the tetrodes'])
    end
    if isempty(f(an).epochs)
        error(['You must define the desired epochs for the filter before filtering the tetrodes'])
    end
    datadir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    tetinfo = loaddatastruct(datadir,animalprefix,'tetinfo');
    adir = f(an).animal{2};
    aprefix = f(an).animal{3};

    for g = 1:length(f(an).epochs)
        if isempty(f(an).epochs{g})
            f(an).eegdata{g} = [];
        end
        for e = 1:size(f(an).epochs{g},1)
	    if (iscell(filterString) & (length(filterString) > 1))
	        for c = 1:size(f(an).data{g}{e},1)
		    % this means that the first element of the filter is the 
		    % name of a function, so we need to run the function
		    cind = [f(an).epochs{g}(e,:) f(an).data{g}{e}(c,:)];
		    f(an).eegdata{g}{e}(c) = feval(filterString{1}, adir, ...
			aprefix, cind, filterString{2}, filterString{3:end}); 
		end
	    else
	        % this is a standard filter
	        f(an).eegdata{g}{e} = evaluatefilter( ...
  tetinfo{f(an).epochs{g}(e,1)}{f(an).epochs{g}(e,2)}, filterString);
	    end
        end
    end
end
