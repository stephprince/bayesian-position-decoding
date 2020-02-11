function lineardayprocess(directoryname,fileprefix,days, varargin)
%LINEARDAYPROCESS(directoryname,fileprefix,days, options)
%
%Runs linearizeposition for all run epochs in each day and saves the data in
%'linpos' in the directoryname folder.  See LINEARIZEPOSITION for the definitions 
%of the options.
%
%directoryname - example 'data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile
%
%days -          a vector of experiment day numbers 
%

lowercasethree = '';

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
        
    end
end


days = days(:)';

for day = days
   linpos = [];
   
   dsz = '';
   if (day < 10)
      dsz = '0';
   end

   eval(['load ',directoryname,fileprefix,'task',dsz,num2str(day), '.mat']);
   eval(['task = ',lowercasethree,'task;'])
   for i = 1:length(task{day})   
      if ((~isempty(task{day}{i})) && (strcmp(task{day}{i}.type,'run')) )
            
            disp(['Day ',num2str(day), ', Epoch ',num2str(i)])
            index = [day i];
            [linpos{day}{i}.statematrix,linpos{day}{i}.segmenttable, linpos{day}{i}.trajwells, linpos{day}{i}.wellSegmentInfo, linpos{day}{i}.segmentInfo] = linearizeposition(directoryname,fileprefix, index, varargin);
            linpos{day}{i}.statematrixfields = ['time startwell endwell segment segheaddir velocity lineardist'];
            linpos{day}{i}.segmenttablefields = ['trajnum segnum segmentID'];

      end
   end
   
   eval([lowercasethree,'linpos = linpos;']);
   eval(['save ',directoryname,fileprefix,'linpos',dsz,num2str(day),' ',lowercasethree,'linpos']);
end
