function [variables, expression] = parsefilterstring(inputstr)

pattern = '\$[\w\d_.]+([^\w\d_]|$)';
variables = regexp(inputstr,pattern,'match');
if isempty(variables)
    error(['Not a valid filter:   ',inputstr]);
end
for i = 1:length(variables)
    variables{i} = regexprep(variables{i},'[^\w\d_.]','');
end
expression = regexprep(inputstr,'\$','structVar.'); %replace all '$' with 'structVar.'
    



% function [variable, expression] = parsefilterstring(inputstr)
% %
% % [variable, expression] = parsefilterstring(inputstr)
% % Parses an input string with the format 'variable: expression'
% 
% colon = strfind(inputstr,':');
% if isempty(colon)
%     parseError(inputstr);
% end
% colon = colon(1);
% 
% variable = inputstr(1:colon-1);
% expression = inputstr(colon+1:end);





