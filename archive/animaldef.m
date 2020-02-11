function animalinfo = animaldef(animalname)

switch lower(animalname)
    case 'barack'
        animalinfo = {'Barack', '/data/asinger/Bar/', 'bar'};
    case 'calvin'
        animalinfo = {'Calvin', '/data/asinger/Cal/', 'cal'};
    case 'dwight'
        animalinfo = {'Dwight', '/data/asinger/Dwi/', 'dwi'};
   
    %mattias's animals also on /data26/
    case 'dudley' %Animal 1, first exposure to first track on recording day 1, track with walls
        animalinfo = {'Dudley', '/data/mkarlsso/Dud/', 'dud'};
    case 'miles' %Animal 2, first exposure to first track on recording day 1, track with walls
        animalinfo = {'Miles', '/data/mkarlsso/Mil/', 'mil'};
    case 'conley' %Animal 3, first exposure to first track on recording day 1
        animalinfo = {'Conley', '/data/mkarlsso/Con/', 'con'};
    case 'ten' %Animal 4, first exposure to first track on recording day 1
        animalinfo = {'Ten', '/data/mkarlsso/Ten/', 'ten'};
    case 'bond' %Animal 5, pretrained on track 1 for 5 days
        animalinfo = {'Bond', '/data/mkarlsso/Bon/', 'bon'};
    case 'frank'  %Animal 6, pretrained on track 1 for 5 days
        animalinfo = {'Frank', '/data/mkarlsso/Fra/', 'fra'}; %also on /data26/, prev on /data19a/
    case 'nine' %Animal 7, pretrained on track 1 for 5 days
        animalinfo = {'Nine', '/data/mkarlsso/Nin/', 'nin'};
    case 'alex' % coded tracks same,  %pretrained on track 1 for 5 days
        animalinfo = {'Alex', '/data/mkarlsso/Ale/', 'ale'};
   
    % Maggie's animals, also on /data13/ and /vortex/
    case 'five'  %CA1
        animalinfo = {'Five', '/data/mcarr/Fiv/', 'Fiv'};
    case 'eight' %CA1
        animalinfo = {'Eight', '/data/mcarr/Eig/', 'Eig'};
    case 'coriander' %CA3 and CA1
        animalinfo = {'Coriander', '/data/mcarr/Cor/', 'Cor'};
    case 'six' % maggie says bad behavior
        animalinfo = {'Six', '/data/mcarr/Six/', 'Six'};
    case 'seven' %maggie says good behavior but very far dorsal, so EEG bleed from DG
        animalinfo = {'Seven', '/data/mcarr/Sev/', 'Sev'};
    otherwise
        error(['Animal ',animalname, ' not defined.']);
end