function mtt_message(message,indent)

% Prints messages, including the name of the function calling
% mtt_message.
%
% INPUT ARGUMENTS:
% message:              string with message text
% indent:               indentation level (0-3)
%                       0 dont print anything
%                       3 indent message by two tabs
% 
% Part of the marine turbulence toolbox:
% https://github.com/MarineDataTools/marine_turbulence_toolbox  

%_____get name of calling function_________________________________________

  s      = dbstack;
  
  file   = s(end).file;
  name   = s(end).name;  


%_____compose message string_______________________________________________

  switch indent
    case 0
      return
    case 1
      str = sprintf('\n*** %s: %s',file,message);
    case 2
      str = sprintf('\t ** %s: %s',file,message);
    case 3
      str = sprintf('\t \t * %s: %s',file,message);
    otherwise
      error('Unkown indentation level');
  end
  
  
%_____print________________________________________________________________
  
  disp(str);

