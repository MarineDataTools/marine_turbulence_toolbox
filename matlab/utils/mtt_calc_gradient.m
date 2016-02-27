function [dxdy] = mtt_calc_gradient(x,y,varargin)
%
% computes central difference gradient delta x / delta y
% input: x,y, vectors of the same length
% output dxdy: delta x/delta y (dxdy(1) and dxdy(end) == NaN)
%        
%
% Part of the marine turbulence toolbox:
% https://github.com/MarineDataTools/marine_turbulence_toolbox

  global mtt_verbosity
  if(~isempty(mtt_verbosity))
    verbosity = mtt_verbosity;
  else
    verbosity = 0;
  end

  for i=1:length(varargin)
    if(strcmpi(lower(varargin{i}),'verbosity'))
      verbosity = varargin{i + 1};
    end
  end

  if(verbosity == 3)
    mtt_message(' ',1)
  end  

    
  if(length(x) ~= length(y))
    error('length of x and y differ, aborting')
    return
  end
  
  dxdy       = NaN(length(x),1);
  ynew       = NaN(length(x),1);

  for i=2:length(x)-1
    dxdy(i,1) = (x(i+1) - x(i-1))/(y(i+1) - y(i-1));
  end

end
