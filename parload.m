function [loadedfile] = parload(filename)
% load the saved .mat file during parfor loop
  eval(['load ' filename]);
  loadedfile = x;
end

