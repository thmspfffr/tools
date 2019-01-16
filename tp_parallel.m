function y = tp_parallel(fn,outdir,start,override)
% Takes input:
% fn: filename (without ending)
% outdir: out directory
% start: 1 or 0
% output: 1 if file exists: continue, otherwise 0

if start == 1
  % start of processing
  if ~override
    if ~exist(sprintf([outdir '%s_processing.txt'],fn)) && ~exist(sprintf([outdir '%s.mat'],fn))
      system(['touch ' outdir sprintf('%s_processing.txt',fn)]);
      y = 0;
    else
      y = 1;
    end
  else
    if ~exist(sprintf([outdir '%s_processing.txt'],fn))
      system(['touch ' outdir sprintf('%s_processing.txt',fn)]);
      y = 0;
    else
      y = 1;
    end
  end
elseif start == 0
  % end
  if exist(sprintf([outdir '%s_processing.txt'],fn)) || exist(sprintf([outdir '%s.mat'],fn))
    delete(sprintf([outdir '%s_processing.txt'],fn))
  end
end
