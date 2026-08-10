function myfprintf(fid,varargin)
    fprintf(varargin{:});
    fprintf(fid,varargin{:});
end