if exist('../expcode', 'dir')
    addpath('../expcode')
else
    error(['To run this script, the "Matrix exponential times' ...
        ' vector" code must be downloaded to the folder ../expcode. ' ...
        'Can be downloaded at: ' ...
        'https://www.mathworks.com/matlabcentral/fileexchange/'...
        '29576-matrix-exponential-times-a-vector'])
end