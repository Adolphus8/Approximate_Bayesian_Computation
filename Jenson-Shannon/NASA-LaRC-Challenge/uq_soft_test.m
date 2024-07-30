% TESTING NASA UQ CHALLENGE SOFTWARE
%
% Execute this script and compare the truth results below with those generated on your computer.
%

% CODED BY:  Luis G. Crespo NIA/NASA-DSCB
% CODED ON:  4 Feb, 2013

clear all;
your_g=x_to_g(p_to_x([(1:21)/21;(21:-1:1)/21;ones(1,21)/2;ones(1,21);zeros(1,21)]));

dat=NaN;
switch mexext
    case 'mexmaci64'  % mac 64-bit
        if exist('g_mexmaci64')==2
            dat=g_mexmaci64;
        else
            disp('Code validation data has not been produced for this architecture')
        end
    case 'mexglx'     % linux 32-bit
        if exist('g_mexglx')==2
            dat=g_mexglx;
        else
            disp('Code validation data has not been produced for this architecture')
        end
    case 'mexa64'     % linux 64-bit
        if exist('g_mexa64')==2
            dat=g_mexa64;
        else
            disp('Code validation data has not been produced for this architecture')
        end
    case 'mexw32'     % Windows 32-bit
        if exist('g_mexw32')==2
            dat=g_mexw32;
        else
            disp('Code validation data has not been produced for this architecture')
        end
    case 'mexw64'     % Windows 64-bit
        if exist('g_mexw64')==2
            dat=g_mexw64;
        else
            disp('Code validation data has not been produced for this architecture')
        end
    otherwise
end

try
    maxrelerr = discrepancy(your_g, dat);
catch err
    display(err.message)
end


