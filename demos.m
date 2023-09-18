function tbxStruct  = demos

% DEMOS Demo list for the Aircraft Library.

% Version 1.12 (R11.1,R18a)
% Giampiero Campa, 29-Feb-2004
% Copyright 2018 The MathWorks, Inc.

if nargout == 0, demo blockset; return; end

tbxStruct.Name='Aircraft Library';
tbxStruct.Type='Blockset';

tbxStruct.Help={
   'Airlib is a library of aircraft models to be '
   'used with Simulink 3.1 (Matlab 5.3) or later.'
};

 tbxStruct.DemoList = { ' Aircraft Library',  'airlib'; ...
                        ' Aircraft Models',  'airlibex'; ...
                        ' General Controller',  'airgk'; ...
                        ' Wind Handling',  'fdcwind'; ...
                        ' Boeing 747',  'b747cl'};
