function [varargout] = IQMstochsim2(model,V,time,varargin)
% IQMstochsim2: Stochastic simulation of IQMmodels that only contain mass
% action kinetics. The simulator is based on the paper:
% Ullah, M., Schmidt, H., Cho, K.-H., Wolkenhauer, O. (2006) Deterministic
% Modelling and Stochastic Simulation of Pathways using MATLAB, 
% IEE Proceedings - Systems Biology, 153(2), 53-60
%
% The IQMmodel needs to have a certain format, explained below.
%
% USAGE:
% ======
% [output] = IQMstochsim(model,V,time)         
% [output] = IQMstochsim(model,V,time,runs)         
% [output] = IQMstochsim(model,V,time,runs,Nsample)         
%
% model:   IQMmodel 
%          There are certain limitations on an IQMmodel used for stochastic
%          simulation. Please read further below.
% V:       Volume of the reaction space (given in Liter)
% time:    End time for simulation 
% runs:    number of realizations(simulation runs) (default: 1)
% Nsample: Each Nsample-th point will be used for output (to save memory) 
%          (default: 100)
%
% Output Arguments:
% =================
% If no output arguments are given, the result of the simulation is
% plotted (if runs>1 only the mean is plotted). Otherwise the output
% argument has the following structure:
%
% output.time:            cell-array with time vectors for the single runs 
% output.speciesdata:     cell-array with simulation data for the single runs
% output.runs:            number of runs
% output.timemean:        ensemble of all time instants in the single runs
% output.speciesdatamean: matrix containing the means of the simulation data
% output.species:         cell-array containing the names of the species
%
% FORMAT OF THE IQMmodel:
% ======================
% IQMmodels that can be used for stochastic simulation need to follow some
% rules:
% 1) All reaction kinetics need to be of mass action type and be defined in
%    the following syntax:     'ReactionName' = 'kineticParameter' * ...
% 2) All reactions have to be irreversible. You can use the 
%    function IQMmakeirreversible to convert your model
% 3) The reactions can at maximum have 2 substrates and 2 products.
% 4) The right hand side of the ODEs needs only to consist of reaction rate
%    terms and eventually stoichiometric coefficients. This is necessary in
%    order to be able to determine the stoichiometric matrix.
%    More information about the required syntax can be found in the help
%    text of the function IQMstoichiometry
% 5) No variables, functions, events, functionsMATLAB are allowed to be
%    present.
% 6) Initial conditions of species are assumed to be given in numbers of
%    molecules

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2008 Henning Schmidt, henning@IQMtoolbox2.org
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle variable input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runs = 1;
Nsample = 100;
shownr = 0;
if nargin == 3,
elseif nargin == 4,
    runs = varargin{1};
elseif nargin == 5,
    runs = varargin{1};
    Nsample = varargin{2};
elseif nargin == 6,
    runs = varargin{1};
    Nsample = varargin{2};
    shownr  = varargin{3};
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert model to MA structure
MA = IQMconvert2MA(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain necessary data for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial numbers of molecules
n0 = MA.initialConditions;
% stoichiometric matrix
D = MA.N;
% kinetic parameters
k = MA.kineticParameters';
% Avogadro's number times volume
NAV = V*1e-9*6.02214199e23;   
% Stoich coeffs of the reactants (obtained via IQMconvert2MA)
L = MA.L;               
% molecularity
K = sum(L);   
% 'particle' rate constant
kp = k./NAV.^(K-1);
% normalise for correct units
kp = kp.*(1e9).^(K-1);         
% stochastic rate constant 
c = kp.*prod(factorial(L));  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stochastic simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ts,Ns,TT,NBAR] = stochIQM(n0,c,D,L,time,runs,Nsample,shownr);
if nargout == 0,
    % Display mean data
    datanames = MA.species;
    IQMplot(TT,NBAR,datanames);
else
    % return results in structure
    output = [];
    output.time = Ts;
    output.speciesdata = Ns;
    output.runs = runs;
    output.timemean = TT;
    output.speciesdatamean = NBAR;
    output.species = MA.species;
    varargout{1} = output;
end
return

























