/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    Info << "Choose the time step." << endl;
    scalar maxDeltaTFactSurf =
        min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL));
    Info << "Calculate deltaTFact." << endl;
    scalar deltaTFactSurf = min(min(maxDeltaTFactSurf, 1.0 + 0.1*maxDeltaTFactSurf), 1.2);
    Info << "Find min value. " << endl;
    // runTime.setDeltaT
    scalar optionalDeltaTSurfFlow =
    (
        min
        (
            deltaTFactSurf*runTime.deltaTValueSurf(),
            maxDeltaTFactSurf
        )
    );

    Info << "Delta t from the subsurface flow is " << optionalDeltaTSurfFlow << endl;
}

// ************************************************************************* //
