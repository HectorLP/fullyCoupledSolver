/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Description
    Timestep for subsurface flow solver 

\*---------------------------------------------------------------------------*/
if (adjustTimeStep)
{
    scalarField dSdt = mag(-fvc::div(phib))().internalField() / eps.internalField();
    scalar dSdTmax = gMax(dSdTmax + SMALL);

    //-Choose between Coats and Courant number as driving number
    scalar CFLUse = -1;
    scalar maxDeltaTFactSubFlow = -1;
    switch(static_cast<int>(CFLint))
    {
        case 0:
            #include "CoatsNo.H"
            break;
        // case 1:
        //     #include "ToddNo.H"
        //     break;
        case 1:
            #include "CourantNoSub.H"
            break;
    }
    scalar deltaTFactSubFlow = min(min(maxDeltaTFactSubFlow, 1.0 + 0.1 * maxDeltaTFactSubFlow), 1.2);

    scalar optionalDeltaTSubflow = (
            min(
                dtForS, 
                min(
                    maxDeltaTFactSubFlow * runTime.deltaValue(),
                    maxDeltaTFactSubFlow
                )
            )    
    );

    Info << "Delta t from the subsurface flow is " << optionalDeltaTSubflow << endl;
}
