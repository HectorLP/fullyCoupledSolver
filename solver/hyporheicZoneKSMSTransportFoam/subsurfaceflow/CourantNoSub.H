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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNuma = 0.0;
scalar meanCoNumb = 0.0;

if (meshSubsurface.nInternalFaces())
{
    scalarField sumPhia
        (
            fvc::surfaceSum(mag(phia))().internalField()
        );
    scalarField sumPhib
        (
            fvc::surfaceSum(mag(phib))().internalField()
        );

    CoNum = max(0.5*gMax(sumPhia/meshSubsurface.V().field())*runTime.deltaTValue(),0.5*gMax(sumPhib/meshSubsurface.V().field())*runTime.deltaTValue());

    meanCoNuma = 
        0.5*(gSum(sumPhia)/gSum(meshSubsurface.V().field()))*runTime.deltaTValue();
    meanCoNumb = 
        0.5*(gSum(sumPhib)/gSum(meshSubsurface.V().field()))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNuma << "," << meanCoNumb << " max: " << CoNum << endl;

CFLUse = CoNum;
maxDeltaTFact = maxCo/(CFLUse + SMALL);

// ************************************************************************* //
