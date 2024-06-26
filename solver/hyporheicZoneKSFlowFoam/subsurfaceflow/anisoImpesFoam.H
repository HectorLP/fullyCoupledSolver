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

Application
    anisoImpesFoam

Description
    Transient solver for incompressible two-phase flow (Darcy's law) in porous media
    using the IMPES method (IMplicit Pressure Explicit Saturation).
    Permeability is anisotropic (K == volTensorField)

Developers
    P. Horgue, C. Soulaine, J. Franc, R. Guibert and G. Debenest
    "An open-source toolbox for multiphase flow in porous media"

\*---------------------------------------------------------------------------*/

// #include "fvCFD.H"
// #include "incompressiblePhase.H"
// #include "capillarityModel.H"
// #include "relativePermeabilityModel.H"
// #include "sourceEventFile.H"
// // #include "outputEventFile.H"
// #include "patchEventFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// using namespace Foam;

// int main(int argc, char *argv[])
// {
    // #include "setRootCase.H"
    // #include "createTime.H"
    // #include "createMesh.H"
    // #include "createTimeControls.H"
    // #include "readGravitationalAcceleration.H"
    // #include "createFields.H"
    // #include "createSbFields.H"
    // #include "readTimeControls.H"
    // #include "readEvent.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Info<< "\nStarting time loop\n" << endl;

    // while (runTime.run())
    // {
    
    if (outputEventIsPresent) outputEvent.update(runTime.timeOutputValue());
    if (sourceEventIsPresent) sourceEvent.update(runTime.timeOutputValue());
    if (patchEventIsPresent) patchEvent.update(runTime.timeOutputValue());
    // #include "setDeltaT.H"

    // runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;
    Info << "Computer Source term (not inlcuded.)" << endl;
    #include "computeSourceTerm.H"
    Info << "Calculate saturation update." << endl;
    //- Solve saturation equation (explicit)
    #include "SEqn.H"
    Info << "Update the hydraulic properties." << endl;
    #include "updateSbProperties.H"

    //- Solve pressure equation (implicit)
    Info << "Update pressure field from the new hydraulic properties." << endl;
    Info << "Start to calculate the pressure field." << endl;
    #include "pEqn.H"
    Info << "Calculate the tracer transport." << endl;
    // #include "CSubEqn.H"
    // #include "YiSub.H"
    // Info << "Output the updated fields and so on." << endl;
    // runTime.write();

    // #include "eventWrite.H"

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
//     }

//     Info<< "End\n" << endl;

//     return 0;
// }

// ************************************************************************* //
