
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
    
    // if (outputEventIsPresent) outputEvent.update(runTime.timeOutputValue());
    // if (sourceEventIsPresent) sourceEvent.update(runTime.timeOutputValue());
    // if (patchEventIsPresent) patchEvent.update(runTime.timeOutputValue());
    // #include "setDeltaT.H"

    // runTime++;

    Info<< "Time = " << runTime.timeName() << nl << endl;
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
    #include "CSubEqn.H"
    #include "YiSub.H"
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
