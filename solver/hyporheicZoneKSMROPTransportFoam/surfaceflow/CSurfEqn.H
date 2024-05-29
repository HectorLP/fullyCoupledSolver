{
    surfaceScalarField phiH = D_surf * (1 - Hc) / fvc::interpolate(alpha1 + Hc * \
        (1 - alpha1)) * fvc::snGrad(alpha1) * meshSurface.magSf();
    
    surfaceScalarField fluxDir2 = fvc::snGrad(alpha1) * meshSurface.magSf();
    surfaceScalarField alphaUp2 = upwind<scalar>(meshSurface, fluxDir2).interpolate(alpha1);
    surfaceScalarField phiHUp = D_surf / (alphaUp2 + (1 - alphaUp2) * Hc) * fvc::snGrad(alpha1) * \
        meshSurface.magSf();

    surfaceScalarField alphaDown2 = downwind<scalar>(meshSurface, fluxDir2).interpolate(alpha1);
    surfaceScalarField phiHDown = -Hc * D_surf / (alphaDown2 + (1 - alphaDown2) * Hc) * \
        fvc::snGrad(alpha1) * meshSurface.magSf();

    volScalarField SpC = (C_Surf - C_Surf.oldTime()) / runTime.deltaT();

    while (simple.correctNonOrthogonal())
    {
        //compoute C equation
        fvScalarMatrix CSurfEqn
        (
            fvm::ddt(C_Surf)
            - fvm::laplacian(D_surf, C_Surf)
            ==
            SpC
        );

        if (phiHScheme == "Gauss linear")
        {
            CSurfEqn += fvm::div(phiH, C_Surf, "div(phiH, C_Surf)");
        }
        else if (phiHScheme == "Gauss upwind")
        {
            CSurfEqn += fvm::div(phiHUp, C_Surf, "div(phiH, C_Surf)")
            + fvm::div(phiHDown, C_Surf, "div(phiH, C_Surf)");
        }
        else 
        {
            Info << "div(phiH, C_Surf) should be equal to Gauss linear or Gauss upwind " << endl 
                << abort(FatalError); 
        }

        if (simple.finalNonOrthogonalIter())
        {
            CSurfEqn.solve(meshSurface.solver("CSurfFinal"));
        }
        else 
        {
            CSurfEqn.solve(meshSurface.solver("C_Surf"));
        }
    }
    Info << "The tracer concentration is " << C_Surf.weightedAverage(meshSurface.V()).value() 
        << " Min(C_Surf) = " << gMin(C_Surf.internalField()) << " Max(Yi) = " << gMax(C_Surf.internalField())
        << endl;
    
}