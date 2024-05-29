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

    forAll(solutionSpeciesSurf, i)
    {
        volScalarField& YiSurf = YSurf[i];
        //compute diffusion coefficient


        volScalarField SpYiSurf = (YiSurf - YiSurf.oldTime()) / runTime.deltaT();

        while (simple.correctNonOrthogonal())
        {
            //compoute C equation
            fvScalarMatrix YiSurfEqn
            (
                fvm::ddt(YiSurf)
                - fvm::laplacian(D_surf, YiSurf)
                ==
                SpYiSurf
            );

            if (phiHScheme == "Gauss linear")
            {
                YiSurfEqn += fvm::div(phiH, YiSurf, "div(phiH, YiSurf)");
            }
            else if (phiHScheme == "Gauss upwind")
            {
                YiSurfEqn += fvm::div(phiHUp, YiSurf, "div(phiH, YiSurf)")
                + fvm::div(phiHDown, C_Surf, "div(phiH, YiSurf)");
            }
            else 
            {
                Info << "div(phiH, C_Surf) should be equal to Gauss linear or Gauss upwind " << endl 
                    << abort(FatalError); 
            }

            if (simple.finalNonOrthogonalIter())
            {
                YiSurfEqn.solve(meshSurface.solver("YiSurfFinal"));
            }
            else 
            {
                YiSurfEqn.solve(meshSurface.solver("YiSurf"));
            }
        }
        Info << "The tracer concentration is " << YiSurf.weightedAverage(meshSurface.V()).value() 
            << " Min(YiSurf) = " << gMin(YiSurf.internalField()) << " Max(Yi) = " << gMax(YiSurf.internalField())
            << endl;
    }
 }