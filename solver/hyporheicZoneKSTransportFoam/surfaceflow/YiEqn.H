{
    const speciesTable& solutionSpecies = speciesMixture.species();

    forAll(solutionSpecies, i)
    {
        volScalarField& Yi = speciesMixture.Y(i);

        //compute diffusion coefficient
        surfaceScalarField D = speciesMixture.DmY(i);

        //calculate Henry's transfer flux
        //surfaceScalarField phiH = speciesMixture.phiH(i);
        surfaceScalarField phiHUp = speciesMixture.phiHUp(i);
        surfaceScalarField phiHDown = speciesMixture.phiHDown(i);

        volScalarField SpY = (Yi-Yi.oldTime())/runTime.deltaT();

        while (simple.correctNonOrthogonal())
        {
            //compute Y equation
            fvScalarMatrix YiEqn
            (
                fvm::ddt(Yi)
                - fvm::laplacian(D,Yi)
                //+ fvm::div(phiH,Yi,"div(phiH,Yi)")
                + fvm::div(phiHUp,Yi,"div(phiH,Yi)")
                + fvm::div(phiHDown,Yi,"div(phiH,Yi)")
                ==
                SpY
            );

            //solve equations
            if (simple.finalNonOrthogonalIter())
            {
                YiEqn.solve(meshSurface.solver("YiFinal"));
            }
            else
            {
                YiEqn.solve(meshSurface.solver("Yi"));
            }
        }

        Info<< "Species concentration = "
        << Yi.weightedAverage(meshSurface.V()).value()
        << "  Min(Yi) = " << gMin(Yi.internalField())
        << "  Max(Yi) = " << gMax(Yi.internalField())
        << endl;
    }

    //calculate mass transfer flux
    speciesMixture.correct();
 }