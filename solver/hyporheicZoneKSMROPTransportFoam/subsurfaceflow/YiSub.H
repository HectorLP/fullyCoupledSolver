{
    //- update dispersion tensor coefficient
    //DeffModel->correct(C_Sub,Ub,Sb);

    forAll(solutionSpecies, i)
    {
        volScalarField& YiSub = Y_Sub[i];
        Info << "Calculate effective diffusion coefficient." << endl;
        volScalarField DeffYi = Sb * D_sub;
        Info << "Calculate the ADE matrix and solve it." << endl;
        while(simpleSubsurface.correctNonOrthogonal())
        {
            fvScalarMatrix YiSubEqn
            (
                Sb * fvm::ddt(YiSub)
                + YiSub * fvc::ddt(Sb)
                - fvm::laplacian(DeffYi, YiSub, "laplacian(DeffYi,YiSub)")
                + fvm::div(phiSub, YiSub, "div(phiSub,YiSub)")
                ==
                fvOptions(YiSub)
            );
        
            YiSubEqn.relax();
            // // fvOptions.constrain(CSubEqn);
            YiSubEqn.solve();
            fvOptions.correct(YiSub);
            //C_Sub.correctBoundaryConditions();
            YiSub.max(0.0);
        }
        Info << "Concentration Min(YiSub) = " << min(YiSub).value() << " Max(YiSub)" << \
            max(YiSub).value() << endl;
    }
    Info << "Solve the reactions." << endl;
    rmSub->reactionStep(runTime.deltaT());
}