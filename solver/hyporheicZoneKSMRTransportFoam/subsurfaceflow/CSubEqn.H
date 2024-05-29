{
    //- update dispersion tensor coefficient
    //DeffModel->correct(C_Sub,Ub,Sb);
    Info << "Calculate effective diffusion coefficient." << endl;
    volScalarField DeffC = Sb * D_sub;
    Info << "Calculate the ADE matrix and solve it." << endl;
    while(simpleSubsurface.correctNonOrthogonal())
    {
        fvScalarMatrix CSubEqn
        (
            Sb * fvm::ddt(C_Sub)
            + C_Sub * fvc::ddt(Sb)
            - fvm::laplacian(DeffC, C_Sub, "laplacian(DeffC, C_Sub)")
            + fvm::div(phiSub, C_Sub, "div(phiSub,C_Sub)")
            ==
            fvOptions(C_Sub)
        );
    
        CSubEqn.relax();
        // // fvOptions.constrain(CSubEqn);
        CSubEqn.solve();
        fvOptions.correct(C_Sub);
        //C_Sub.correctBoundaryConditions();
        C_Sub.max(0.0);
    }
    Info << "Concentration Min(C_Sub) = " << min(C_Sub).value() << " Max(C_Sub)" << \
        max(C_Sub).value() << endl;
}