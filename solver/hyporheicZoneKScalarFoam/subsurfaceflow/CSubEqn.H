{
    //- update dispersion tensor coefficient
    DeffModel->correct(C_Sub,Ub,Sb);

    volScalarField DeffC = Sb * D_sub;

    fvScalarMatrix CSubEqn
    (
        R * Sb * fvm::ddt(C_Sub)
        + R * C_Sub * fvc::ddt(Sb)
        + fvm::div(phiSub, C_Sub, "div(phiSub, C_Sub)")
        - fvm::laplacian(DeffC, C_Sub, "laplacian(DeffC, C_Sub)")
        ==
        - R * Sb * fvm::Sp(lambda, C_Sub)
    );
    CSubEqn.solve();

    Info << "Concentration Min(C_Sub) = " << min(C_Sub).value() << " Max(C_Sub)" << \
        max(C_Sub).value() << endl;
}