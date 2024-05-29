{

    fvScalarMatrix pSubEqn
        (
            fvm::laplacian(-Mf, pSub) + fvc::div(phiG)
            // capillary term
            + fvc::div(phiPc)*activateCapillarity
            ==
            // event source terms
            - sourceTerm
        );
    Info << "FInish constructing the matrix for solving the pressure field." << endl;
    pSubEqn.solve();
    Info << "Finish solving the pressure field and construct pressure flux." << endl;
    phiP = pSubEqn.flux();
    Info << "Get the pressure flux and reconstruct the entire flux." << endl;
    phiSub = phiP+phiG+phiPc*activateCapillarity;
    Info << "Get the whole flux and calculate the velocity." << endl;
    USub = fvc::reconstruct(phiSub);
    Info << "Update the velocity on the boundary." << endl;
    USub.correctBoundaryConditions();
    Info << "Update the phase flux a and b." << endl;
    phib == Fbf*phiSub+(((Lbf-Fbf*Lf)&g) & meshSubsurface.Sf())+(1-Fbf)*phiPc*activateCapillarity;
    Info << "Update the flux of phase a." << endl;
    phia == phiSub - phib;
    Info << "Caluclate the velocity of phase b." << endl;
    Ub = fvc::reconstruct(phib);
    Info << "Cacluate the velocity of phase a." << endl;
    Ua = USub-Ub;
    Info << "Update the velocity of phase a at the boundary." << endl;
    Ua.correctBoundaryConditions();
    Info << "Update the velocity of phase b at the boundary." << endl;
    Ub.correctBoundaryConditions();  

}
