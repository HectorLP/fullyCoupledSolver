{
    // Info << "Update the boundary condition with Ua" << endl;
    Info << "Update the boundary condition with Ub." << endl;
    Ub.correctBoundaryConditions();
    Info << "Update path field, depending on Ua or Ub." << endl;
    Ua.correctBoundaryConditions();

    Info << "Update path field, depending on Ua or Ub." << endl;
    forAll(meshSubsurface.boundary(),patchi)
    {
        if (isA< fixedValueFvPatchField<vector> >(Ua.boundaryField()[patchi]))
        {
            phia.boundaryFieldRef()[patchi] = Ua.boundaryField()[patchi] & meshSubsurface.Sf().boundaryField()[patchi];
        }
        if (isA< fixedValueFvPatchField<vector> >(Ub.boundaryField()[patchi]))
        {
            phib.boundaryFieldRef()[patchi] = Ub.boundaryField()[patchi] & meshSubsurface.Sf().boundaryField()[patchi];
        }
    }
    Info << "Solve the saturation equation." << endl;
    fvScalarMatrix SbEqn
        (
            eps*fvm::ddt(Sb) + fvc::div(phib) 
            ==
            // event source terms
            - sourceTerm
        );
    Info << "Finish constructing equation for solving." << endl;
    SbEqn.solve();
    // int ncells = meshSubsurface.cells().size();
    // for (int i = 0; i < ncells; ++i)
    // {
    //     if (Sb[i] < 0)
    //     {
    //         Sb[i] = 0;
    //     }
    // }
    Info << "Saturation b " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;

}
