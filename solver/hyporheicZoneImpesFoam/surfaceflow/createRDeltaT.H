bool LTS = fv::localEulerDdt::enabled(meshSurface);

tmp<volScalarField> trDeltaT;

if (LTS)
{
    Info<< "Using LTS" << endl;

    trDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTName,
                runTime.timeName(),
                meshSurface,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            meshSurface,
            dimensionedScalar("one", dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
}
