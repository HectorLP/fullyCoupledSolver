{
    //Calculate the compression coefficient
    //Direction of interfacial flux
    surfaceScalarField fluxDir = fvc::snGrad(alpha1) * meshSurface.magSf();
    //upwind and downwind alpha1
    surfaceScalarField alphaUp = upwind<scalar>(meshSurface, fluxDir).interpolate(alpha1);
    surfaceScalarField alphaDown = downwind<scalar>(meshSurface, fluxDir).interpolate(alpha1);
    //upwind and downwind C
    surfaceScalarField CSurfUp = upwind<scalar>(meshSurface, fluxDir).interpolate(C_Surf);
    surfaceScalarField CSurfDown = downwind<scalar>(meshSurface, fluxDir).interpolate(C_Surf);

    dimensionedScalar sgn = sign(max(alphaDown * CSurfDown) - max((1. - alphaUp) * CSurfUp));

    //normal compression coefficient 
    surfaceScalarField deltaCSurf1 = max(-max(C_Surf), min(max(C_Surf), (CSurfDown - CSurfUp) / \
        (alphaDown - alphaUp + 1e-4)));

    surfaceScalarField deltaCSurf2 = max(-max(C_Surf), min(max(C_Surf), CSurfDown / (alphaDown + \
        (1. - alphaDown) * Hc) - Hc * CSurfUp / (alphaUp + (1 - alphaUp) * Hc)));
    surfaceScalarField compressionCoeff = sgn * max(mag(deltaCSurf1), mag(deltaCSurf2));

    word alpharScheme("div(phirb, alpha)");
    word CSurfScheme("div(phi, C_Surf)");

    surfaceScalarField phic(1.0 * mag(phi/meshSurface.magSf()));
    //phic = min(0.001*phic, max(phic));
    surfaceScalarField phir = phic * mixture.nHatf();

    volScalarField& tmpCSurf = C_Surf;
    scalar maxCSurf = max(gMax(tmpCSurf), gMax(tmpCSurf.boundaryField())) + 1e-30;
    tmpCSurf.oldTime() == tmpCSurf.oldTime()/maxCSurf;
    tmpCSurf == tmpCSurf / maxCSurf;
    surfaceScalarField phiCSurfComp = fvc::flux
    (
        -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
        alpha1,
        alpharScheme
    );

    // surfaceScalarField phiCSurf = 
    // fvc::flux
    // (
    //     phi,
    //     tmpCSurf,
    //     CSurfScheme
    // )
    // +phiCSurfComp * compressionCoeff;

    // MULES::explicitSolve(tmpCSurf, phi, phiCSurf, 1.0, 0);
    tmp<surfaceScalarField> tCSurfPhiUn
    (
        fvc::flux
        (
            phi,
            tmpCSurf,
            CSurfScheme
        )
        + phiCSurfComp * compressionCoeff
    );

    {
        surfaceScalarField tCSurfPhi10 = tCSurfPhiUn;
        
        MULES::explicitSolve
        (
            geometricOneField(),
            tmpCSurf,
            phi,
            tCSurfPhi10,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField()
        );
    }

    tmpCSurf.oldTime() == tmpCSurf.oldTime() * maxCSurf;
    tmpCSurf == tmpCSurf * maxCSurf;
}