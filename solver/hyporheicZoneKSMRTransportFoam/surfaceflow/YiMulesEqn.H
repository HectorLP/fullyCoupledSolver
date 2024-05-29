// {
// 	word alpharScheme("div(phirb,alpha)");
// 	word YiScheme("div(phi,Yi)");

//     // Standard face-flux compression coefficient
//     surfaceScalarField phic(speciesMixture.cYi()*mag(phi/meshSurface.magSf()));

//     surfaceScalarField phir(phic*mixture.nHatf());

//     const speciesTable& solutionSpecies = speciesMixture.species();

//     forAll(solutionSpecies, i)
// 	{
//         volScalarField& Yi = speciesMixture.Y(i);

//         scalar maxYi = max(gMax(Yi),gMax(Yi.boundaryField()));

//         Yi.oldTime() == Yi.oldTime()/maxYi;
//         Yi == Yi/maxYi;

// 		surfaceScalarField phiComp = fvc::flux
//         (
//             -fvc::flux(-phir, alpha2, alpharScheme),
//             alpha1,
//             alpharScheme
//         );

//         tmp<surfaceScalarField> tYiPhi1Un
//         (
//             fvc::flux
//             (
//                 phi,
//                 Yi,
//                 YiScheme
//             )
// 		+   phiComp*speciesMixture.compressionCoeff(i)
//         );

//         {
//             surfaceScalarField YiPhi10 = tYiPhi1Un;

//             MULES::explicitSolve
//             (
//                 geometricOneField(),
//                 Yi,
//                 phi,
//                 YiPhi10,
//                 zeroField(),
//                 zeroField(),
//                 oneField(),
//                 zeroField()
//             );
//         }

//         Yi.oldTime() == Yi.oldTime()*maxYi;
//         Yi == Yi*maxYi;
//     }
//  }
{
    //Calculate the compression coefficient
    //Direction of interfacial flux
    surfaceScalarField fluxDirY = fvc::snGrad(alpha1) * meshSurface.magSf();
    //upwind and downwind alpha1
    surfaceScalarField alphaUpY = upwind<scalar>(meshSurface, fluxDirY).interpolate(alpha1);
    surfaceScalarField alphaDownY = downwind<scalar>(meshSurface, fluxDirY).interpolate(alpha1);
    surfaceScalarField phic(1.0 * mag(phi/meshSurface.magSf()));
    //phic = min(0.001*phic, max(phic));
    surfaceScalarField phir = phic * mixture.nHatf();
    //upwind and downwind C
    forAll(solutionSpeciesSurf, i)
    {
        volScalarField& YiSurf = YSurf[i];
        Info << "The chemical is " << YiSurf.name() << endl;
        surfaceScalarField YiSurfUp = upwind<scalar>(meshSurface, fluxDirY).interpolate(YiSurf);
        surfaceScalarField YiSurfDown = downwind<scalar>(meshSurface, fluxDirY).interpolate(YiSurf);

        dimensionedScalar sgn = sign(max(alphaDownY * YiSurfDown) - max((1. - alphaUpY) * YiSurfUp));

        //normal compression coefficient 
        surfaceScalarField deltaYiSurf1 = max(-max(YiSurf), min(max(YiSurf), (YiSurfDown - YiSurfUp) / \
            (alphaDownY - alphaUpY + 1e-4)));

        surfaceScalarField deltaYiSurf2 = max(-max(YiSurf), min(max(YiSurf), YiSurfDown / (alphaDownY + \
            (1. - alphaDownY) * Hc) - Hc * YiSurfUp / (alphaUpY + (1 - alphaUpY) * Hc)));
        surfaceScalarField compressionCoeffYi = sgn * max(mag(deltaYiSurf1), mag(deltaYiSurf2));

        word alpharScheme("div(phirb, alpha)");
        word YiSurfScheme("div(phi, YiSurf)");



        volScalarField& tmpYiSurf = YiSurf;
        scalar maxYiSurf = max(gMax(tmpYiSurf), gMax(tmpYiSurf.boundaryField())) + 1e-30;
        tmpYiSurf.oldTime() == tmpYiSurf.oldTime()/maxYiSurf;
        tmpYiSurf == tmpYiSurf / maxYiSurf;
        surfaceScalarField phiYiSurfComp = fvc::flux
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
        tmp<surfaceScalarField> tYiSurfPhiUn
        (
            fvc::flux
            (
                phi,
                tmpYiSurf,
                YiSurfScheme
            )
            + phiYiSurfComp * compressionCoeffYi
        );

        {
            surfaceScalarField tYiSurfPhi10 = tYiSurfPhiUn;
            
            MULES::explicitSolve
            (
                geometricOneField(),
                tmpYiSurf,
                phi,
                tYiSurfPhi10,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
        }

        tmpYiSurf.oldTime() == tmpYiSurf.oldTime() * maxYiSurf;
        tmpYiSurf == tmpYiSurf * maxYiSurf;
    }
}