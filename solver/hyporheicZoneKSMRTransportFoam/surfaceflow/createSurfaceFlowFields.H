#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface
);

#include "createSurfaceFlowPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

volScalarField gradalpha1
(
    IOobject
    (
        "gradalpha1",
        runTime.timeName(),
        meshSurface,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(fvc::grad(alpha1))
);

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

Info << "Copy the saturation field." << endl;
// volScalarField alpha1Tmp(mixture.alpha1());
volScalarField alpha1Tmp
(
    IOobject
    (
        "alpha1Tmp",
        runTime.timeName(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshSurface
);

// volScalarField CSurf //surface scalar 
// (
//     IOobject
//     (
//         "C",
//         runTime.timeName(),
//         meshSurface,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     meshSurface
// );

Info << "Finish copying the saturation field." << endl;
//Whether to solve the scalar C (non-reactive tracer) in surface water flow
// Switch CSurfTracerSwitch
// (
//     mixture.lookup("CSurftracerSwitch")
// );

// //Whether to solve reactive transport process in surface flow
// Switch ChemSurfSwitch
// (
//     mixture.lookup("CheSurfSwitch")
// );

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        meshSurface,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        meshSurface,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


Info<< "\nReading g for surface flow" << endl;
uniformDimensionedVectorField gSurf
(
    IOobject
    (
        "gSurf",
        runTime.constant(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


Info<< "\nReading hRef" << endl;
uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        runTime.constant(),
        meshSurface,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar("hRef", dimLength, 0)
);


Info<< "Calculating field g.h\n" << endl;
dimensionedScalar ghRefSurf
(
    mag(gSurf.value()) > SMALL
    ? gSurf & (cmptMag(gSurf.value())/mag(gSurf.value()))*hRef
    : dimensionedScalar("ghRef", gSurf.dimensions()*dimLength, 0)
);
volScalarField ghSurf("ghSurf", (gSurf & meshSurface.C()) - ghRefSurf);
surfaceScalarField ghfSurf("ghfSurf", (gSurf & meshSurface.Cf()) - ghRefSurf);
// #include "readGravitationalAcceleration.H"
// #include "readhRef.H"
// #include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        meshSurface,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*ghSurf
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p + rho*ghSurf;
}

meshSurface.setFluxRequired(p_rgh.name());
meshSurface.setFluxRequired(alpha1.name());

// PtrList<volScalarField>& surAlpha = alpha1;

// multivariateSurfaceInterpolationScheme<scalar>::fieldTable surfAlphaFields;

// surfAlphaFields.add(surAlpha);

// PtrList<volScalarField> alphaFromSub = surAlpha;

// surfAlphaFields.add(alphaFromSub);

// Info << "Reading the thermoPhysicalProperties.\n" << endl;

// IOdictionary thermoPhysicalProperties
// (
//     IOobject
//     (
//         "thermoPhysicalProperties",
//         runTime.constant(),
//         meshSurface,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     )
// );

// Info << "Create chemical species mixture.\n" << endl;
// twoPhaseMultiComponentTransportMixture<inertMultiComponentMixture, inertMultiComponentMixture> speciesMixture(meshSurface, alpha1);
volScalarField C_Surf
(
    IOobject
    (
        "C_Surf",
        runTime.timeName(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface
);

IOdictionary transportProperties_surf
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
dimensionedScalar D_surf(transportProperties_surf.lookup("D_surf"));
dimensionedScalar Hc(transportProperties_surf.lookup("Hc"));

//******************************************************************//
// scheme for phiH
//******************************************************************//

//******************************************************************//
// scheme for phiH
//******************************************************************//
OStringStream buf;

const tokenList& tokens = meshSurface.divScheme("div(phiH,YiSurf)");

label len = tokens.size();

forAll (tokens, tok)
{
    buf << meshSurface.divScheme("div(phiH,YiSurf)")[tok];
    if (--len)
    {
        buf << ' ';
    }
}

string phiHScheme = buf.str();

//define the reactive transport part
// Info << "Create chemical species mixture.\n" << endl;
// twoPhaseMultiComponentTransportMixture<inertMultiComponentMixture, inertMultiComponentMixture> speciesMixture(meshSurface, alpha1);

Info << "Reading solution component dictionary for surface flow." << endl;
const dictionary& solutionComponentDictSurf = transportProperties_surf.subDict("solutionSpecies");

speciesTable solutionSpeciesSurf
(
    solutionComponentDictSurf.toc()
);

Info << "Create the fluid solution concentration field.\n" << endl;
PtrList<volScalarField> YSurf(solutionSpeciesSurf.size());

forAll(solutionSpeciesSurf, i)
{
    Info << "Read the chemicals in surface water: " << solutionSpeciesSurf[i] << endl;
    IOobject header
    (
        solutionSpeciesSurf[i],
        meshSurface.time().timeName(),
        meshSurface,
        IOobject::NO_READ
    );

    if (header.typeHeaderOk<volScalarField>(true))
    {
        Info << "Read the chemical concentration from file: " << solutionSpeciesSurf[i] << endl;
        YSurf.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    solutionSpeciesSurf[i],
                    meshSurface.time().timeName(),
                    meshSurface,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                meshSurface
            )
        );
    }
    else
    {
        volScalarField YSurfdefault
        (
            IOobject
            (
                "YSurfdefault",
                meshSurface.time().timeName(),
                meshSurface,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            meshSurface
        );

        YSurf.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    solutionSpeciesSurf[i],
                    meshSurface.time().timeName(),
                    meshSurface,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                YSurfdefault
            )
        );
    }

}

//For surface species
wordList surfaceSpeciesNamesSurf;
if (transportProperties_surf.found("surfaceSpecies"))
{
    Info << "Reading surface species dictionary" << endl;
    const dictionary& surfaceSpeciesDictSurf = transportProperties_surf.subDict("surfaceSpecies");
    Info << "Reading surface species name\n" << endl;
    surfaceSpeciesNamesSurf = surfaceSpeciesDictSurf.toc();
}

speciesTable surfaceSpeciesSurf
(
    surfaceSpeciesNamesSurf
);

Info << "Create surface composition vector\n" << endl;
// PtrList<volScalarField> sYSurf(surfaceSpeciesSurf.size());

// forAll(surfaceSpeciesSurf, i)
// {
//     IOobject header
//     (
//         surfaceSpeciesSurf[i],
//         meshSurface.time().timeName(),
//         meshSurface,
//         IOobject::NO_READ
//     );

//     // check if field exists and can be read
//     if (header.typeHeaderOk<volScalarField>(true))
//     {
//         sYSurf.set
//         (
//             i,
//             new volScalarField
//             (
//                 IOobject
//                 (
//                     surfaceSpeciesSurf[i],
//                     meshSurface.time().timeName(),
//                     meshSurface,
//                     IOobject::MUST_READ,
//                     IOobject::AUTO_WRITE
//                 ),
//                 meshSurface
//             )
//         );
//     }
//     else
//     {
//         volScalarField Ydefault
//         (
//             IOobject
//             (
//                 "Ydefault",
//                 meshSurface.time().timeName(),
//                 meshSurface,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             ),
//             meshSurface
//         );

//         sYSurf.set
//         (
//             i,
//             new volScalarField
//             (
//                 IOobject
//                 (
//                     surfaceSpeciesSurf[i],
//                     meshSurface.time().timeName(),
//                     meshSurface,
//                     IOobject::NO_READ,
//                     IOobject::AUTO_WRITE
//                 ),
//                 Ydefault
//             )
//         );
//     }
// }

PtrList<volScalarField> sYSurf(surfaceSpeciesSurf.size());

forAll(surfaceSpeciesSurf, i)
{
    IOobject header
    (
        surfaceSpeciesSurf[i],
        meshSurface.time().timeName(),
        meshSurface,
        IOobject::NO_READ
    );

    // check if field exists and can be read
    if (header.typeHeaderOk<volScalarField>(true))
    {
        sYSurf.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    surfaceSpeciesSurf[i],
                    meshSurface.time().timeName(),
                    meshSurface,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                meshSurface
            )
        );
    }
    else
    {
        volScalarField Ydefault
        (
            IOobject
            (
                "Ydefault",
                meshSurface.time().timeName(),
                meshSurface,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            meshSurface
        );

        sYSurf.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    surfaceSpeciesSurf[i],
                    meshSurface.time().timeName(),
                    meshSurface,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                Ydefault
            )
        );
    }
}

Info << "Reading surface masters dictionary" << endl;
wordList surfaceMastersNamesSurf;
if (transportProperties_surf.found("surfaceMasters"))
{
    const dictionary& surfacesDictSurf = transportProperties_surf.subDict("surfaceMasters");
    Info << "Reading surface masters \n" << endl;
    surfaceMastersNamesSurf = surfacesDictSurf.toc();
}


speciesTable surfaceMastersSurf
(
    surfaceMastersNamesSurf
);

Info<< "Reading field Surf\n" << endl;
volScalarField SurfSurf
(
    IOobject
    (
        "SurfSurf",
        runTime.timeName(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface
);

Info << "Reading kinetic phases dictionary" << endl;
wordList kineticPhasesNamesSurf;
if (transportProperties_surf.found("kineticPhases"))
{
    Info << "Reading kinetic phases dictionary" << endl;
    const dictionary& kineticPhasesDictSurf = transportProperties_surf.subDict("kineticPhases");
    Info << "Reading kinetic phases names\n" << endl;
    kineticPhasesNamesSurf = kineticPhasesDictSurf.toc();
}

speciesTable kineticPhasesSurf
(
    kineticPhasesNamesSurf
);

// PtrList<volScalarField> RSurf(kineticPhasesSurf.size());
// forAll(kineticPhasesSurf, i)
// {
//     RSurf.set
//     (
//         i,
//         new volScalarField
//         (
//             IOobject
//             (
//                 "RSurf_" + kineticPhasesSurf[i],
//                 meshSurface.time().timeName(),
//                 meshSurface,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             meshSurface,
//             dimensionedScalar("RSurf_"+kineticPhasesSurf[i], dimensionSet(0,-2,-1,0,1,0,0), 0.0)
//         )
//     );
// }

Info << "The kinetic phase surf is " << kineticPhasesSurf.size() << endl;
PtrList<volScalarField> RSurf(kineticPhasesSurf.size());
forAll(kineticPhasesSurf, i)
{
    RSurf.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "RSurf_" + kineticPhasesSurf[i],
                meshSurface.time().timeName(),
                meshSurface,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            meshSurface,
            dimensionedScalar("RSurf_"+kineticPhasesSurf[i], dimensionSet(0,-2,-1,0,1,0,0), 0.0)
        )
    );
}

wordList selectedOutputNamesSurf;
if (transportProperties_surf.found("selectedOutput"))
{
    Info << "Reading selected output" << endl;
    const dictionary& selectedOutputDictSurf = transportProperties_surf.subDict("selectedOutput");
    Info << "Reading solution species name\n" << endl;
    selectedOutputNamesSurf = selectedOutputDictSurf.toc();
}


// Info << "Create selected output vector\n" << endl;
// PtrList<volScalarField> sOutSurf(selectedOutputNamesSurf.size());

// forAll(selectedOutputNamesSurf, i)
// {
//     sOutSurf.set
//     (
//         i,
//         new volScalarField
//         (
//             IOobject
//             (
//                 selectedOutputNamesSurf[i],
//                 meshSurface.time().timeName(),
//                 meshSurface,
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             meshSurface
//         )
//     );
// }

Info << "Create selected output vector\n" << endl;
PtrList<volScalarField> sOutSurf(selectedOutputNamesSurf.size());
Info << "The selected output from surface is " << selectedOutputNamesSurf.size() << endl;

forAll(selectedOutputNamesSurf, i)
{
    sOutSurf.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                selectedOutputNamesSurf[i],
                meshSurface.time().timeName(),
                meshSurface,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            meshSurface
        )
    );
}

//surface potential
volScalarField psiSurf
(
    IOobject
    (
        "psiSurf",
        runTime.timeName(),
        meshSurface,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface,
    dimensionedScalar("psiSurf", dimensionSet(1,2,-3,0,0,-1,0), 0.0)
);

//ionic strength
volScalarField ISurf
(
    IOobject
    (
        "ISurf",
        runTime.timeName(),
        meshSurface,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface,
    dimensionedScalar("ISub", dimMoles/dimVolume, 0.0)
);

volScalarField alpha1Surf
(
    IOobject
    (
        "alpha1Surf",
        runTime.timeName(),
        meshSubsurface
    ),
    meshSurface,
    scalar(1.0)
);


#include "createSurfaceFlowMRF.H"
