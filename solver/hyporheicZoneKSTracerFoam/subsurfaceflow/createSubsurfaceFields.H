//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////
Info<< "\nReading g for subsurface flow" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

word CFL =  runTime.controlDict().lookupOrDefault<word>("CFL", "Coats");
if (adjustTimeStep) adjustTimeStep=true; // to remove warnings at compilation
int CFLint = -1;
if (CFL == "Coats") CFLint = 0;
else if (CFL == "Courant") CFLint = 1; 
else
{
    FatalErrorIn
        (
            "in createFields.H"
        )
        << "CFL condition unrecongnized : Coats and Courant available" 
            << exit(FatalError);
}

Info<< "Reading field p\n" << endl;
volScalarField pSub
(
    IOobject
    (
        "pSub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface
);

Info<< "Reading field Sb \n" << endl;
volScalarField Sb
(
    IOobject
    (
        "Sb",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface
);

scalar dSmax(runTime.controlDict().lookupOrDefault<scalar>("dSmax",0.));

//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Sbmin(transportProperties.lookupOrDefault("Sbmin",dimensionedScalar("Sbmin",dimless,0)));
dimensionedScalar Sbmax(transportProperties.lookupOrDefault("Sbmax",dimensionedScalar("Sbmax",dimless,1.0)));
dimensionedScalar D_sub(transportProperties.lookup("D_sub"));

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////

autoPtr<incompressiblePhase> phasea = incompressiblePhase::New(meshSubsurface,transportProperties,"a");
volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& mua = phasea->mu();

autoPtr<incompressiblePhase> phaseb = incompressiblePhase::New(meshSubsurface,transportProperties,"b");
volVectorField& Ub = phaseb->U();
surfaceScalarField& phib = phaseb->phi();
const dimensionedScalar& rhob = phaseb->rho();
const dimensionedScalar& mub = phaseb->mu();	

dimensionedScalar Mmu(mub/mua);//viscosity ratio

// Relative permeability model 
autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New("krModel",transportProperties,Sb);

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info<< "Reading Porous Media Properties \n" << endl;

// Porosity	
dimensionedScalar epsScalar(transportProperties.lookupOrDefault("eps",dimensionedScalar("",dimless,0.)));
   
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshSubsurface
);

// Intrinsic permeability       
Info<< "Reading field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface
);

// interpolated permeability
surfaceScalarField Kf(fvc::interpolate(K,"K"));

Info<< "Reading porous transport properties" << endl;
dictionary porousTransport(transportProperties.subDict("porousTransport"));
word phaseName(porousTransport.lookupOrDefault<word>("phaseName","notDefined"));
dimensionedScalar Kd(porousTransport.lookupOrDefault("Kd",dimensionedScalar("Kd_default",dimensionSet(-1,3,0,0,0),0.)));
dimensionedScalar rs(porousTransport.lookupOrDefault("rs",dimensionedScalar("rs_default",dimensionSet(1,-3,0,0,0),0.)));
dimensionedScalar epsTotal(porousTransport.lookupOrDefault("epsTotal",dimensionedScalar("epsTotal_default",dimless,0.)));
dimensionedScalar lambda(porousTransport.lookupOrDefault("lambda",dimensionedScalar("lambda_default",dimensionSet(0,0,-1,0,0),0.)));

//Whether to solve the scalar C (non-reactive tracer)
Switch CSubTracerSwitch
(
    transportProperties.lookup("CSubTracerSwitch")
);

//Whether to solve reactive transport process in the subsurface
Switch ChemSubSwitch
(
    transportProperties.lookup("CheSubSwitch")
);

//Whether to include the pressure gradient of surface water to the flow calculation
Switch gradSurfPressureGradientSwitch
(
    transportProperties.lookup("gradSurfPressureGradientSwitch")
);

Info << nl << "porousTransport parameters" << nl << "{" << endl;
Info << "    " << "phaseName " << " : " << phaseName << endl;
Info << "    " << Kd.name() << " : " << Kd.value() << endl;
Info << "    " << rs.name() << " : " << rs.value() << endl;
Info << "    " << epsTotal.name() << " : " << epsTotal.value() << endl;
Info << "    " << lambda.name() << " : " << lambda.value() << endl;
Info << "} \n" << endl;

//- creation of dispersion model
autoPtr<dispersionModels> DeffModel= dispersionModels::New("DeffModel",transportProperties,meshSubsurface);

Info << "Computing retard coefficient R" << endl;
volScalarField R ("R", 1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps);
// volScalarField R
// (
//     IOobject
//     (

//     ),
//     1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps.value();
// );
/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info<< "Reading field USub\n" << endl;
volVectorField USub
(
    IOobject
    (
        "USub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Ua + Ub
);

// #include "createPhi.H"
// surfaceScalarField phiP = phi;

Info<< "Reading/calculating face flux field phiSub\n" << endl;

surfaceScalarField phiSub
(
    IOobject
    (
        "phiSub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(USub)
);
surfaceScalarField phiP = phiSub;
/////////////////////////////////////////////////////////////////////////////
////////////////////////// CAPILLARY PRESSURE ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////

scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",0.));
autoPtr<capillarityModel> pcModel = capillarityModel::New("pcModel",transportProperties,Sb);

volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        meshSubsurface,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    meshSubsurface,
    dimensionedScalar("",dimless/dimTime,0)
);

//- Looking for 2nd order time-scheme for source term
bool isBackwardScheme = false;
if (word(meshSubsurface.ddtScheme("source")) == "backward")
{
    isBackwardScheme = true;
}

Info<< "Reading field C\n" << endl;
volScalarField C_Sub
(
    IOobject
    (
        "C_Sub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface
);

// Info << "Reading solution component dictionary" << endl;
// const dictionary& solutionComponentDict = transportProperties.subDict("solutionSpecies");
// Info << "Reading solution component name\n" << endl;
// speciesTable solutionSpecies
// (
//     solutionComponentDict.toc()
// );
// Info << "Reading each solution component's diffusion coefficient." << endl;
// PtrList<dimensionedScalar> DY(solutionSpecies.size());

// forAll(solutionComponentDict, i)
// {
//     const dictionary& subdict = solutionComponentDict.subDict(solutionSpecies[i]);
//     DY.set
//     (
//         i, 
//         new dimensionedScalar(subdict.lookup("D"))
//     );
// }


// Info << "Create the fluid solution concentration field.\n" << endl;
// PtrList<volScalarField> Y_Sub(solutionSpecies.size());

// forAll(solutionSpecies, i)
// {
//     IOobject header
//     (
//         solutionSpecies[i],
//         meshSubsurface.time().timeName(),
//         meshSubsurface,
//         IOobject::NO_READ
//     );

//     if (header.typeHeaderOk<volScalarField>(true))
//     {
//         Y_Sub.set
//         (
//             i,
//             new volScalarField
//             (
//                 IOobject
//                 (
//                     solutionSpecies[i],
//                     meshSubsurface.time().timeName(),
//                     meshSubsurface,
//                     IOobject::MUST_READ,
//                     IOobject::AUTO_WRITE
//                 ),
//                 meshSubsurface
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
//                 meshSubsurface.time().timeName(),
//                 meshSubsurface,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             ),
//             meshSubsurface
//         );

//         Y_Sub.set
//         (
//             i,
//             new volScalarField
//             (
//                 IOobject
//                 (
//                     solutionSpecies[i],
//                     meshSubsurface.time().timeName(),
//                     meshSubsurface,
//                     IOobject::NO_READ,
//                     IOobject::AUTO_WRITE
//                 ),
//                 Ydefault
//             )
//         );
//     }
// }