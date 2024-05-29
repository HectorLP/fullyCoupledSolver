/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "basicPorousReactionsModels.H"
#include "fvcDdt.H"

//* * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(basicPorousReactionsModels, 0);
    defineRunTimeSelectionTable(basicPorousReactionsModels, dictionary);
} 

//* * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //
Foam::basicPorousReactionsModels::basicPorousReactionsModels
(
    const fvMesh& mesh,
    const dictionary& dict 
)
:
    mesh_(mesh),
    mineralList_(dict.lookup("mineral")),
    componentList_(dict.lookup("solutionComponent")),
    Ys_(mineralList_.size()),
    inertMineral_
    (
        IOobject
        (
            "inertMineral",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("inertMineral", dimless, 0.0),
        "zeroGradient"
    ),
    eps_ 
    (
        IOobject
        (
            "eps",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("eps",dimless,1.0),
        "zeroGradient"
    ),
    eps0_
    (
        IOobject
        (
            "eps0", 
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("eps0", dimless, 1.0),
        "zeroGradient"
    ),
    rhol_
    (
        dict.lookup("rhol")
    ),
    rhos_(mineralList_.size()),
    dMinvdRho_
    (
        IOobject
        (
            "dMinvdRho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("dMinvdRho",dimless/dimTime, 0.0),
        "zeroGradient"
    ),
    porousMedia_(mineralList_.size()),
    absolutePermeabilityModelPtr_
    (
        absolutePermeabilityModel::New(mesh, dict)
    ),

    dispersionTensorModelPtr_
    (
        dispersionTensorModel::New(mesh, dict)
    ),
    phiName_(dict.lookupOrDefault<word>("phi","phi")),
    phi_(mesh.lookupObject<surfaceScalarField>(phiName_))
{
    forAll(mineralList_, s)
    {
        word currentMineral = mineralList_[s];
        Info << "The current mineral is " << currentMineral << endl;

        Ys_.set
        (
            s,
            new volScalarField
            (
                IOobject
                (
                    "Ys."+mineralList_[s],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
        Ys_[s].write();

        rhos_.set
        (
            s,
            new dimensionedScalar
            (
                dict.subDict(currentMineral+"Properties").lookup("rhos")
            )
        );

        porousMedia_.set
        (
            s,
            new porousModel
            (
                mesh,
                mineralList_[s],
                Ys_[s],
                dict 
            )
        );
    }
    updatePorosity();
}

// -------------------------------------------------------------------------//
void Foam::basicPorousReactionsModels::updatePorosity()
{
    eps_=0.0*eps_;
    forAll(mineralList_,s)
    {
        eps_+=Ys_[s];
    }
    eps_=1.-eps_-inertMineral_;
}

void Foam::basicPorousReactionsModels::updateMinvdRho()
{
    dMinvdRho_=0.0*dMinvdRho_;
    forAll(mineralList_,s)
    {
        dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s]) * (1./rhol_-1./rhos_[s]);
    }
}
