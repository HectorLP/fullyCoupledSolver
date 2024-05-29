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
Class
    Foam::porousReactionsModels
Description
    geochemicalModel module for geometric porous media properties including
      - porosity
      - absolute permeability models
      - surface area models
SourceFiles
    porousReactionsModels.C
Bases on C. Soulaine's work

\*---------------------------------------------------------------------------*/

#ifndef porousReactionsModels_H
#define porousReactionsModels_H

#include "volFields.H"
#include "surfaceInterpolate.H"

#include "basicPorousReactionsModels.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class porousReactionsModels definition
\*---------------------------------------------------------------------------*/

class porousReactionsModels
{
    protected:
        const fvMesh &mesh_;

        dictionary porousReactionsModelsDict_;

        autoPtr<basicPorousReactionsModels> basicPorousReactionsModelsPtr_;

        porousReactionsModels(const porousReactionsModels&);

        PtrList<volScalarField>& Yf_;
    
    public:
        //constructors
        //-construct from components
        // porousReactionsModels
        // (
        //     const fvMesh& mesh,
        //     const dictionary& dict 
        // );

        // porousReactionsModels
        // (
        //     const fvMesh& mesh,
        //     const volScalarField& Ys,
        //     const dictionary& dict
        // );

        porousReactionsModels
        (
            const fvMesh& mesh,
            const dictionary& dict,
            PtrList<volScalarField>& Yf
        );

        ~porousReactionsModels()
        {}

        //-member functions
        const volScalarField& absolutePermeability() const 
        {
            return basicPorousReactionsModelsPtr_->absolutePermeability();
        }

        tmp<volScalarField> inversePermeability() const
        {
            return basicPorousReactionsModelsPtr_->inversePermeability();
        }
        const surfaceScalarField & Kf() const
        {
            return basicPorousReactionsModelsPtr_->Kf();
        }
        // const surfaceTensorField & KT() const 
        // {
        //     return basicPorousReactionsModelsPtr_->KT();
        // }
        const volScalarField & porosity() const 
        {
            return basicPorousReactionsModelsPtr_->porosity();
        }
        const volScalarField & dMinvdRho() const 
        {
            return basicPorousReactionsModelsPtr_->dMinvdRho();
        }

        const volTensorField & effectiveDispersionTensor() const 
        {
            return basicPorousReactionsModelsPtr_->effectiveDispersionTensor();
        }

        void update()
        {
            basicPorousReactionsModelsPtr_->updateFluidComposition();
            basicPorousReactionsModelsPtr_->updateMineralDistribution();
            basicPorousReactionsModelsPtr_->updatePorosity();
            basicPorousReactionsModelsPtr_->updateSurfaceArea();
            basicPorousReactionsModelsPtr_->updatePermeability();
            basicPorousReactionsModelsPtr_->updateDispersionTensor();
            basicPorousReactionsModelsPtr_->updateMinvdRho();
            Info << "Geochemical components are updated." << endl;
        }

        void updateFluidCompositionOnlyReactions()
        {
            basicPorousReactionsModelsPtr_->updateFluidCompositionOnlyReactionsNew();
        }

        void updateMineralsParameters()
        {
            basicPorousReactionsModelsPtr_->updateMineralDistribution();
            // basicPorousReactionsModelsPtr_->updatePorosity();
            basicPorousReactionsModelsPtr_->updateSurfaceArea();
            // basicPorousReactionsModelsPtr_->updatePermeability();
            // basicPorousReactionsModelsPtr_->updateDispersionTensor();
            basicPorousReactionsModelsPtr_->updateMinvdRho();
        }

        // PtrList<volScalarField>& Y()
        // {
        //     basicPorousReactionsModelsPtr_->Y();
        // }        

        // const PtrList<volScalarField>& Y_cons()
        // {
        //     basicPorousReactionsModelsPtr_->Y();
        // }

        // volScalarField& Y(const word& specieName)
        // {
        //     basicPorousReactionsModelsPtr_->Y(specieName);
        // }

        // const volScalarField& Y_icons(const word& specieName) const
        // {
        //     basicPorousReactionsModelsPtr_->Y(specieName);
        // }

        // const volScalarField& Y_icons(const label i) const
        // {
        //     basicPorousReactionsModelsPtr_->Y(i);
        // }

        // volScalarField& Y(const label& i)
        // {
        //     basicPorousReactionsModelsPtr_->Y(i);
        // }
};
}

#endif