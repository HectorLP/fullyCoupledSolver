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
    Foam::basicPorousReactionsModels
Description
    basicPorousReactionsModels module
SourceFiles
    basicPorousReactionsModels.C
Authors
    17/02/2020 : C. Soulaine
\*---------------------------------------------------------------------------*/

#ifndef basicPorousReactionsModels_H
#define basicPorousReactionsModels_H

#include "volFields.H"
#include "speciesTable.H"
#include "runTimeSelectionTables.H"
#include "porousSinglePhaseTransportModels.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class basicPorousReactionsModels Declaration
\*---------------------------------------------------------------------------*/
class basicPorousReactionsModels
{
    protected:
        const fvMesh &mesh_;

        speciesTable mineralList_;
        
        speciesTable componentList_;

        PtrList<volScalarField> Ys_;

        volScalarField inertMineral_;

        volScalarField eps_;
        volScalarField eps0_;

        dimensionedScalar rhol_;
        
        PtrList<dimensionedScalar> rhos_;

        volScalarField dMinvdRho_;

        PtrList<porousModel> porousMedia_;

        autoPtr<absolutePermeabilityModel> absolutePermeabilityModelPtr_;

        autoPtr<dispersionTensorModel> dispersionTensorModelPtr_;

        word phiName_;

        const surfaceScalarField &phi_;

        // chemical concentration fields
        PtrList<volScalarField> Y_;

        basicPorousReactionsModels(const basicPorousReactionsModels&);
    
    public:
        TypeName("basicPorousReactionsModels");

        declareRunTimeSelectionTable
        (
            autoPtr,
            basicPorousReactionsModels,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict 
            ),
            (mesh, dict)
        );

        //
        static autoPtr<basicPorousReactionsModels> New 
        (
            const fvMesh& mesh,
            const dictionary& dict 
        );

        //construct from components

        basicPorousReactionsModels
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        virtual ~basicPorousReactionsModels()
        {}

        //member functions
        virtual void updatePorosity();

        tmp<volScalarField> porosity() const
        {
            return eps_;
        }

        tmp<volScalarField> absolutePermeability() const 
        {
            return absolutePermeabilityModelPtr_->absolutePermeability();
        }

        tmp<volScalarField> inversePermeability() const 
        {
            return absolutePermeabilityModelPtr_->inversePermeability();
        }

        tmp<surfaceScalarField> Kf() const 
        {
            return absolutePermeabilityModelPtr_->Kf();
        }

        void updatePermeability()
        {
            absolutePermeabilityModelPtr_->updatePermeability();
        }

        void updateSurfaceArea()
        {
            forAll(porousMedia_, s)
            {
                porousMedia_[s].updateSurfaceArea();
            }
        }

        tmp<volTensorField> effectiveDispersionTensor() const 
        {
            return dispersionTensorModelPtr_->effectiveDispersionTensor();
        }

        void updateDispersionTensor()
        {
            dispersionTensorModelPtr_->updateDispersionTensor();
        }

        tmp<volScalarField> dMinvdRho() const
        {
            return dMinvdRho_;
        }

        virtual void updateMinvdRho();

        virtual void updateFluidComposition() = 0;

        virtual void updateMineralDistribution() = 0;

        virtual void updateFluidCompositionOnlyReactions() = 0;

        inline PtrList<volScalarField>& Y();

        inline const PtrList<volScalarField>& Y() const;

        inline volScalarField& Y(const label i);

        inline const volScalarField& Y(const label i) const;

        inline volScalarField& Y(const word& specieName);

        inline const volScalarField& Y(const word& specieName) const;

};
}

#include "basicPorousReactionsModelsI.H"

#endif 
