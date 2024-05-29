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
    Foam::twoPhaseMixture
Description
    phreeqcRM module
SourceFiles
    phreeqc.C
\*---------------------------------------------------------------------------*/

#ifndef phreeqcRMModels_H
#define phreeqcRMModels_H

#include "volFields.H"
#include "speciesTable.H"
#include "runTimeSelectionTables.H"
#include "basicPorousReactionsModels.H"
#include "PhreeqcRM.h"

namespace Foam
{
namespace porousReactionsModels
{
class phreeqcRMModels : public basicPorousReactionsModels
{

    dictionary phreeqcDict_;
    protected:
    //protected member data
        dictionary transportPropertiesDict_;

        Switch activateUpdatePorosity_;

        Switch setComponentH2O_;

        scalar numThreads_;
        label nthread_;
        int nxyz_;

        scalar strangSteps_;

        word phreeqcInputFile_;
        word phreeqcDataBase_;

        PtrList<dictionary> mineralSubDict_;

        PtrList<Switch> activatePhaseEquilibrium_;

        PtrList<dimensionedScalar> Vm_;

        volScalarField linkDomainZonesToPhreeqcSol_;

        Switch cvODE;
        double cvODETol_;
        Switch use_SNIA;

        word UName_;

        //number of compnonent in the water
        int nComponent; 

        const volVectorField& U_;
        volScalarField pH_;

        PhreeqcRM phreeqc_;

        IRM_RESULT status;

        //Private member functions
    private:
        std::string generateEqPhasesInputString();

        std::string generateKineticsInputString();

        void updateKineticsParameters();

        void initializeMineralDistribution();

        void initializeFluidComposition();

        void activateSelectedOutput();

        void updatepH();

        void updatePorosityPhreeqc();

        phreeqcRMModels(const phreeqcRMModels&);
    
    public:

        TypeName("phreeqcRMModels");

        //Declare run-time constructor selection table 

        phreeqcRMModels
        (const fvMesh& mesh,
         const dictionary& dict
        );

        virtual ~phreeqcRMModels()
        {}

        //Member functions

        virtual void updateFluidComposition();

        virtual void updateMineralDistribution();

        virtual void updatePorosity();

        virtual void updateFluidCompositionOnlyReactions();

        virtual void updateFluidCompositionFromFlow(PtrList<volScalarField>&);

};

}

}

#endif