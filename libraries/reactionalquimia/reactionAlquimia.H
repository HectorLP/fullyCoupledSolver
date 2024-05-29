/*-----------------------------------------------------------------------------*\
Class
    Foam::reactionAlquimia module
 \*-----------------------------------------------------------------------------*/

#ifndef reactionAlquimia_H
#define reactionAlquimia_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"


#include "fvCFD.H"
#include "fvPatch.H"
#include "fvOptions.H"
#include "simpleControl.H"

//#include "pestc.h"

#include "alquimia_containers.h"
#include "alquimia_interface.h"
#include "alquimia_memory.h"
#include "alquimia_util.h"

namespace Foam
{
    class reactionAlquimia
    {
        protected:
            void initialize();
            void initializeInitialCondition();
            void initializeBoundaryCondition();
            int ncells;

            //double* pH;
            // double* mineralReactionRate;
            // double* saturation;
             double mineral_volume_fraction;
        public:
            AlquimiaProperties* chem_properties;
            AlquimiaState* chem_state;
            AlquimiaInterface chem;
            void* chem_engine_alquimia;
            AlquimiaEngineStatus chem_status;
            AlquimiaAuxiliaryData* chem_aux_data;
            AlquimiaAuxiliaryOutputData* chem_aux_output;
            AlquimiaSizes chem_sizes;
            AlquimiaProblemMetaData chem_metadata;
            AlquimiaGeochemicalCondition* chem_ic_geochem;
            AlquimiaGeochemicalCondition* chem_bc_geochem;
            AlquimiaState* chem_bc_state;
            AlquimiaAuxiliaryData* chem_aux_bc_data;
            AlquimiaAuxiliaryOutputData* chem_aux_bc_output;
            AlquimiaEngineFunctionality chem_engine_functionality;

            double* originalMineralVolumeFraction;
            double** accumulatedMineralVolumeChange;
            const char* chem_engine_ca;
            const char* chem_input_ca;

            const fvMesh& mesh_of;

            scalar temperature_of;
            scalar criticalValue_of;

            double waterDensity_of;

            word chem_engine_of;
            word crunchflow_input_file_of;
            bool initializeFromOpenFOAM;
            wordList& initialConditionList_of;
            wordList& boundaryConditionList_of;
            wordList& pHValue_of;
            wordList& mineralPhases_of;

            speciesTable primarySpecies_of;
            PtrList<volScalarField>& Y_primary_of;
            PtrList<volScalarField>& pH_of;
            volScalarField& indicator_of;
            volScalarField& bcIndicator_of;

            reactionAlquimia(const fvMesh&,
                            scalar,
                            scalar,
                            word,
                            word,
                            bool, 
                            wordList&,
                            wordList&,
                            wordList&,
                            const speciesTable&,
                            wordList&,
                            PtrList<volScalarField>&,
                            PtrList<volScalarField>&,
                            volScalarField&,
                            volScalarField&,
                            const dimensionedScalar&
            );
            virtual ~reactionAlquimia();
            void initializeMineralsBoundary();
            void updateReactiveBulkSurfaceArea();
            void keepMineralVolumeFraction();
            void overwriteChemicalYMultiphase(volScalarField&);
            void overwriteChemicalY();
            void keepPrimaryConcPositive();
            void reactionTimeStep(dimensionedScalar);
            void reactionTimeStepSaturation(dimensionedScalar, volScalarField&);
            void reactionTimeStepMovementofFace(PtrList<dimensionedScalar>&, volScalarField&);
    };
}

#endif
