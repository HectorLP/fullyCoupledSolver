/*-----------------------------------------------------------------------------*\
Definitions of reactionAlquimia class with its functions
\*-----------------------------------------------------------------------------*/
#include "petsc.h"
#include "reactionAlquimia.H"
#include "fvcAverage.H"

/*-----------------------------------------------------------------------------*\
Function to initialize the boundary condition and initial condition from CrunchFLow
through Alquimia
/*-----------------------------------------------------------------------------*/
void Foam::reactionAlquimia::initialize()
{
    char help[] = "Integrate CrunchFlow to ADE solver of OpenFOAM.";
    char *a  =  " ";//" ";
    char **ar = &a;
    char ***argv = &ar;
    int numLine = 1;
    int *argc = &numLine;
    PetscErrorCode ierr;
    ierr = PetscInitialize(argc, argv, (char*)0, help);
    ierr = PetscInitializeFortran();
    int ncells = mesh_of.cells().size();
    originalMineralVolumeFraction = new double [ncells];
    //get number of cells
    //patchBoundary = patchesNameList;
    // Info << "The number of initial condition names in CrunchFlow (including NULL) " << initialConditionList_of.size() << endl;
    // Info << "The number of boundary names in CrunchFlow (including NULL) "  << boundaryConditionList_of.size() << endl;
    /*Convert string to const char*/
    // Info << "Define parameters used in CrunchFlow." << endl;
    std::string chem_engine_crunch;
    chem_engine_crunch = chem_engine_of;
    std::cout << chem_engine_crunch << "\n";
    chem_engine_ca = chem_engine_crunch.c_str();
    // printf("The chemical engine is %s.\n", chem_engine_ca);
    std::string chem_input_crunch;
    chem_input_crunch = crunchflow_input_file_of;
    chem_input_ca = chem_input_crunch.c_str();
    //Start to initialize the elements from alquimia
    // Info << "Start to set up the chemistry engine." << endl;
    AllocateAlquimiaEngineStatus(&chem_status);
    // Info << "Start to set up the interface for chemistry engine." << endl;
    CreateAlquimiaInterface(chem_engine_ca, &chem, &chem_status);
    // Info << "Finish the interface setup." << endl;
    if (chem_status.error != 0)
    {
        Info << "Error happened to set up the chemistry engine." << endl;
        alquimia_error(chem_status.message);
        std::abort();
    }
    // Info << "Start to set up the functionality" << endl;
    bool hand_off_ca = true;
    AlquimiaEngineFunctionality chem_engine_functionality;
    chem.Setup(chem_input_ca, hand_off_ca, &chem_engine_alquimia, \
                &chem_sizes, &chem_engine_functionality, \
                &chem_status);
    if (chem_status.error != 0)
    {
        Info << "Error happened to set up storage for engine." << endl;
        Info << chem_status.message << endl;
        std::abort();
    }
    // Info << "Start the Metadata part." << endl;
    AllocateAlquimiaProblemMetaData(&chem_sizes, &chem_metadata);
    //Metadata
    // Info << "Get chem problem data." << endl;
    chem.GetProblemMetaData(&chem_engine_alquimia, &chem_metadata, \
                            &chem_status);
    if (chem_status.error != 0)
    {
        alquimia_error("Error happened for metadata: %s.", chem_status.message);
        std::abort();
    }
    //initialize the accumulated mineral volume change for each mineral
    accumulatedMineralVolumeChange = new double* [chem_sizes.num_minerals];
    for (int i = 0; i < chem_sizes.num_minerals; ++i)
    {
      accumulatedMineralVolumeChange[i] = new double [ncells];
    }
    for (int i = 0; i < chem_sizes.num_minerals; ++i)
    {
      for (long j = 0; j < ncells; ++j)
      {
        accumulatedMineralVolumeChange[i][j] = 0.0;
      }
    }
    //Memory for initial condition
    initializeInitialCondition();
    initializeBoundaryCondition();
}

/*Initialize the initial condition of geochemistry for the simulation domain*/
void Foam::reactionAlquimia::initializeInitialCondition()
{
  static const double aqueous_pressure = 101325.0;
  int ncells = mesh_of.cells().size();
  // Info << "Set memory for initial condition." << endl;
  chem_ic_geochem = new AlquimiaGeochemicalCondition[initialConditionList_of.size()];
  chem_properties = new AlquimiaProperties [ncells];
  chem_state = new AlquimiaState [ncells];
  chem_aux_data = new AlquimiaAuxiliaryData [ncells];
  chem_aux_output = new AlquimiaAuxiliaryOutputData [ncells];
  // Info << "The total number of cells is " << ncells << endl;
  forAll(mesh_of.cells(), celli)
  {
      AllocateAlquimiaState(&chem_sizes, &chem_state[celli]);
      AllocateAlquimiaProperties(&chem_sizes, &chem_properties[celli]);
      AllocateAlquimiaAuxiliaryData(&chem_sizes, &chem_aux_data[celli]);
      AllocateAlquimiaAuxiliaryOutputData(&chem_sizes, &chem_aux_output[celli]);
  }
  forAll(initialConditionList_of, i)
  {
      std::string ic_name = initialConditionList_of[i];  //icInCrunchFlow_of[i];
      const char* tmp_ic = ic_name.c_str();
      // Info << "The name for the initial boundary is " << ic_name << endl;
      if (tmp_ic != NULL)
      {
          AllocateAlquimiaGeochemicalCondition(strlen(tmp_ic), 0, \
                                      0, &chem_ic_geochem[i]);
          strcpy(chem_ic_geochem[i].name, tmp_ic);
      }
      else
      {
          chem_ic_geochem[i].name = NULL;
      }
      //delete  tmp_ic;
  }
  forAll(mesh_of.cells(), i)
  {
      chem_properties[i].volume = mesh_of.V()[i];
      chem_properties[i].saturation = 1.0;
      chem_state[i].water_density = waterDensity_of;
      chem_state[i].temperature = temperature_of;
      chem_state[i].aqueous_pressure = aqueous_pressure;
      int tmpIndicator = indicator_of[i];
      // Info << "The cell index is " << i << endl;
      // Info << "The indicator for the initial condition is " << tmpIndicator << endl;
      // Info << "The name of the initial condition is " << chem_ic_geochem[tmpIndicator].name << endl;
      chem.ProcessCondition(&chem_engine_alquimia, &chem_ic_geochem[tmpIndicator], \
          &chem_properties[i], &chem_state[i], &chem_aux_data[i], &chem_status);
      if (chem_status.error != 0)
      {
          Info << "Error happens to the process of initialization with area: " << initialConditionList_of[tmpIndicator] << endl;
          chem_status.message;
          std::abort();
      }
      chem.GetAuxiliaryOutput(&chem_engine_alquimia, &chem_properties[i], &chem_state[i], \
      &chem_aux_data[i], &chem_aux_output[i], &chem_status);
  }
  if (initializeFromOpenFOAM == true)
  {
      // Info << "Overwrite the initial condition" << endl;
      // std::cin.get();
      keepPrimaryConcPositive();
  }
}

/*Initialize the boundary condition of geochemistry except walls having minerals for dissolution*/
void Foam::reactionAlquimia::initializeBoundaryCondition()
{
  //Memory for boundary condition
  int ncells = mesh_of.cells().size();
  chem_bc_geochem = new AlquimiaGeochemicalCondition[boundaryConditionList_of.size()];
  chem_bc_state = new AlquimiaState [boundaryConditionList_of.size()];
  chem_aux_bc_data = new AlquimiaAuxiliaryData [boundaryConditionList_of.size()];
  chem_aux_bc_output = new AlquimiaAuxiliaryOutputData [boundaryConditionList_of.size()];
  // Info << "Set up boundary condition in geochemical part." << endl;
  forAll(boundaryConditionList_of, i)
  {
      std::string bc_name = boundaryConditionList_of[i];
      const char* tmp_bc = bc_name.c_str();
      // Info << "The name for the boundary is " << bc_name << endl;
      if (tmp_bc != NULL)
      {
          AllocateAlquimiaGeochemicalCondition(strlen(tmp_bc), 0, \
                                  0, &chem_bc_geochem[i]);
          strcpy(chem_bc_geochem[i].name, tmp_bc);
      }
      else
      {
          chem_bc_geochem[i].name = NULL;
      }
      AllocateAlquimiaState(&chem_sizes, &chem_bc_state[i]);
      AllocateAlquimiaAuxiliaryData(&chem_sizes, &chem_aux_bc_data[i]);
      AllocateAlquimiaAuxiliaryOutputData(&chem_sizes, &chem_aux_bc_output[i]);
  }
  // Info << "Initial reaction part of the boundary in CrunchFlow through Alquimia." << endl;
  static const double aqueous_pressure = 101325.0;
  forAll(boundaryConditionList_of, i)
  {
      chem_bc_state[i].water_density = waterDensity_of;
      chem_bc_state[i].temperature = temperature_of;
      //chem_bc_state[i].porosity = porosity_ca;
      chem_bc_state[i].aqueous_pressure = aqueous_pressure;
      if (chem_bc_geochem[i].name != NULL)
      {
          chem.ProcessCondition(&chem_engine_alquimia, &chem_bc_geochem[i], \
                          &chem_properties[ncells - 1], &chem_bc_state[i], \
                          &chem_aux_bc_data[i], &chem_status);
          if (chem_status.error != 0)
          {
              Info << "Has error with setting up a boundary condition named: " << boundaryConditionList_of[i] << endl;
              Info << chem_status.message << endl;
              chem_status.error;
          }
          chem.GetAuxiliaryOutput(&chem_engine_alquimia, &chem_properties[ncells - 1], \
                                &chem_bc_state[i], &chem_aux_bc_data[i], &chem_aux_bc_output[i], \
                                &chem_status);
      }
      else
      {
          for(int c = 0; c < chem_bc_state[i].total_mobile.size; ++c)
          {
              chem_bc_state[i].total_mobile.data[c] = chem_state[ncells - 1].total_mobile.data[c];
          }
      }
  }
  // Info << "Finish the initialization for the boundaries' chemical settings." << endl;
}

/*constructor of the class*/
Foam::reactionAlquimia::reactionAlquimia
(
    const fvMesh& mesh, scalar temperature, scalar criticalValue, \
    word chem_engine, \
    word crunchflow_input_file, \
    bool initializeFromOF, \
    wordList& initialConditionList, wordList& boundaryConditionList, wordList & pHValue, \
    const speciesTable& primarySpecies, wordList& mineralPhases,
    PtrList<volScalarField>& Y_primary, \
    PtrList<volScalarField>& pH, volScalarField& indicator, volScalarField& bcIndicator, \
    const dimensionedScalar& waterDensity
):
mesh_of(mesh),
initialConditionList_of(initialConditionList),
boundaryConditionList_of(boundaryConditionList),
pHValue_of(pHValue),
primarySpecies_of(primarySpecies),
mineralPhases_of(mineralPhases),
Y_primary_of(Y_primary),
pH_of(pH),
indicator_of(indicator),
bcIndicator_of(bcIndicator)
{
    temperature_of = temperature;
    criticalValue_of = criticalValue;
    chem_engine_of = chem_engine;
    crunchflow_input_file_of = crunchflow_input_file;
    initializeFromOpenFOAM = initializeFromOF;
    waterDensity_of = waterDensity.value();
    initialize();
}

/*Destroy the class for ending the simulation*/
Foam::reactionAlquimia::~reactionAlquimia()
{
    /*Delete parameters of Alquimia for using CrunchFlow*/
    int ncells = mesh_of.cells().size();
    // Info << "Delete parameters for Alquimia and shut it down." << endl;
    //Destroy chemistry data
    for (int i = 0; i < initialConditionList_of.size(); ++i)
    {
        if (&chem_ic_geochem[i].name != NULL)
        {
            FreeAlquimiaGeochemicalCondition(&chem_ic_geochem[i]);
        }
    }
    for (int i = 0; i < boundaryConditionList_of.size(); ++i)
    {
        if (&chem_bc_geochem[i].name != NULL)
        {
            FreeAlquimiaGeochemicalCondition(&chem_bc_geochem[i]);
        }
        FreeAlquimiaState(&chem_bc_state[i]);
        FreeAlquimiaAuxiliaryData(&chem_aux_bc_data[i]);
    }
    delete [] originalMineralVolumeFraction;
    delete [] accumulatedMineralVolumeChange;
    delete [] chem_bc_geochem;
    delete [] chem_bc_state;
    delete [] chem_aux_bc_data;
    delete [] chem_aux_bc_output;
    delete [] chem_ic_geochem;

    for (int i = 0; i < ncells; ++i)
    {
        FreeAlquimiaState(&chem_state[i]);
        FreeAlquimiaProperties(&chem_properties[i]);
        FreeAlquimiaAuxiliaryData(&chem_aux_data[i]);
        FreeAlquimiaAuxiliaryOutputData(&chem_aux_output[i]);
    }
    delete [] chem_state;
    delete [] chem_properties;
    delete [] chem_aux_data;
    delete [] chem_aux_output;
    FreeAlquimiaProblemMetaData(&chem_metadata);

    //Destroy chemistry engine
    chem.Shutdown(&chem_engine_alquimia, &chem_status);
    FreeAlquimiaEngineStatus(&chem_status);
    delete [] chem_engine_ca;
    delete [] chem_input_ca;
    PetscErrorCode ierr;
    ierr = PetscFinalize();
}

//initialize the cell which has mineral solid phase and the cell which belongs  to the pore space
//The mineral volume fraction and bulk surface area of cell having minerals are non-zero
void Foam::reactionAlquimia::initializeMineralsBoundary()
{
    //get all the types of the boundary
    // Info << "Initialize the bulk surface area for the reactive solid phase." << endl;
    // std::vector<long> cellHavingPatch;
    forAll(mesh_of.boundaryMesh(), patchI)
    {
        string patchName = mesh_of.boundaryMesh()[patchI].name();
        // Info << "The patch index is " << patchI << endl;
        // Info << "The patch name is " << patchName << endl;
        label patchIDs = mesh_of.boundaryMesh().findPatchID(patchName);
        const polyPatch& myPatches = mesh_of.boundaryMesh()[patchIDs];
        double surfaceArea, cellVolume;
        //In 2D, 'frontandback' patch covers the entire domain, so it is not included here to determine
        //whether the mineral volume fraction and bulk surface area of a cell should be changed or not
        if (myPatches.size() < mesh_of.cells().size() and myPatches.type()== "wall")
        {
          forAll(myPatches, faces)
          {
              int tmpBC_indicator = bcIndicator_of.boundaryFieldRef()[patchIDs][faces];
              int cellIndex = mesh_of.boundary()[patchIDs].faceCells()[faces];
              // Info << "The cell index is " << cellIndex << endl;
              // Info << "The boundary indicator is " << tmpBC_indicator << endl;
              int tmpCount = 0;
              if (tmpBC_indicator >= 0)
              {
                  for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
                  {
                      // Info << "The mineral volume is " << chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j] << endl;
                      if (chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j] > 0.0)
                      {
                          tmpCount += 1;
                      }
                  }
                  if (tmpCount > 0)
                  {

                      // cellHavingPatch.push_back(cellIndex);
                      // Info << "The cell has mineral" << endl;
                      surfaceArea = mesh_of.magSf().boundaryField()[patchIDs][faces];
                      cellVolume = mesh_of.V()[cellIndex];
                      chem.ProcessCondition(&chem_engine_alquimia, &chem_bc_geochem[tmpBC_indicator], \
                                      &chem_properties[cellIndex], &chem_state[cellIndex], \
                                      &chem_aux_data[cellIndex], &chem_status);
                      for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
                      {
                          chem_state[cellIndex].mineral_specific_surface_area.data[j] = surfaceArea / cellVolume;
                          //chem_state[cellIndex].mineral_volume_fraction.data[j] = chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j];
                      }
                      forAll(primarySpecies_of, j)
                      {
                          Y_primary_of[j].boundaryFieldRef()[patchIDs][faces] = \
                          chem_bc_state[tmpBC_indicator].total_mobile.data[j] * 1000.;
                          // Info << "The concentration value is " << Y_primary_of[j].boundaryFieldRef()[patchIDs][faces] << endl;
                      }
                  }
                  else
                  {
                      // Info << "The cell does not have any minerals." << endl;
                      for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
                      {
                          chem_state[cellIndex].mineral_specific_surface_area.data[j] = 0.0;
                          chem_state[cellIndex].mineral_volume_fraction.data[j] = 0.0;
                      }
                  }
              }
          }
        }
    }
    // Info << "Finish updating the bulk surface and mineral volume fraction values." << endl;
    //for the other areas in the domain, the volume fraction and bulk surface area should be 0
    // forAll(mesh_of.cells(), i)
    // {
    //     if (std::count(cellHavingPatch.begin(), cellHavingPatch.end(), i))
    //     {
    //         continue;
    //     }
    //     else
    //     {
    //         for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
    //         {
    //             chem_state[i].mineral_specific_surface_area.data[j] = 0.0;
    //             chem_state[i].mineral_volume_fraction.data[j] = 0.0;
    //         }
    //     }
    // }
}


//update the reactive bulk surface area of the cell having the solid phase patch
//by the calculation A = patch_area/cell_volume.
void Foam::reactionAlquimia::updateReactiveBulkSurfaceArea()
{
    //get all the types of the boundary
    // Info << "Update the bulk surface area for the reactive solid phase." << endl;
    // std::vector<long> cellHavingPatch;
    forAll(mesh_of.boundaryMesh(), patchI)
    {
        string patchName = mesh_of.boundaryMesh()[patchI].name();
        // Info << "The patch index is " << patchI << endl;
        // Info << "The patch name is " << patchName << endl;
        label patchIDs = mesh_of.boundaryMesh().findPatchID(patchName);
        const polyPatch& myPatches = mesh_of.boundaryMesh()[patchIDs];
        double surfaceArea, cellVolume;
        if (myPatches.size() < mesh_of.cells().size() and myPatches.type()== "wall")
        {
          forAll(myPatches, faces)
          {

              int tmpBC_indicator = bcIndicator_of.boundaryFieldRef()[patchIDs][faces];
              int cellIndex = mesh_of.boundary()[patchIDs].faceCells()[faces];
              // Info << "The cell index is " << cellIndex << endl;
              // Info << "The boundary indicator is " << tmpBC_indicator << endl;
              int tmpCount = 0;
              if (tmpBC_indicator >= 0)
              {
                  for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
                  {
                      if (chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j] > 0.0)
                      {
                          tmpCount += 1;
                      }
                  }
                  if (tmpCount > 0)
                  {
                      // cellHavingPatch.push_back(cellIndex);
                      surfaceArea = mesh_of.magSf().boundaryField()[patchIDs][faces];
                      cellVolume = mesh_of.V()[cellIndex];
                      for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
                      {
                          chem_state[cellIndex].mineral_specific_surface_area.data[j] = surfaceArea / cellVolume;
                      }
                  }
                  else
                  {
                      for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
                      {
                          chem_state[cellIndex].mineral_specific_surface_area.data[j] = 0.0;
                      }
                  }
              }
           }
        }
    }
    // Info << "Finish updating the bulk surface values." << endl;
    //for the other areas in the domain, the volume fraction and bulk surface area should be 0
    // forAll(mesh_of.cells(), i)
    // {
    //     if (std::count(cellHavingPatch.begin(), cellHavingPatch.end(), i))
    //     {
    //         continue;
    //     }
    //     else
    //     {
    //         for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
    //         {
    //             chem_state[i].mineral_specific_surface_area.data[j] = 0.0;
    //         }
    //     }
    // }
}

/*The mineral volume fraction should be always positive to avoid zero value for the reactions*/
void Foam::reactionAlquimia::keepMineralVolumeFraction()
{
    //get all the types of the boundary
    // Info << "Keep the positive mineral volume fraction for the reactive solid phase." << endl;
    // std::vector<long> cellHavingPatch;
    forAll(mesh_of.boundaryMesh(), patchI)
    {
        string patchName = mesh_of.boundaryMesh()[patchI].name();
        // Info << "The patch index is " << patchI << endl;
        // Info << "The patch name is " << patchName << endl;
        label patchIDs = mesh_of.boundaryMesh().findPatchID(patchName);
        const polyPatch& myPatches = mesh_of.boundaryMesh()[patchIDs];
        double surfaceArea, cellVolume;
        if (myPatches.size() < mesh_of.cells().size() and myPatches.type()== "wall")
        {
          forAll(myPatches, faces)
          {
              int tmpBC_indicator = bcIndicator_of.boundaryFieldRef()[patchIDs][faces];
              int cellIndex = mesh_of.boundary()[patchIDs].faceCells()[faces];
              // Info << "The cell index is " << cellIndex << endl;
              // Info << "The boundary indicator is " << tmpBC_indicator << endl;
              int tmpCount = 0;
              if (tmpBC_indicator >= 0)
              {
                  for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
                  {
                      // Info << "The mineral volume is " << chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j] << endl;
                      if (chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j] > 0.0)
                      {
                          tmpCount += 1;
                      }
                  }
                  if (tmpCount > 0)
                  {
                      // cellHavingPatch.push_back(cellIndex);
                      surfaceArea = mesh_of.magSf().boundaryField()[patchIDs][faces];
                      cellVolume = mesh_of.V()[cellIndex];
                      for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
                      {
                          //Calculated the volume change due to the mineral dissolution
                          double tmpOriginalVF = chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j];
                          double tmpVF = chem_state[cellIndex].mineral_volume_fraction.data[j];
                          accumulatedMineralVolumeChange[j][cellIndex] += (tmpOriginalVF - tmpVF) * cellVolume;
                          // chem_state[cellIndex].mineral_specific_surface_area.data[j] = surfaceArea / cellVolume;
                          chem_state[cellIndex].mineral_volume_fraction.data[j] = chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j];
                      }
                  }
                  else
                  {
                      for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
                      {
                          // chem_state[cellIndex].mineral_specific_surface_area.data[j] = 0.0;
                          chem_state[cellIndex].mineral_volume_fraction.data[j] = 0.0;
                      }
                  }
              }
          }
        }
    }
    // Info << "Finish updating mineral volume fraction values." << endl;
    //for the other areas in the domain, the volume fraction and bulk surface area should be 0
    // forAll(mesh_of.cells(), i)
    // {
    //     if (std::count(cellHavingPatch.begin(), cellHavingPatch.end(), i))
    //     {
    //         continue;
    //     }
    //     else
    //     {
    //         for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
    //         {
    //             // chem_state[i].mineral_specific_surface_area.data[j] = 0.0;
    //             chem_state[i].mineral_volume_fraction.data[j] = 0.0;
    //         }
    //     }
    // }
}

/*Function to overwrite the primary chemical species concentration on the OpenFoam
boundary and concentration values in the simulation domain.*/
void Foam::reactionAlquimia::overwriteChemicalY()
{
  // Info << "Start to overwrite primary species concentrations of each cell from CrunchFlow." << endl;
  forAll(mesh_of.cells(), i)
  {
      forAll(primarySpecies_of, j)
      {
	  if (chem_metadata.positivity.data[j] == 1 and chem_state[i].total_mobile.data[j]  <= 0)
	 {
	     Y_primary_of[j][i] = 1.0e-20;
	 }
	 else
	 {
        Y_primary_of[j][i] = chem_state[i].total_mobile.data[j] * 1000.0;
	 }
      }

  }
  //Overwrite the value on the boundary
  //Use mineral volume fraction to differentiate the reactive surface boundary from others (inlet/outlet)
  // Info << "Start to overwrite the boundary condition of primary species from CrunchFlow." << endl;
  forAll(mesh_of.boundaryMesh(), patchI)
  {
      string patchName = mesh_of.boundaryMesh()[patchI].name();
      // Info << "The boundary name in OpenFOAM " << patchName << endl;
      label patchIDs = mesh_of.boundaryMesh().findPatchID(patchName);
      const fvPatch& myPatches = mesh_of.boundary()[patchIDs];
      //Info << "The current patch size is " << myPatches.size() << endl;
      if (myPatches.size() < mesh_of.cells().size() and myPatches.type()== "patch")
      {
        forAll(myPatches, faces)
        {
            //Info << "The index of face is " << faces << endl;
            int tmpBC_indicator = bcIndicator_of.boundaryFieldRef()[patchIDs][faces];
            //Info << "The boundary indicator is " << tmpBC_indicator << endl;
            int tmpCountMineralV = 0;
            for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
            {
                // Info << "The volume fraction is " << chem_state[cellIndex].mineral_volume_fraction.data[j] << endl;
                if (chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j] > 0.0)
                {
                    tmpCountMineralV += 1;
                }
            }
            if (tmpCountMineralV == 0)
            {
              //int tmpBC_indicator = bcIndicator_of[cellIndex];
              // Info << "The boundary indicator is " << tmpBC_indicator << endl;
              if (tmpBC_indicator >= 0)
              {
                  forAll(primarySpecies_of, j)
                  {
                      Y_primary_of[j].boundaryFieldRef()[patchIDs][faces] = \
                      chem_bc_state[tmpBC_indicator].total_mobile.data[j] * 1000.;
                      // Info << "The concentration value is " << Y_primary_of[j].boundaryFieldRef()[patchIDs][faces] << endl;
                  }

                  if (pHValue_of.size() > 0)
                  {
                      // Info << "The pH value on the boundary is " << chem_aux_bc_output[tmpBC_indicator].pH << endl;
                      // Info << "The size of scalar field is " << pH_of[0].size() << endl;
                      pH_of[0].boundaryFieldRef()[patchIDs][faces] = chem_aux_bc_output[tmpBC_indicator].pH;
                  }
              }
            }
          }
        }
      }
  // Info << "Finish overwriting the concentrations of primary species on the boundaries." << endl;
}

void Foam::reactionAlquimia::overwriteChemicalYMultiphase(volScalarField& alpha1_of)
{
  // Info << "Start to overwrite primary species concentrations of each cell from CrunchFlow." << endl;
  forAll(mesh_of.cells(), i)
  {
      forAll(primarySpecies_of, j)
      {
	  if (chem_metadata.positivity.data[j] == 1 and chem_state[i].total_mobile.data[j]  <= 0)
	 {
	     Y_primary_of[j][i] = 1.0e-15;
	 }
	 else
	 {
             Y_primary_of[j][i] = chem_state[i].total_mobile.data[j] * 1000.0;
	 }
      }

  }
  //Overwrite the value on the boundary
  //Use mineral volume fraction to differentiate the reactive surface boundary from others (inlet/outlet)
  // Info << "Start to overwrite the boundary condition of primary species from CrunchFlow." << endl;
  forAll(mesh_of.boundaryMesh(), patchI)
  {
      string patchName = mesh_of.boundaryMesh()[patchI].name();
      // Info << "The boundary name in OpenFOAM " << patchName << endl;
      label patchIDs = mesh_of.boundaryMesh().findPatchID(patchName);
      const fvPatch& myPatches = mesh_of.boundary()[patchIDs];
      //Info << "The current patch size is " << myPatches.size() << endl;
      if (myPatches.size() < mesh_of.cells().size() and myPatches.type()== "patch")
      {
        forAll(myPatches, faces)
        {
            //Info << "The index of face is " << faces << endl;
            int tmpBC_indicator = bcIndicator_of.boundaryFieldRef()[patchIDs][faces];
            //Info << "The boundary indicator is " << tmpBC_indicator << endl;
            int tmpCountMineralV = 0;
            for (unsigned int j = 0; j < chem_sizes.num_minerals; ++j)
            {
                // Info << "The volume fraction is " << chem_state[cellIndex].mineral_volume_fraction.data[j] << endl;
                if (chem_bc_state[tmpBC_indicator].mineral_volume_fraction.data[j] > 0.0)
                {
                    tmpCountMineralV += 1;
                }
            }
            if (tmpCountMineralV == 0)
            {
              //int tmpBC_indicator = bcIndicator_of[cellIndex];
              // Info << "The boundary indicator is " << tmpBC_indicator << endl;
              if (tmpBC_indicator >= 0)
              {
                  forAll(primarySpecies_of, j)
                  {
                      Y_primary_of[j].boundaryFieldRef()[patchIDs][faces] = \
                      chem_bc_state[tmpBC_indicator].total_mobile.data[j] * 1000. * alpha1_of.boundaryFieldRef()[patchIDs][faces];
                      // Info << "The concentration value is " << Y_primary_of[j].boundaryFieldRef()[patchIDs][faces] << endl;
                  }

                  if (pHValue_of.size() > 0)
                  {
                      // Info << "The pH value on the boundary is " << chem_aux_bc_output[tmpBC_indicator].pH << endl;
                      // Info << "The size of scalar field is " << pH_of[0].size() << endl;
                      pH_of[0].boundaryFieldRef()[patchIDs][faces] = chem_aux_bc_output[tmpBC_indicator].pH * alpha1_of.boundaryFieldRef()[patchIDs][faces];
                  }
              }
            }
          }
        }
      }
  // Info << "Finish overwriting the concentrations of primary species on the boundaries." << endl;
}

/*Keep the primary species concentration positive with the standard from the array
of positivity.*/
void::Foam::reactionAlquimia::keepPrimaryConcPositive()
{
  // Info << "Convert the primary species concentration to the dimension used in CrunchFlow and keep them positive." << endl;
  forAll(mesh_of.cells(), i)
  {
      forAll(primarySpecies_of, j)
      {
          chem_state[i].total_mobile.data[j] = Y_primary_of[j][i] / 1000.;
          // Info << "The concentration is " << Y_primary_of[j][i] << endl;
          if (chem_metadata.positivity.data[j] == 1 and Y_primary_of[j][i] < 0.0)
          {
              //Info << "This " << i << "th primary species concentration must be positive." << endl;
              //Info << "The primary posivity is defined as " << chem_metadata.positivity.data[j] << endl;
              chem_state[i].total_mobile.data[j] = 1.0e-15;
          }
      }
  }
  // Info << "Overwrite the part of alquimia." << endl;
  // std::cin.get();
}

/*Function to do the reaction calculation for the single phase flow*/
void Foam::reactionAlquimia::reactionTimeStep(dimensionedScalar deltaT)
{
  // Info << "Solve reactions with CrunchFlow in single phase flow." << endl;
  // Info << "Operator split method is used here." << endl;
  // Info << "The time step is " << deltaT << endl;
  double timeStep;
  timeStep = deltaT.value();
  keepPrimaryConcPositive();
  forAll(mesh_of.cells(), i)
  {
      chem_state[i].water_density = waterDensity_of;
      chem_state[i].porosity = 1.0;
      chem.ReactionStepOperatorSplit(&chem_engine_alquimia, timeStep, \
                          &chem_properties[i], &chem_state[i], \
                          &chem_aux_data[i], &chem_status);
      //Info << "Solver information after doing the reaction calculation." << endl;
      bool converged = chem_status.converged;

      if (chem_status.error != 0 or converged != 1)
      {
          // status = chem_status.error;
          Info << "The converge status is " << converged << endl;
          Info << "The number of rhs evaluations is " << chem_status.num_rhs_evaluations << endl;
          Info << "The number of jacobian evaluations is " << chem_status.num_jacobian_evaluations << endl;
          Info << "The number of newton_iterations is " << chem_status.num_newton_iterations << endl;
          Info << "The solver for reaction failed." << endl;
          Info << chem_status.message << endl;
          std::abort();
      }
      chem.GetAuxiliaryOutput(&chem_engine_alquimia, &chem_properties[i], &chem_state[i], \
      &chem_aux_data[i], &chem_aux_output[i], &chem_status);
    }
    keepPrimaryConcPositive();
}

/*Function to do the reaction calculation for the multiphase phase flow condition*/
void Foam::reactionAlquimia::reactionTimeStepSaturation(dimensionedScalar deltaT, \
                              volScalarField& alpha1_of)
{
    // Info << "Solve reactions with CrunchFlow" << endl;
    // Info << "Operator split method is used here." << endl;
    // Info << "The time step is " << deltaT << endl;
    double timeStep;
    timeStep = deltaT.value();
    keepPrimaryConcPositive();
    forAll(mesh_of.cells(), i)
    {
        chem_state[i].water_density = waterDensity_of;
        chem_state[i].porosity = 1.0;
        if (alpha1_of[i] >= criticalValue_of)
        {
            chem_properties[i].saturation = alpha1_of[i];
            chem.ReactionStepOperatorSplit(&chem_engine_alquimia, timeStep, \
                                &chem_properties[i], &chem_state[i], \
                                &chem_aux_data[i], &chem_status);
            //Info << "Solver information after doing the reaction calculation." << endl;
            bool converged = chem_status.converged;

            if (chem_status.error != 0 or converged != 1)
            {
                // status = chem_status.error;
                Info << "The converge status is " << converged << endl;
                Info << "The number of rhs evaluations is " << chem_status.num_rhs_evaluations << endl;
                Info << "The number of jacobian evaluations is " << chem_status.num_jacobian_evaluations << endl;
                Info << "The number of newton_iterations is " << chem_status.num_newton_iterations << endl;
                Info << "The solver for reaction failed." << endl;
                Info << chem_status.message << endl;
                std::abort();
            }
            chem.GetAuxiliaryOutput(&chem_engine_alquimia, &chem_properties[i], &chem_state[i], \
            &chem_aux_data[i], &chem_aux_output[i], &chem_status);
        }
        else
        {
            forAll(primarySpecies_of, j)
            {
                chem_state[i].total_mobile.data[j] = 1.0e-20;
                chem_aux_output[i].mineral_reaction_rate.data[j] = 0.0;
            }
            for (int j = 0; j < chem_sizes.num_minerals; ++j)
            {
              chem_aux_output[i].mineral_saturation_index.data[j] = 1;
            }
            chem_aux_output[i].pH = 7.0;
        }
    }
}

/*Function to get the reaction rate for the dynamic mesh from the reaction calculation.*/
void Foam::reactionAlquimia::reactionTimeStepMovementofFace(PtrList<dimensionedScalar>& molarVolume, \
                                                volScalarField& movementRate)
{
    // Info << "Start to calculate the movement because of reactions." << endl;
    dimensionedScalar One("One",dimensionSet(0,-3,-1,0,1,0,0), 1.0);
    dimensionedScalar Volume("Volume", dimensionSet(0, 3, 0, 0, 0, 0, 0), 1.0);
    dimensionedScalar Area("Area", dimensionSet(0, 2, 0, 0, 0, 0, 0), 1.0);
    // reactionTimeStepSaturation(deltaT, alpha1_of);
    forAll(mesh_of.cells(), i)
    {
        dimensionedScalar tmpRateSum("tmpRateSum", dimensionSet(0, 3, -1, 0, 0, 0, 0), 0.0);
        dimensionedScalar tmpArea("tmpArea", dimensionSet(0, 2, 0, 0, 0, 0, 0), 0.0);
        // double tmpRateSum = 0.0;
        // double tmpArea = 0.0;
        forAll(mineralPhases_of, j)
        {
            double tmpReactionRate = fabs(chem_aux_output[i].mineral_reaction_rate.data[j]);
            double tmpSpecificArea = chem_state[i].mineral_specific_surface_area.data[j];
            tmpRateSum = tmpRateSum + (fabs(chem_aux_output[i].mineral_reaction_rate.data[j]))* One * mesh_of.V()[i] * Volume / molarVolume[j];
            tmpArea = tmpArea + chem_state[i].mineral_specific_surface_area.data[j] * mesh_of.V()[i] * Area;
            // tmpArea += chem_state[i].mineral_specific_surface_area.data[j] * mesh_of.V()[i];
        }
        if (tmpArea.value() != 0)
        {
          movementRate[i] = (tmpRateSum / tmpArea).value();
        }
    }
    //fvPatchField<scalar> tmpMPatch = movementRate.boundaryField()[patch().index()];
    /*Need to calculate the value of movement then transit the value from the cell center to the patch/face*/

}
