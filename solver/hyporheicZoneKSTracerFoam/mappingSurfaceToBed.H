
        //map the surface bottom pressure to the sediment region top patch
        Info << "Map the pressure surface pressure to the top of the sediment." << endl;
        if(flowSwitch) 
        {
           pSub.boundaryFieldRef()[subSurfaceTopPatchID] ==
             surfaceToBedInterpolate.faceInterpolate<scalar>
             (
                p.boundaryField()[surfaceBottomPatchID] 
             );
        }


        //map the surface bottom scalar to the sediment region top patch
        Info << "Map the tracer concentration to the top of sediments." << endl;
        subsurfaceTopFvPatchScalarField& tempCSubp =
              dynamic_cast<subsurfaceTopFvPatchScalarField&>(C_Sub.boundaryFieldRef()[subSurfaceTopPatchID]);

        tempCSubp.CFromSurf() =
            surfaceToBedInterpolate.faceInterpolate<scalar>
              (
                 C_Surf.boundaryField()[surfaceBottomPatchID]
              );
   
         //map the surface bottom nonconservative scalar to the sediment region top patch
        // forAll(solutionSpecies, i)
        // {
        //     //copy the conc of chemicals in subsurface
        //     volScalarField& tmpYSub = geochemicalReactions.Y(i);
        //     //copy the conc of chemicals in surface
        //     volScalarField& tmpYSurf = speciesMixture.Y(i);
        //     chesubsurfaceTopFvPatchScalarField& tmpSubYiPatch = 
        //         dynamic_cast<chesubsurfaceTopFvPatchScalarField&>(tmpYSub.boundaryFieldRef()[subSurfaceTopPatchID]);
            
        //     tmpSubYiPatch.YiFromSurf() = 
        //         surfaceToBedInterpolate.faceInterpolate<scalar>
        //         (
        //             tmpYSurf.boundaryField()[surfaceBottomPatchID]
        //         );
        // }
        // {
        //    forAll(Y, i)
        //      {
        //         if (composition.active(i))
        //         {
        //            //volScalarField& subYi =Y[i];

        //            YiFromSurf[i].boundaryFieldRef()[subSurfaceTopPatchID] ==
        //                  surfaceToBedInterpolate.faceInterpolate<scalar>
        //                  (
        //                  surY[i].boundaryField()[surfaceBottomPatchID]
        //                  );
        //            //forAll(surY, i)
        //          //{ 
        //        //if(surcomposition.active(i))
        //     //{       
        //       //volScalarField& surYi =surY[i]; 
        //       chesubsurfaceTopFvPatchScalarField& tempsubYip =
        //           dynamic_cast<chesubsurfaceTopFvPatchScalarField&>(Y[i].boundaryFieldRef()[subSurfaceTopPatchID]);

        //       tempsubYip.YiFromSurf() =
        //          surfaceToBedInterpolate.faceInterpolate<scalar>
        //          (
        //           surY[i].boundaryField()[surfaceBottomPatchID]
        //          );
        //        }
        //      }
        // }


      //map the saturation from the surface bottom saturation to the sediment region top patch 
    Info << "Map the saturation from the surface bottom to the top of the sediment top. Step 1:" << endl;
    subsurfaceTopFvPatchScalarField& tmpSbSubp = 
            dynamic_cast<subsurfaceTopFvPatchScalarField&> (Sb.boundaryFieldRef()[subSurfaceTopPatchID]);
    Info << "Map the saturation to the top of the sediment. Step 2:" << endl;

    tmpSbSubp.CFromSurf() = 
            surfaceToBedInterpolate.faceInterpolate<scalar>
            (
               alpha1.boundaryField()[surfaceBottomPatchID]
            );
    Info << "Finish the map." << endl;
    int topPatchNum = Sb.boundaryFieldRef()[subSurfaceTopPatchID].size();
    for (int i = 0; i < topPatchNum; ++i)
    {
        if (Sb.boundaryFieldRef()[subSurfaceTopPatchID][i] < 0.9)
        {
            Sb.boundaryFieldRef()[subSurfaceTopPatchID][i] = Sbmin.value();
        }
        else
        {
            Sb.boundaryFieldRef()[subSurfaceTopPatchID][i] = Sbmax.value();
        }
    } 

    // forAll(solutionSpeciesSurf, i)
    // {
    //     Info << "Map the chemical " << i << endl;
    //     volScalarField tmpYiSurf = YSurf[i];
    //     Info << "subsurface top." << endl;
    //     int nCellTop = Y_Sub[i].boundaryFieldRef()[subSurfaceTopPatchID].size();
    //     Info << "The number of cell on the top boundary:" << nCellTop << endl;
    //     chesubsurfaceTopFvPatchScalarField& tmpSubYiTop = 
    //         dynamic_cast<chesubsurfaceTopFvPatchScalarField&>(Y_Sub[i].boundaryFieldRef()[subSurfaceTopPatchID]);
    //     Info << "Change the concentration on the top." << endl;
    //     tmpSubYiTop.YiFromSurf() = 
    //         surfaceToBedInterpolate.faceInterpolate<scalar>
    //         (
    //             tmpYiSurf.boundaryField()[surfaceBottomPatchID]
    //         );
    // }
//   Info << "CSurf.boundaryField()[surfaceBottomPatchID] = " << CSurf.boundaryField()[surfaceBottomPatchID] << endl;
//   Info << "tempCSubp.CFromSurf = " << tempCSubp.CFromSurf() << endl;    
//   Info << "interpolate return = " <<  surfaceToBedInterpolate.faceInterpolate<scalar>(CSurf.boundaryField()[surfaceBottomPatchID]) << endl;



