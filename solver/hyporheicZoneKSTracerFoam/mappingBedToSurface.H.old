
        //map the subsurface top scalar to the surface bed bottom
        Info << "Map the concentration from sediment to riverbed." << endl;
        surfaceBottomFvPatchScalarField& tempCSurfp =
              dynamic_cast<surfaceBottomFvPatchScalarField&>
                      (C_Surf.boundaryFieldRef()[surfaceBottomPatchID]);


        tempCSurfp.CFromSub() =
          bedToSurfaceInterpolate.faceInterpolate<scalar>
          (
             C_Sub.boundaryField()[subSurfaceTopPatchID]
          );
          
        //map the subsurface saturation scalar to the surface bed bottom
        Info << "Construct the normal vector for the surface." << endl;
        const surfaceVectorField Sn_sub(meshSubsurface.Sf()/meshSubsurface.magSf());
        Info << "Get the vector from the sediment." << endl;
        const vectorField& Sn_sub_f = Sn_sub.boundaryField()[subSurfaceTopPatchID];

        // Info << "Map the chemical concentration to the surface flow boundary." << endl;
        // forAll(solutionSpecies, i)
        // {
        //     volScalarField& tmpYSurf1 = speciesMixture.Y(i);
        //     volScalarField& tmpYSub1 = geochemicalReactions.Y(i);
        //     chesurfaceBottomFvPatchScalarField& tmpSurfYiPatch = 
        //         dynamic_cast<chesurfaceBottomFvPatchScalarField&>(tmpYSurf1.boundaryFieldRef()[surfaceBottomPatchID]);
            
        //     tmpSurfYiPatch.YiFromSub() = 
        //         bedToSurfaceInterpolate.faceInterpolate<scalar>
        //         (
        //             tmpYSub1.boundaryField()[subSurfaceTopPatchID]
        //         );
        // }

        //map the subsurface top velocity to the surface bed bottom
        //Using Darcy velocity for consistency in fluid flux across 
        //the interface. Since USub is the pore-water velocity,
        //it needs to be multiplied by porosity to get Darcy velocity.  
        //May not work if the upwelling/downwelling is significant!
        //Only use the normla component may help because the only thing matters
        //is the normal velocity for flux.
        Info << "Map the subsurface velocity to the surface water." << endl;
        if(velocityNormalCouplingSwitch)
        {
            //Use the patch normal component.
            Info << "Use normal component." << endl;
            U.boundaryFieldRef()[surfaceBottomPatchID] ==
               bedToSurfaceInterpolate.faceInterpolate<vector>
               (
                  (Sn_sub_f & USub.boundaryField()[subSurfaceTopPatchID]*Sb.boundaryField()[subSurfaceTopPatchID])*Sn_sub_f
               );
        }
        else
        {
            //Use the whole vector.
            Info << "Use the whole velocity." << endl;
            U.boundaryFieldRef()[surfaceBottomPatchID] ==
               bedToSurfaceInterpolate.faceInterpolate<vector>
               (
                  USub.boundaryField()[subSurfaceTopPatchID]*Sb.boundaryField()[subSurfaceTopPatchID]
               );
        }
       Info << "Finish the velocity mapping from the subsurface." << endl;
        // surfaceBottomFvPatchScalarField& tempCSurfp =
        //       dynamic_cast<surfaceBottomFvPatchScalarField&>
        //               (CSurf.boundaryFieldRef()[surfaceBottomPatchID]);
        Info << "Map the saturation value from the top of the sediment to the bottom of the surface." << endl;
        surfaceBottomFvPatchScalarField& tmpAlpha1Surf = 
                dynamic_cast<surfaceBottomFvPatchScalarField&> (alpha1Tmp.boundaryFieldRef()[surfaceBottomPatchID]);
        Info << "Finish the tmporary pointer for surface saturation." << endl;
        tmpAlpha1Surf.CFromSub() = 
            bedToSurfaceInterpolate.faceInterpolate<scalar>
            (
                Sb.boundaryField()[subSurfaceTopPatchID]
            );
        Info << "Finish the transfer of the saturation from sediment to river." << endl;
        Info << "Copy the boundary value to the real field." << endl;
        alpha1.boundaryFieldRef()[surfaceBottomPatchID] = alpha1Tmp.boundaryFieldRef()[surfaceBottomPatchID];
         //map the subsurface top nonconservatice scalar to the surface bed bottom
        // forAll(surY, i)
        //   {
        //         if (surcomposition.active(i))
        //        {
        //                 //volScalarField& surYi =surY[i];
        //            YiFromSub[i].boundaryFieldRef()[surfaceBottomPatchID] ==
        //                         bedToSurfaceInterpolate.faceInterpolate<scalar>
        //                         (
        //                                 Y[i].boundaryField()[subSurfaceTopPatchID]
        //                         );

        //         //volScalarField& Yi =Y[i];
        //                 chesurfaceBottomFvPatchScalarField& tempYiSurfp =
        //             dynamic_cast<chesurfaceBottomFvPatchScalarField&>(surY[i].boundaryFieldRef()[surfaceBottomPatchID]);

        //                 tempYiSurfp.YiFromSub() =
        //                         bedToSurfaceInterpolate.faceInterpolate<scalar>
        //                         (
        //                                 Y[i].boundaryField()[subSurfaceTopPatchID]
        //                         );
        //       }
        //   }
                

        //Info << "U.boundaryField = " << U.boundaryField()[surfaceBottomPatchID] << endl; 
