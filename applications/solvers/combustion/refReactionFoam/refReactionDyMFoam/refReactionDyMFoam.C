/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    flameletReactionDyMFoam

Description
    Solver for combustion with chemical reactions using density based
    thermodynamics package.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
//#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "basicSprayCloud.H"
#include "psiFlameletModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    //createTime.H: 
    // - declares 'runTime' object from Foam::Time class
    // - 'runTime' object reads file 'controlDict' in system
    #include "createTime.H"
    //createMesh.H: 
    // - prints out runTime for creating mesh
    // - declare 'mesh' object from Foam::fvMesh
    //#include "createMesh.H"
    #include "createDynamicFvMesh.H"
    //createControl.H:
    // - switch for pressure correction algorithm: PISO-PIMPLE-SIMPLE
    #include "createControl.H"
    #include "createControls.H"
    //createTimeControls.H:
    // - read time step control from controlDict through 'runTime' object
    // - declares and initializes the global variables: adjustTimeStep, maxCo, and maxDeltaT
    // - this is included in createControls.H
    //#include "createTimeControls.H"
    //createFields.H:
    // - exists in the same directory
    //
    #include "createFields.H" 
    #include "createFieldRefs.H"
    #include "createRhoUf.H"
    #include "createFvOptions.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // - readTimeControls.H is included in readControls.H
        //#include "readTimeControls.H"
        #include "readControls.H"

        #include "trackMaxScalars.H"
        #include "computeFuelMassFraction.H"
        {
            // Store divrhoU from the previous time-step/mesh for the correctPhi
            volScalarField divrhoU
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );

            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"

            runTime++;

            Info<< "Time = " << runTime.timeName() << nl << endl;

            // Store momentum to set rhoUf for introduced faces.
            volVectorField rhoU("rhoU", rho*U);

            // Do any mesh changes
            mesh.update();

            if (mesh.changing() && correctPhi)
            {
                // Calculate absolute flux from the mapped surface velocity
                phi = mesh.Sf() & rhoUf;

                #include "correctPhi.H"

                // Make the fluxes relative to the mesh-motion
                fvc::makeRelative(phi, rho, U);
            }
        }

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }

        parcels.evolve();

        #include "runRIFsolver.H"

        #include "rhoEqn.H"

        //- updateSh(): Update heat release rate quantity (enthalpy quantity) from RIF solution
        //- This should be done before using 'combustion->Sh()'.
        //- updateSh() function can't be called within Sh() function which is const type function.
        combustion->updateSh();
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "ZiEqn.H"
            #include "ZvarEqn.H"
            // species field is transported when RIF is inactive
            #include "YEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
        #include "updateScalars.H"

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

