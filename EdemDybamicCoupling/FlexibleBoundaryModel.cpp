// Copyright DEM Solutions
// Mark Cook
// April 2009

// Modified Bill Edwards
// April 2013

// Updated to use IEDEMCoupling.h (EDEM 2017.2)
// May 2017

// Abhi thesis
// Updated for the implementation of coupling interface
// for flexible boundary condition of Electrode model
 

// EDEM Coupling include files
#define _CRT_SECURE_NO_WARNINGS
#include "IEDEMCoupling.h"

// Additional libraries 
#include <iostream>
#include <time.h>
#include <fstream>

using namespace std;
using namespace NApiEDEM;

class CGeometry // A simple geometry class
{
public:

    // ID
    int         id;
    
    // Geometry properties
    double      mass;
    C3x3Matrix	momentOfInertia;
	C3x3Matrix	initialMomentOfInertia;

    // External forces
    C3dVector	force;
    C3dVector	torque;
    
    // Geometry accelerations & velocities
    C3dVector   acceleration;
    C3dVector   angularAcceleration;
    C3dVector   velocity;
    C3dVector   angularVelocity;
    
    // Geometry position
    C3dPoint    originalCenterOfMass;
    C3dPoint	centerOfMass;
    C3x3Matrix	orientation;
    C3dVector	totalTranslation;
} ;

// Define the ending simulation time
const double ENDTIME = 2200000;
// Define the number of time-steps EDEM runs for between data exchanges
const int TIME_STEP_RATIO = 2;
// This controls whether output is written to the console or not
const bool CONSOLEOUT = true;
// This controls whether output is written to a file or not
const bool FILEOUT = false;


int main()
{
    // Simulation time-step
    double dt;

    // Simulation settings
    double simTime = 0.0; // simulation start time will be set to this should be relevant to the saved timestep if not 0
    double couplingTime = 50; // the calculation of plate dynamics will start from this timestep
    double SpringInitiationTime = couplingTime + 5050; // pressure control starts from this time
    
    // Declare geometries
    CGeometry box;
    
    
    // Mass of the moving geometry
    box.mass = 3e7;
    double BoxCSArea = 2500; //In [m²]
    double startCM = 40; // in [m] point of start of the spring initiation global coordinate position

    // Pressure settings
    double CompactionPressure = 1e4; // in pa with proper scaling 
    double MaxPressure = 2.001e4; // in pa with proper scaling 
    
    bool applyForce = false;
    double AppliedZForce = 0.0;
    double AppliedDamp = 0.0;

    // constant pressure for the compaction phase
    double compactionInitialVel = 0.001;
    double CompactionForce = CompactionPressure * BoxCSArea; //later this value will be updated once the seperator reaches startCM position 
    double compactionDampCoeff = 0.1 * CompactionForce/compactionInitialVel;

    // stiffness estimation
    C3dVector Maxforce = C3dVector(0, 0, -(MaxPressure * BoxCSArea));
    C3dVector MaxforceDisplacement = C3dVector(0, 0, 10); // disp. values to compute stiffness
    double stiffness = abs(Maxforce.z()) / MaxforceDisplacement.z();
    double SpringForce;

    // Alternative damping method: non viscous damping
    double Dampingfraction = 0.7;
    double DampingForce;

    // velocity capping during Pressure control time
    double VelCap = 0.0009;

///////////////////////////// Coupling Initialisation /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
    
    IEDEMCoupling coupling; // Create an instance of the EDEM Coupling to use
    
    // Initialise the coupling
    if (!coupling.initialiseCoupling())
    {
        std::cout << "Can't intialise the EDEM Coupling Client" << endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "EDEM Coupling Client initialised" << endl << "Connecting to EDEM..." << endl;

    // Connect to EDEM
    if(!coupling.connectCoupling())
    {
        std::cout << "Could not connect to EDEM" << endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Connection to EDEM successful" << endl;


////////////////////////// Obtain Simulation Variables /////////////////////////////
////////////////////////////////////////////////////////////////////////////////////    

    // Set the EDEM time to zero
    
    // If you have your options set to load last time-step in the EDEM
    // Simulator it will cause problems with reloading the last saved 
    // time-step.

    // Either disable the option or save your EDEM deck at time-step zero
    // before simulating
    coupling.getEDEMTime(simTime);

    if (simTime >= SpringInitiationTime)
    {
        // if the simulation starts after the spring initiation
        // the starting point for the spring is not recorded so set it to 0
        coupling.setEDEMTime(0);
    }

    coupling.getEDEMTimeStep(dt);

    dt *= TIME_STEP_RATIO;
    
    // Get the geometry ID
    if(coupling.getGeometryId("TopPlate", box.id))
    {
        std::cout << "Found the geometry" << endl;
    }

    if (simTime > 0) //Initialize variable if continuing simulation
    {
        coupling.getGeometryForces(box.id, box.force, box.torque);
        coupling.getGeometryVelocity(box.id, box.velocity, box.angularVelocity);
        if (abs(box.force.z()) > abs(0.95 * CompactionForce))
        {
            applyForce = true;
        }
    }

    // Get the initial coupling position and orientation from EDEM
    coupling.getGeometryPosition(box.id, box.originalCenterOfMass, box.orientation);
    coupling.getGeometryCoM(box.id,box.centerOfMass);
    coupling.getGeometryTranslation(box.id, box.totalTranslation);
    coupling.getGeometryVelocity(box.id, box.velocity, box.angularVelocity);
    
    
    // Prepare data write out file if required
    if(FILEOUT == true)
    {
        //open file
        fstream DataOut( "Data_Output.csv", ios::out|ios::app );

        // Add date and time
        time_t rawtime;
        time ( &rawtime );
        DataOut << ctime (&rawtime) << endl;

        // Add headers
        DataOut	<< "Time,Force(X),Force(Y),Force(Z),"
                << "Vel(X),Vel(Y),Vel(Z),"
                << endl;

        DataOut.close();
    }


///////////////////////////////Main Simulation Loop ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

    while (simTime < ENDTIME)
    {
        simTime += dt;
        
        if (simTime > couplingTime)
        {
            // Get geometry forces
            coupling.getGeometryForces(box.id, box.force, box.torque);
            coupling.getGeometryTranslation(box.id, box.totalTranslation);
            
            if (simTime <= SpringInitiationTime)
            {
                AppliedZForce = CompactionForce;
                AppliedDamp = -box.velocity.z() * compactionDampCoeff;

                startCM = box.centerOfMass.z();
            }

            if (simTime > SpringInitiationTime) 
            {
                // spring force
                SpringForce = ((box.centerOfMass.z() - startCM) * stiffness);

                // with non-viscous damping 
                AppliedDamp = -1 * abs(box.force.z() - CompactionForce - SpringForce) * Dampingfraction * (box.velocity.z() / abs(box.velocity.z()));

                // calculate new compaction pressure to be applied on the geometry
                AppliedZForce = CompactionForce + SpringForce;
                
            }

            if (applyForce)
            {
                // Update force values
                box.force.setX(0); // Only 1 degree of freedom
                box.force.setY(0); // Only 1 degree of freedom
                box.force += C3dVector(0, 0, -AppliedZForce); // Add the compressive force
                box.force += C3dVector(0, 0, AppliedDamp); // Add damping force

                // Acceleration
                box.acceleration = box.force / box.mass;

                // Velocity
                box.velocity += box.acceleration * dt;

                // velocity capping
                if (abs(box.velocity.z()) > compactionInitialVel)
                {
                    box.velocity = C3dVector(0, 0, (box.velocity.z() / abs(box.velocity.z())) * VelCap);
                }

                // Position
                box.totalTranslation += box.velocity * dt;
                
                box.centerOfMass = box.originalCenterOfMass + box.totalTranslation;
            }

            else // Initial compression
            {
                box.velocity = C3dVector(0, 0, -compactionInitialVel);
                box.totalTranslation += box.velocity * dt;
                box.centerOfMass = box.originalCenterOfMass + box.totalTranslation;
                
                

                if (abs(box.force.z()) > abs(0.95 * CompactionForce))
                {
                    applyForce = true;
                }
            } 
        }//
               
//////////////////////// Update EDEM with Calculated Motion /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
        
        // Send motion data to EDEM 
        coupling.setGeometryMotion(box.id,
            box.totalTranslation,
            box.orientation,
            box.velocity,
            box.angularVelocity,
            dt); // Action time should be a multiple
                 // of the EDEM time-step

        // Tell EDEM to perform the simulation step and pause the simulation when done
        coupling.simulate(dt, ENDTIME);



//////////////////////////////////// Reporting //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
        if (CONSOLEOUT == true)
        {
            if (simTime <= couplingTime)
            {
                std::cout << "-----------Before Coupling------------" << endl;
                std::cout << "Time = " << simTime << endl;
                std::cout << "Force on the geometry - " << box.force.x() << " " << box.force.y() << " " << box.force.z() << endl;
                std::cout << "Geometry Translation - " << box.totalTranslation.x() << " " << box.totalTranslation.y() << " " << box.totalTranslation.z() << endl;
            }
            if (simTime > couplingTime && simTime <= SpringInitiationTime)
            {
                std::cout << "-----------compaction phase-----------" << endl;
                std::cout << "Time = " << simTime << endl;
                std::cout << "Force on the geometry - " << box.force.x() << " " << box.force.y() << " " << box.force.z() << endl;
                std::cout << "Compaction Force reorded - " << AppliedZForce << endl;
                std::cout << "Damping applied - " << AppliedDamp << endl;
                std::cout << "Velocity of geometry - " << box.velocity.x() << " " << box.velocity.y() << " " << box.velocity.z() << endl;
                std::cout << "Geometry Translation - " << box.totalTranslation.x() << " " << box.totalTranslation.y() << " " << box.totalTranslation.z() << endl;
                std::cout << "Geometry position - " << box.centerOfMass.x() << " " << box.centerOfMass.y() << " " << box.centerOfMass.z() << endl;

            }
            if (simTime > SpringInitiationTime)
            {
                std::cout << "-----------Spring-damper phase--------" << endl;
                std::cout << "Time = " << simTime << endl;
                std::cout << "Force on the geometry - " << box.force.x() << " " << box.force.y() << " " << box.force.z() << endl;
                std::cout << "Compaction Force reorded - " << AppliedZForce << endl;
                std::cout << "Damping applied - " << AppliedDamp << endl;
                std::cout << "Velocity of geometry - " << box.velocity.x() << " " << box.velocity.y() << " " << box.velocity.z() << endl;
                std::cout << "Geometry Translation - " << box.totalTranslation.x() << " " << box.totalTranslation.y() << " " << box.totalTranslation.z() << endl;
                std::cout << "Geometry position - " << box.centerOfMass.x() << " " << box.centerOfMass.y() << " " << box.centerOfMass.z() << endl;
            }

        }

        if (FILEOUT == true)
        {
            //open file
            fstream DataOut( "Data_Output.csv", ios::out|ios::app );

            std::cout << simTime << ",";

            std::cout << box.force.x() << "," << box.force.y() << "," << box.force.z() << ",";

            std::cout << box.velocity.x() << "," << box.velocity.y() << "," << box.velocity.z() << ",";

            //close file
            DataOut.close();
        }
    }

    return 0; // Program exited normally
}
