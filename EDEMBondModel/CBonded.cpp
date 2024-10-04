//Last edited on 21 jun 2023
//Updated to Contact Model Version 3.3
//Has Damping
//HertzMindlin removed
//Has 'break in compression' as prefs file option
//Bond formation between multiple particles // added : Abhi
//uses elastic modulus and length of the bond to estimate the stiffness // added : Abhi
//uses stiffness ratio to differentiate between normal and shear stiffness // added : Abhi

#include "CBonded.h"
#include "ICustomPropertyManagerApi_1_0.h"
#include "HelpersV3_0_0.h"
#include "time.h"
#include <cstring>
#include <fstream>
#include <vector> // Abhi : included for the vector operations

using namespace std;
using namespace NApi;
using namespace NApiCore;
using namespace NApiCm;
using namespace NApiHelpersV3_0_0;

//Add the name for the preference file
const string CBonded::PREFS_FILE = "bonded_particle_prefs.txt";

//Add the name for the openCl file
const string CBonded::GPU_FILE = "BondedParticle_v3_1_0";

//Add the name for the custom property
const string CBonded::BOND_STATUS = "BondStatus";
const string CBonded::BOND_NORMAL_FORCE = "NormBondForce";
const string CBonded::BOND_TANGENTIAL_FORCE = "TangBondForce";
const string CBonded::BOND_NORMAL_TORQUE = "NormBondTorque";
const string CBonded::BOND_TANGENTIAL_TORQUE = "TangBondTorque";
const string CBonded::BOND_SEPERATION = "BondSeperation";
const string CBonded::BOND_PREFS = "BondProps";

CBonded::CBonded() :
    m_requestedBondTime(0.0)
{
    ;
}

CBonded::~CBonded()
{

}

void CBonded::getPreferenceFileName(char prefFileName[FILE_PATH_MAX_LENGTH])
{
    //Copy in the name of the preference file
    //Example code
    strncpy(prefFileName, PREFS_FILE.c_str(), FILE_PATH_MAX_LENGTH);
}

bool CBonded::isThreadSafe()
{
    //Thread safe
    return true;
}

bool CBonded::usesCustomProperties()
{
    //Does not use custom properties
    return true;
}

NApi::EPluginModelType CBonded::getModelType()
{
    return EPluginModelType::eOptional;
}

NApi::EPluginExecutionChainPosition CBonded::getExecutionChainPosition()
{
    return EPluginExecutionChainPosition::eAfterBasePos;
}

void CBonded::setFilePath(const char simFile[])
{
    ;
}

void CBonded::getGpuFileName(char nameFile[NApi::FILE_PATH_MAX_LENGTH])
{
    strncpy(nameFile, GPU_FILE.c_str(), NApi::FILE_PATH_MAX_LENGTH);
}

bool CBonded::setup(NApiCore::IApiManager_1_0& apiManager,
    const char                 prefFile[],
    char                       customMsg[NApi::ERROR_MSG_MAX_LENGTH])
{
    ifstream prefsFile(prefFile);

    if (!prefsFile)
    {
        return false;
    }
    else
    {
        srand(time(NULL));
        string description;

        //time at which the bonds form
        prefsFile >> description >> m_requestedBondTime;

        //torque feedback to be applied to particels as well as bonds? not applied in EDEM internal model
        prefsFile >> description >> torque_feedback;

        //want to break in compression as well as tension?
        prefsFile >> description >> break_compression;

        //bond damping value (0-1)
        prefsFile >> description >> damping_coeff;

        //bond rotation friction (typcially similar values to edem rolling friction)
        prefsFile >> description >> rotation_coeff;

        //scale of particle physical radius wrt contact radius
        //if you have a very large contact radius may only want to bond to particels within x% of radius
        prefsFile >> description >> m_contactRadiusScale;
        m_contactRadiusScale = (m_contactRadiusScale + 100) / 100;        //convert to usable format

        //read in bond type info to memory
        double nElasticModulus,
            nPlasticModulus,
            nStiffnessRatio,
            nNormalStrength,
            nShearStrength,
            nBondDiskScale,
            nYieldStress;
        string str;

        while (prefsFile)
        {
            prefsFile >> str
                >> description >> nElasticModulus
                >> description >> nPlasticModulus
                >> description >> nStiffnessRatio
                >> description >> nBondDiskScale
                >> description >> nYieldStress
                >> description >> nNormalStrength
                >> description >> nShearStrength;


            // Abhi: Added a vector to store the particle types and iterate the bond parameters
            //__________________________________________________________________________________________________
            vector<string> PTypeArray; // array to store string values of paticle types
            vector<string> PTypeArray2;

            string::size_type i = 0; // starting location to check string of particle types and initialised to 0
            string::size_type j(str.find(':')); // ending location (where the : delimiter is found)

            while (j != string::npos) // loop untill the end of string
            {
                string surf = str.substr(i, 1); // substring (start position, length)
                PTypeArray.push_back(surf);

                i = j + 1;
                j = str.find(':', i);

            }
            string surfEnd = str.substr(i); // adding the last Particle type to the array 
            PTypeArray.push_back(surfEnd);

            PTypeArray2 = PTypeArray;
            //__________________________________________________________________________________________________

            SBondParameters bondParams;

            bondParams.m_nElasticModulus = nElasticModulus;
            bondParams.m_nPlasticModulus = nPlasticModulus;
            bondParams.m_nStiffnessRatio = nStiffnessRatio;
            bondParams.m_nNormalStrength = nNormalStrength;
            bondParams.m_nShearStrength = nShearStrength;
            bondParams.m_nBondDiskScale = nBondDiskScale;
            bondParams.m_nYieldStress = nYieldStress; // Abhi: (ADDED yield stress to the bond properties)
            
            if (nNormalStrength < nYieldStress || nShearStrength < (nYieldStress*0.5))
            {
                return false;
            }
            if (nElasticModulus < nPlasticModulus)
            {
                return false;
            }

            // Abhi: assigning the bond parameters to all the possible combinations of particle types without repetition
            //__________________________________________________________________________________________________
            for (const string& surfA : PTypeArray)
            {
                for (const string& surfB : PTypeArray2)
                {

                    m_bondParameters.addBondParameters(surfA, surfB, bondParams);

                }

                PTypeArray2.erase(PTypeArray2.begin()); //deleting the first element to avoid repetition
            }
            //__________________________________________________________________________________________________
        }

    }

    return true;
}

bool CBonded::starting(NApiCore::IApiManager_1_0& apiManager, int numThreads)
{
    NApiCore::ICustomPropertyManagerApi_1_0* contactCustomPropertyManager
        = static_cast<NApiCore::ICustomPropertyManagerApi_1_0*>(apiManager.getApi(eContactCustomPropertyManager, 1, 0));

    iBOND_STATUS = contactCustomPropertyManager->getPropertyIndex(BOND_STATUS.c_str());
    iBOND_PREFS = contactCustomPropertyManager->getPropertyIndex(BOND_PREFS.c_str());;
    iBOND_NORMAL_FORCE = contactCustomPropertyManager->getPropertyIndex(BOND_NORMAL_FORCE.c_str());;
    iBOND_TANGENTIAL_FORCE = contactCustomPropertyManager->getPropertyIndex(BOND_TANGENTIAL_FORCE.c_str());;
    iBOND_NORMAL_TORQUE = contactCustomPropertyManager->getPropertyIndex(BOND_NORMAL_TORQUE.c_str());;
    iBOND_TANGENTIAL_TORQUE = contactCustomPropertyManager->getPropertyIndex(BOND_TANGENTIAL_TORQUE.c_str());;
    iBOND_SEPERATION = contactCustomPropertyManager->getPropertyIndex(BOND_SEPERATION.c_str());;

    return true;
}

void CBonded::stopping(NApiCore::IApiManager_1_0& apiManager)
{
    ;
}

NApi::ECalculateResult CBonded::calculateForce(int                                            threadID,
    const NCalcForceTypesV3_0_0::STimeStepData& timeStepData,
    const NCalcForceTypesV3_0_0::SDiscreteElement& element1,
    NApiCore::ICustomPropertyDataApi_1_0* elem1CustomProperties,
    const NCalcForceTypesV3_0_0::SDiscreteElement& element2,
    NApiCore::ICustomPropertyDataApi_1_0* elem2CustomProperties,
    NApiCore::ICustomPropertyDataApi_1_0* contactCustomProperties,
    NApiCore::ICustomPropertyDataApi_1_0* simulationCustomProperties,
    const NCalcForceTypesV3_0_0::SInteraction& interaction,
    const NCalcForceTypesV3_0_0::SContact& contact,
    NApiHelpersV3_0_0::CSimple3DVector& tangentialPhysicalOverlap,
    NCalcForceTypesV3_0_0::SContactResult& contactResults)
{
    //Setup bond status props
    const double* m_BondStatus = contactCustomProperties->getValue(iBOND_STATUS);
    double* m_BondStatusDelta = contactCustomProperties->getDelta(iBOND_STATUS);

    const double* m_BondSep = contactCustomProperties->getValue(iBOND_SEPERATION);
    double* m_BondSepDelta = contactCustomProperties->getDelta(iBOND_SEPERATION);

    //props for static values
    const double* m_BondPrefs = contactCustomProperties->getValue(iBOND_PREFS);

    //this is only recorded once, can later be used to check if bond is in compression or tension
    double* m_BondPrefsDelta = contactCustomProperties->getDelta(iBOND_PREFS);


    //check if in bond formation radius
    CSimple3DVector nPos1(element1.position);
    CSimple3DVector nPos2(element2.position);

    //single sphere only if using this!
    double nDistance = (nPos1 - nPos2).length();//center to center particle distance



    // Some temporary variables
    CSimple3DVector F_bond_n,
        F_bond_t,
        F_bond_n_damp,
        F_bond_t_damp,
        F_total_n,
        F_total_t,
        T_bond_1,
        T_bond_2,
        T_total_1,
        T_total_2,
        F_bond_np,
        FNormCurrent,
        FTangCurrent,
        TNormCurrent,
        TTangCurrent,
        F_bond_tp,
        dF_n_p,
        dF_t_p;

    /************************************************************************/
    /* Making / obtaining bonds from records                                */
    /************************************************************************/
    //only particle-particle


    if (true == element2.isSphere)
    {
        //are particles within bonding seperation distance? 
        //this allows you to have a very large contact radius but only bond to neighbouring elements
        //e.g only bond to particles within 10% of radius

        double nBondSeperation = nDistance - (element2.physicalRadius * m_contactRadiusScale) - (element1.physicalRadius * m_contactRadiusScale);
        
        //create bonds if:
        // form bond only if the timestep is greater than the bond formation step and below 3 timesteps 
        // this is ro avoid the bond formation throughout the simulation and limit to only few timesteps
        if ((timeStepData.time > m_requestedBondTime && 
            timeStepData.time < m_requestedBondTime + (3 * timeStepData.timeStep)) && 
            m_BondStatus[0] == 0.0 && 
            nBondSeperation < 0.0)
        {
            
            //change bond status to set that it exists
            m_BondStatusDelta[0] += 1.0;

            //record center to center seperation
            m_BondSepDelta[0] += nDistance;


            SBondParameters bondParams = m_bondParameters.getBondParameters(element1.type, element2.type);

            //normal stiffness in the elastic region
            m_BondPrefsDelta[0] = bondParams.m_nElasticModulus / (element1.physicalRadius + element2.physicalRadius);
            
            // shear stiffness in the elastic region
            m_BondPrefsDelta[1] = (bondParams.m_nElasticModulus / (element1.physicalRadius + element2.physicalRadius))/ bondParams.m_nStiffnessRatio; 

            //this creates a step function for bond strenghts, could be implemented using sin/cos for a smoother function

            //use smallest radius of two particles for contacts
            if (element1.physicalRadius <= element2.physicalRadius)
            {
                m_BondPrefsDelta[2] = element1.physicalRadius * bondParams.m_nBondDiskScale;
                //1 off calc for AREA
                m_BondPrefsDelta[5] = PI * element1.physicalRadius * element1.physicalRadius * bondParams.m_nBondDiskScale * bondParams.m_nBondDiskScale;
            }
            else
            {
                m_BondPrefsDelta[2] = element2.physicalRadius * bondParams.m_nBondDiskScale;
                //1 off calc for AREA if elem2 is smaller
                m_BondPrefsDelta[5] = PI * element2.physicalRadius * element2.physicalRadius * bondParams.m_nBondDiskScale * bondParams.m_nBondDiskScale;
            }

            m_BondPrefsDelta[3] = bondParams.m_nNormalStrength;
            m_BondPrefsDelta[4] = bondParams.m_nShearStrength;
            m_BondPrefsDelta[6] = bondParams.m_nYieldStress;

            //Normal stiffness in the plastic region
            m_BondPrefsDelta[7] = bondParams.m_nPlasticModulus / (element1.physicalRadius + element2.physicalRadius);

            //Shear stiffness in the plastic region considering incompressible material
            //stiffness ratio (kappa) corresponding to poissons ratio 0.5 is 3
            m_BondPrefsDelta[8] = (bondParams.m_nPlasticModulus / (element1.physicalRadius + element2.physicalRadius)) / 3.0;



        }
    }




    /************************************************************************/
    /* Predefining useful information on contact                            */
    /************************************************************************/

    //The unit vector from element 1 to the contact point */
    CSimple3DVector unitCPVect = contact.contactPoint - element1.position;
    unitCPVect.normalise();

    /* Calculate the relative velocities of the elements*/
    CSimple3DVector relVel = element1.velocityAtContactPoint - element2.velocityAtContactPoint;
    CSimple3DVector relVel_n = unitCPVect * unitCPVect.dot(relVel);
    CSimple3DVector relVel_t = relVel - relVel_n;

    double nEquivMass = element1.mass * element2.mass /
        (element1.mass + element2.mass); //equivMass(element1.mass, element2.mass);





/************************************************************************/
/* Bonded Particle Calculation                                          */
/************************************************************************/

    if (m_BondStatus[0] == 1.0)
    {

        const double* m_NForce = contactCustomProperties->getValue(iBOND_NORMAL_FORCE);
        double* m_NForceDelta = contactCustomProperties->getDelta(iBOND_NORMAL_FORCE);
        const double* m_TForce = contactCustomProperties->getValue(iBOND_TANGENTIAL_FORCE);
        double* m_TForceDelta = contactCustomProperties->getDelta(iBOND_TANGENTIAL_FORCE);
        const double* m_NTorque = contactCustomProperties->getValue(iBOND_NORMAL_TORQUE);
        double* m_NTorqueDelta = contactCustomProperties->getDelta(iBOND_NORMAL_TORQUE);
        const double* m_TTorque = contactCustomProperties->getValue(iBOND_TANGENTIAL_TORQUE);
        double* m_TTorqueDelta = contactCustomProperties->getDelta(iBOND_TANGENTIAL_TORQUE);

        double J = 0.5 * PI * m_BondPrefs[2] * m_BondPrefs[2] * m_BondPrefs[2] * m_BondPrefs[2];
        double I = J / 2.0;
        
        //Abhi: (Added)
        double YieldPoint;

        // yield point value edit :Abhi (added) store the initial value of the yield point taken from the text file
        YieldPoint = m_BondPrefs[6];

        //___________________________________________________________________________________________
        // Abhilash EDIT elastoplastic bond model
        // INITIALLY EXTRACT THE VALUES OF FORCE AND TORQUE TO TEST FOR THE STRESS CONDITION

        FNormCurrent = CSimple3DVector(m_NForce[0], m_NForce[1], m_NForce[2]);
        FTangCurrent = CSimple3DVector(m_TForce[0], m_TForce[1], m_TForce[2]);
        TNormCurrent = CSimple3DVector(m_NTorque[0], m_NTorque[1], m_NTorque[2]);
        TTangCurrent = CSimple3DVector(m_TTorque[0], m_TTorque[1], m_TTorque[2]);


        //check if bond is in compression or tension
        double sep = nDistance - m_BondSep[0];

        double StressCurrent;

        if (sep <= 0.0 && break_compression != 1)
        {
            //if in compression only consider the moments when breaking bonds
            StressCurrent = TTangCurrent.length() / I * m_BondPrefs[2];

        }
        else
        {
            //else bond will break in compression as well as tension
            StressCurrent = -FNormCurrent.length() / m_BondPrefs[5] + TTangCurrent.length() / I * m_BondPrefs[2];
        }

        double TorqueCurrent = FTangCurrent.length() / m_BondPrefs[5] + TNormCurrent.length() / J * m_BondPrefs[2];

        m_BondPrefsDelta[9] += StressCurrent - m_BondPrefs[9];
        m_BondPrefsDelta[10] += TorqueCurrent - m_BondPrefs[10];

        CSimple3DVector ndOverlap_n = relVel_n * timeStepData.timeStep;
        CSimple3DVector ndOverlap_t = relVel_t * timeStepData.timeStep;

        CSimple3DVector dF_n = -ndOverlap_n * m_BondPrefs[0] * m_BondPrefs[5];
        CSimple3DVector dF_t = -ndOverlap_t * m_BondPrefs[1] * m_BondPrefs[5];
        CSimple3DVector dF_n_P = -ndOverlap_n * m_BondPrefs[7] * m_BondPrefs[5];
        CSimple3DVector dF_t_P = -ndOverlap_t * m_BondPrefs[8] * m_BondPrefs[5];

        CSimple3DVector angle = (element1.angVel - element2.angVel) * timeStepData.timeStep;
        CSimple3DVector normalAngle = unitCPVect * unitCPVect.dot(angle);
        CSimple3DVector tangAngle = angle - normalAngle;

        CSimple3DVector dT_n = -normalAngle * J * m_BondPrefs[1];
        CSimple3DVector dT_t = -tangAngle * I * m_BondPrefs[0];

        CSimple3DVector dT_n_p = -normalAngle * J * m_BondPrefs[8];
        CSimple3DVector dT_t_p = -tangAngle * I * m_BondPrefs[7];


        if (fabs(StressCurrent) <= (m_BondPrefs[3]) && fabs(TorqueCurrent) <= (m_BondPrefs[4]))
        {
            // if bonds are not broken update the forces and torques

            F_bond_n = CSimple3DVector(m_NForce[0], m_NForce[1], m_NForce[2]);
            F_bond_t = CSimple3DVector(m_TForce[0], m_TForce[1], m_TForce[2]);
            CSimple3DVector T_n = CSimple3DVector(m_NTorque[0], m_NTorque[1], m_NTorque[2]);
            CSimple3DVector T_t = CSimple3DVector(m_TTorque[0], m_TTorque[1], m_TTorque[2]);

            
            if (fabs(StressCurrent) < YieldPoint && fabs(TorqueCurrent) < YieldPoint*0.5) // Abhi : (added)
            {

                F_bond_n = F_bond_n + dF_n;
                F_bond_t = F_bond_t + dF_t;
                T_n = T_n + dT_n;
                T_t = T_t + dT_t;

                m_NForceDelta[0] += dF_n.getX();
                m_NForceDelta[1] += dF_n.getY();
                m_NForceDelta[2] += dF_n.getZ();

                m_TForceDelta[0] += dF_t.getX();
                m_TForceDelta[1] += dF_t.getY();
                m_TForceDelta[2] += dF_t.getZ();

                // confirmation required to use torque
                m_NTorqueDelta[0] += dT_n.getX();
                m_NTorqueDelta[1] += dT_n.getY();
                m_NTorqueDelta[2] += dT_n.getZ();

                m_TTorqueDelta[0] += dT_t.getX();
                m_TTorqueDelta[1] += dT_t.getY();
                m_TTorqueDelta[2] += dT_t.getZ();

            }
            else // Abhi: (Added)
            {

                F_bond_n = F_bond_n + dF_n_P;
                F_bond_t = F_bond_t + dF_t_P;
                T_n = T_n + dT_n_p;
                T_t = T_t + dT_t_p;

                m_NForceDelta[0] += dF_n_P.getX();
                m_NForceDelta[1] += dF_n_P.getY();
                m_NForceDelta[2] += dF_n_P.getZ();

                m_TForceDelta[0] += dF_t_P.getX();
                m_TForceDelta[1] += dF_t_P.getY();
                m_TForceDelta[2] += dF_t_P.getZ();

                // confirmation required to use torque
                m_NTorqueDelta[0] += dT_n_p.getX();
                m_NTorqueDelta[1] += dT_n_p.getY();
                m_NTorqueDelta[2] += dT_n_p.getZ();

                m_TTorqueDelta[0] += dT_t_p.getX();
                m_TTorqueDelta[1] += dT_t_p.getY();
                m_TTorqueDelta[2] += dT_t_p.getZ();

                // abhi : adding the incremental value to the existing yield point
                // updating the yield point
                if (fabs(TorqueCurrent) >= YieldPoint * 0.5)
                {
                    if (fabs(StressCurrent) >= YieldPoint)
                    {
                        if (fabs(StressCurrent) < (2 * fabs(TorqueCurrent)))
                        {
                            m_BondPrefsDelta[6] = (2 * fabs(TorqueCurrent)) - YieldPoint;
                        }
                        else
                        {
                            m_BondPrefsDelta[6] = fabs(StressCurrent) - YieldPoint;
                        }
                    }
                    else
                    {
                        m_BondPrefsDelta[6] = fabs(StressCurrent) - YieldPoint;
                    }
                }
                else
                {
                    m_BondPrefsDelta[6] = fabs(StressCurrent) - YieldPoint;
                }
                
            }

            T_bond_1 = T_n + T_t;
            T_bond_2 = -T_bond_1;

        }
        else
        {
            F_bond_n = CSimple3DVector();
            F_bond_t = CSimple3DVector();
            m_BondStatusDelta[0] += 1.0;//set as broken
        }





        //__________________________________________________________________________________________
        /* // Original code


        CSimple3DVector ndOverlap_n = relVel_n * timeStepData.timeStep;
        CSimple3DVector ndOverlap_t = relVel_t * timeStepData.timeStep;

        CSimple3DVector dF_n = -ndOverlap_n * m_BondPrefs[0] * m_BondPrefs[5];
        CSimple3DVector dF_t = -ndOverlap_t * m_BondPrefs[1] * m_BondPrefs[5];

        CSimple3DVector angle = (element1.angVel - element2.angVel) * timeStepData.timeStep;
        CSimple3DVector normalAngle = unitCPVect * unitCPVect.dot(angle);
        CSimple3DVector tangAngle = angle - normalAngle;

        CSimple3DVector dT_n = - normalAngle * J * m_BondPrefs[1];
        CSimple3DVector dT_t = - tangAngle * I * m_BondPrefs[0];

        F_bond_n = CSimple3DVector(m_NForce[0],m_NForce[1],m_NForce[2]);
        F_bond_n = F_bond_n + dF_n;

        F_bond_t = CSimple3DVector(m_TForce[0],m_TForce[1],m_TForce[2]);
        F_bond_t = F_bond_t + dF_t;

        CSimple3DVector T_n = CSimple3DVector(m_NTorque[0],m_NTorque[1],m_NTorque[2]);
        T_n = T_n + dT_n;

        CSimple3DVector T_t = CSimple3DVector(m_TTorque[0],m_TTorque[1],m_TTorque[2]);
        T_t = T_t + dT_t;

        //check if bond is in compression or tension
        double sep = nDistance - m_BondSep[0];

        double maxStress;

        // check for damage
        if(sep <= 0.0 && break_compression != 1)
        {
            //if in compression only consider the moments when breaking bonds
            maxStress = T_t.length() / I * m_BondPrefs[2];

        }
        else
        {
            //else bond will break in compression as well as tension
            maxStress = -F_bond_n.length() / m_BondPrefs[5] + T_t.length() / I * m_BondPrefs[2];
        }

        double maxTorque = F_bond_t.length() / m_BondPrefs[5] + T_n.length() / J * m_BondPrefs[2];


        if(fabs(maxStress) <= (m_BondPrefs[3]) && fabs(maxTorque) <= (m_BondPrefs[4]))
        {
            //if bonds are not broken
            //update forces on bonds
            //N.B normal forces are still added in compression
            //BUT are not used in fracture calcs unless specified in prefs file
            //only when in tension (-ve normal)
            T_bond_1 = T_n + T_t;
            T_bond_2 = -T_bond_1;

            m_NForceDelta[0] += dF_n.getX();
            m_NForceDelta[1] += dF_n.getY();
            m_NForceDelta[2] += dF_n.getZ();

            m_TForceDelta[0] += dF_t.getX();
            m_TForceDelta[1] += dF_t.getY();
            m_TForceDelta[2] += dF_t.getZ();

            m_NTorqueDelta[0] += dT_n.getX();
            m_NTorqueDelta[1] += dT_n.getY();
            m_NTorqueDelta[2] += dT_n.getZ();

            m_TTorqueDelta[0] += dT_t.getX();
            m_TTorqueDelta[1] += dT_t.getY();
            m_TTorqueDelta[2] += dT_t.getZ();
        }
        else
        {
            F_bond_n = CSimple3DVector();
            F_bond_t = CSimple3DVector();
            m_BondStatusDelta[0] += 1.0;//set as broken
        }

        */

    }

    /****************************************/
    /* Normal and Tangential Bond Force Damping */
    /****************************************/
    //Simple equation: damping coeff * stiffness per unit area * area * rel_vel -
    //This is not used due to difficulty of setting damping coefficient

    //F_bond_n_damp = relVel_n *damping_coeff * m_BondPrefs[0] * m_BondPrefs[5];
    //F_bond_t_damp = relVel_t *damping_coeff * m_BondPrefs[1] * m_BondPrefs[5];

    //Refrence for damping equations implemented below:
    //Y. Guo, J. Curtis, C. Wassgren, W. Ketterhagan, B. Hancock. Granular shear flows of flexible rod-like particles.
    //Proceedings of Powders & Grains 2013, 8-12 July 2013, Sydney, Australia.
    
    // introducing plastic or elastic condition by comparing with recorded yield point
    // bondprefs[6] is yield point and is always positive but the current stress in bondprefs[9] is made absolute
    if (fabs(m_BondPrefs[9]) < m_BondPrefs[6]) 
    {
        F_bond_n_damp = relVel_n * damping_coeff * pow((m_BondPrefs[0] * m_BondPrefs[5] * 2.0 * nEquivMass), 0.5);
        F_bond_t_damp = relVel_t * damping_coeff * pow((m_BondPrefs[1] * m_BondPrefs[5] * 2.0 * nEquivMass), 0.5);
    }
    else 
    {
        F_bond_n_damp = relVel_n * damping_coeff * pow((m_BondPrefs[7] * m_BondPrefs[5] * 2.0 * nEquivMass), 0.5);
        F_bond_t_damp = relVel_t * damping_coeff * pow((m_BondPrefs[8] * m_BondPrefs[5] * 2.0 * nEquivMass), 0.5);
    }
    

    /****************************************/
    /* Fill in parameters we were passed in */
    /****************************************/
    F_total_n = F_bond_n - F_bond_n_damp;
    F_total_t = F_bond_t - F_bond_t_damp;

    if (torque_feedback == 1)
    {
        T_total_1 += T_bond_1;
        T_total_2 += T_bond_2;
    }

    contactResults.normalForce += F_total_n;
    contactResults.tangentialForce += F_total_t;

    contactResults.additionalTorque1 += T_total_1;
    contactResults.additionalTorque2 += T_total_2;

    /****************************************/
    /* Bond rotational 'friction' based on particle rolling fricton */
    //Not damping but simple implementation to help energy dissipation
    //values typically can range between 0 and 1 (1 is high damping) but better range is 0-0.1, could also have value higher than 1.
    /****************************************/

    //see Hertz-Mindling code for EDEM - rolling friction equation
    if (!isZero(element1.angVel.lengthSquared()))
    {
        CSimple3DVector torque1 = element1.angVel;
        torque1.normalise();
        double coM1toContactPointDistance = (contact.contactPoint - element1.CoM).length();
        torque1 *= -F_total_n.length() * coM1toContactPointDistance * rotation_coeff;
        contactResults.additionalTorque1 += torque1;
        contactResults.usAdditionalTorque1 += torque1;
    }

    if (!isZero(element2.angVel.lengthSquared()))
    {
        CSimple3DVector torque2 = element2.angVel;
        torque2.normalise();
        double coM2toContactPointDistance = (contact.contactPoint - element2.CoM).length();
        torque2 *= -F_total_n.length() * coM2toContactPointDistance * rotation_coeff;
        contactResults.additionalTorque2 += torque2;
        contactResults.usAdditionalTorque2 += torque2;
    }

    //Alternative method for damping moments below
    //Implementation needs to be updated take into account multi-sphere particles
    //implementation below has damping coeff between 10 and 0 not 1 and 0, error in implementation? 
    //Use friction above for damping instead of method below:

    //CSimple3DVector rel_ang_vel  = (angVel1 - angVel2);
    //CSimple3DVector normal_rel_ang_vel = unitCPVect * unitCPVect.dot(rel_ang_vel);
    //CSimple3DVector tang_rel_ang_vel = rel_ang_vel - normal_rel_ang_vel;

    //assume sphere for MOI calcs
    //CSimple3DVector torque1 = normal_rel_ang_vel * rotation_coeff;
    //torque1 *= -pow((elem1MoIX * 2.0 * J * m_BondPrefs[1]) ,0.5);

    //CSimple3DVector torque2 = tang_rel_ang_vel * rotation_coeff;
    //torque2 *= -pow((elem1MoIX * 2.0 * I * m_BondPrefs[0]) ,0.5);

    //calculatedElem1AdditionalTorqueX += torque1.dx() + torque2.dx();
    //calculatedElem1AdditionalTorqueY += torque1.dy() + torque2.dy();
    //calculatedElem1AdditionalTorqueZ += torque1.dz() + torque2.dz();

       //calculatedElem2AdditionalTorqueX += -torque1.dx() - torque2.dx();
    //calculatedElem2AdditionalTorqueY += -torque1.dy() - torque2.dy();
    //calculatedElem2AdditionalTorqueZ += -torque1.dz() - torque2.dz();

//printf("PP torque1 %d - %e,%e,%e\n", elem1Id, calculatedElem1AdditionalTorqueX, calculatedElem1AdditionalTorqueY, calculatedElem1AdditionalTorqueZ);
//printf("PP torque2 %d - %e,%e,%e\n", elem2Id, calculatedElem2AdditionalTorqueX, calculatedElem2AdditionalTorqueY, calculatedElem2AdditionalTorqueZ);

    return eSuccess;
}



unsigned int CBonded::getNumberOfRequiredProperties(const EPluginPropertyCategory category)
{
    if (eContact == category)
    {
        return 7;
    }
    else
    {
        return 0;
    }
}

bool CBonded::getDetailsForProperty(unsigned int                    propertyIndex,
    NApi::EPluginPropertyCategory   category,
    char                            name[NApi::CUSTOM_PROP_MAX_NAME_LENGTH],
    NApi::EPluginPropertyDataTypes& dataType,
    unsigned int& numberOfElements,
    NApi::EPluginPropertyUnitTypes& unitType,
    char                            initValBuff[NApi::BUFF_SIZE])
{
    // Define the characteristic of the properties: name, number of elements, unit
    if (0 == propertyIndex &&
        eContact == category)
    {
        strncpy(name, BOND_STATUS.c_str(), NApi::CUSTOM_PROP_MAX_NAME_LENGTH);
        dataType = eDouble;
        numberOfElements = 1;
        unitType = eNone;

        // Initialization of the property values in this case 0.0 - only one component
        std::ostringstream oss;
        oss << 0.0;
        strncpy(initValBuff, oss.str().c_str(), NApi::BUFF_SIZE);
        return true;
    }
    if (1 == propertyIndex &&
        eContact == category)
    {
        strncpy(name, BOND_NORMAL_FORCE.c_str(), NApi::CUSTOM_PROP_MAX_NAME_LENGTH);
        dataType = eDouble;
        numberOfElements = 3;
        unitType = eForce;

        // Initialization of the property values in this case 0.0 for all components - X,Y,Z
        std::ostringstream oss;
        oss << 0.0 << NApi::delim() << 0.0 << NApi::delim() << 0.0;
        strncpy(initValBuff, oss.str().c_str(), NApi::BUFF_SIZE);
        return true;
    }
    if (2 == propertyIndex &&
        eContact == category)
    {
        strncpy(name, BOND_TANGENTIAL_FORCE.c_str(), NApi::CUSTOM_PROP_MAX_NAME_LENGTH);
        dataType = eDouble;
        numberOfElements = 3;
        unitType = eForce;

        // Initialization of the property values in this case 0.0 for all components - X,Y,Z
        std::ostringstream oss;
        oss << 0.0 << NApi::delim() << 0.0 << NApi::delim() << 0.0;
        strncpy(initValBuff, oss.str().c_str(), NApi::BUFF_SIZE);
        return true;
    }
    if (3 == propertyIndex &&
        eContact == category)
    {
        strncpy(name, BOND_NORMAL_TORQUE.c_str(), NApi::CUSTOM_PROP_MAX_NAME_LENGTH);
        dataType = eDouble;
        numberOfElements = 3;
        unitType = eForce;

        // Initialization of the property values in this case 0.0 for all components - X,Y,Z
        std::ostringstream oss;
        oss << 0.0 << NApi::delim() << 0.0 << NApi::delim() << 0.0;
        strncpy(initValBuff, oss.str().c_str(), NApi::BUFF_SIZE);
        return true;
    }
    if (4 == propertyIndex &&
        eContact == category)
    {
        strncpy(name, BOND_TANGENTIAL_TORQUE.c_str(), NApi::CUSTOM_PROP_MAX_NAME_LENGTH);
        dataType = eDouble;
        numberOfElements = 3;
        unitType = eForce;

        // Initialization of the property values in this case 0.0 for all components - X,Y,Z
        std::ostringstream oss;
        oss << 0.0 << NApi::delim() << 0.0 << NApi::delim() << 0.0;
        strncpy(initValBuff, oss.str().c_str(), NApi::BUFF_SIZE);
        return true;
    }
    if (5 == propertyIndex &&
        eContact == category)
    {
        strncpy(name, BOND_PREFS.c_str(), NApi::CUSTOM_PROP_MAX_NAME_LENGTH);
        dataType = eDouble;
        numberOfElements = 11; // Abhi (Added two more parameters)
        unitType = eNone;

        // Initialization of the property values in this case 0.0 for all components - 10 in total
        std::ostringstream oss;
        oss << 0.0 << NApi::delim() << 0.0 << NApi::delim() << 0.0 << NApi::delim() << 0.0 << NApi::delim() << 0.0
            << NApi::delim() << 0.0 << NApi::delim() << 0.0 << NApi::delim() << 0.0 << NApi::delim() << 0.0 
            << NApi::delim() << 0.0 << NApi::delim() << 0.0; // Abhi (Added two more parameters)
        strncpy(initValBuff, oss.str().c_str(), NApi::BUFF_SIZE);
        return true;
    }
    if (6 == propertyIndex &&
        eContact == category)
    {
        strncpy(name, BOND_SEPERATION.c_str(), NApi::CUSTOM_PROP_MAX_NAME_LENGTH);
        dataType = eDouble;
        numberOfElements = 1;
        unitType = eLength;

        // Initialization of the property values in this case 0.0 - only one component
        std::ostringstream oss;
        oss << 0.0;
        strncpy(initValBuff, oss.str().c_str(), NApi::BUFF_SIZE);
        return true;
    }
    else
    {
        return false;
    }
}

unsigned int CBonded::getPartPartContactParameterData(const char elem1Type[], const char elem2Type[], void* parameterData)
{
    SBondParameters* contactParams = reinterpret_cast<SBondParameters*>(parameterData);
    *contactParams = m_bondParameters.getBondParameters(elem1Type, elem2Type);
    return sizeof(SBondParameters);
}

typedef struct
{
    double m_requestedBondTime;
    int torque_feedback;
    int break_compression;
    double damping_coeff;
    double rotation_coeff;
    double m_contactRadiusScale;
} SSimulationParameterData;

unsigned int CBonded::getSimulationParameterData(void* parameterData)
{
    SSimulationParameterData* simParams = reinterpret_cast<SSimulationParameterData*>(parameterData);
    simParams->m_requestedBondTime = m_requestedBondTime;
    simParams->torque_feedback = torque_feedback;
    simParams->break_compression = break_compression;
    simParams->damping_coeff = damping_coeff;
    simParams->rotation_coeff = rotation_coeff;
    simParams->m_contactRadiusScale = m_contactRadiusScale;
    return sizeof(SSimulationParameterData);
}