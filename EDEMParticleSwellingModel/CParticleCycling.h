#ifndef CREPLACEPBF_H
#define CREPLACEPBF_H

#include <string.h>
#include "Helpers.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "IPluginParticleBodyForceV2_3_0.h"
#include "IParticleManagerApi_1_2.h"
#include "ICustomPropertyManagerApi_1_0.h"


using namespace std;
using namespace NApi;
using namespace NApiCore;
using namespace NApiPbf;

/**
 * This class provides an implementation of IPluginParticleBodyForceV2_3_0
 * That adds and updates a replace flag property to each particle
 * in the system.
 */
class CParticleCycling : public NApiPbf::IPluginParticleBodyForceV2_3_0
{
public:

	/**
	* Name of the preferences file to load information from
	*/
	static const std::string PREFS_FILE;
    /**
    * Name of replace custom property
    */
	static const std::string MAX_SCALE;
	unsigned int iMAX_SCALE;

    // new custom property same as particle max scale
    static const std::string Cycle_Number;  // declare new property name Edit: Abhilash
    unsigned int iCycle_Number;             // variable to store the property value Edit: Abhilash
	
    /**
    * Constructor, does nothing
    */
    CParticleCycling();

    /**
    * Destructor, does nothing
    */
    virtual ~CParticleCycling();

    /**
    * Returns true to indicate plugin is thread safe
    *
    * Implementation of method from IPluginParticleBodyForceV2_3_0
    */
    virtual bool isThreadSafe();

    /**
    * Returns true to indicate the plugin uses custom properties
    *
    * Implementation of method from IPluginParticleBodyForceV2_3_0
    */
    virtual bool usesCustomProperties();

	virtual void getPreferenceFileName(char prefFileName[NApi::FILE_PATH_MAX_LENGTH]);

	virtual bool setup(NApiCore::IApiManager_1_0 & apiManager, const char prefFile[], char customMsg[NApi::ERROR_MSG_MAX_LENGTH]);



    /**
    * Implementation of method from IPluginParticleBodyForceV2_3_0
    */
    virtual NApi::ECalculateResult externalForce(
                                            int          threadId,
                                            double       time,
                                            double       timestep,
                                            int          id,
                                            const char   type[],
                                            double       mass,
                                            double       volume,
                                            double       density,
                                            unsigned int surfaces,
                                            double       posX,
                                            double       posY,
                                            double       posZ,
                                            double       velX,
                                            double       velY,
                                            double       velZ,
                                            double       angVelX,
                                            double       angVelY,
                                            double       angVelZ,
                                            const double orientation[9],
                                            NApiCore::ICustomPropertyDataApi_1_0* particlePropData,
                                            NApiCore::ICustomPropertyDataApi_1_0* simulationPropData,
                                            double&      calculatedForceX,
                                            double&      calculatedForceY,
                                            double&      calculatedForceZ,
                                            double&      calculatedTorqueX,
                                            double&      calculatedTorqueY,
                                            double&      calculatedTorqueZ);

    /**
    * Returns 1 to indicate the plugin wishes to register 1 property
    *
    * Implementation of method from IPluginParticleBodyForceV2_3_0
    */
    virtual unsigned int getNumberOfRequiredProperties(
                                const NApi::EPluginPropertyCategory category);

    /**
    * Returns details for our custom properties
    *
    * Implementation of method from IPluginParticleBodyForceV2_3_0
    */
	virtual bool getDetailsForProperty(
		unsigned int                    propertyIndex,
		NApi::EPluginPropertyCategory   category,
		char                            name[NApi::CUSTOM_PROP_MAX_NAME_LENGTH],
		NApi::EPluginPropertyDataTypes& dataType,
		unsigned int&                   numberOfElements,
		NApi::EPluginPropertyUnitTypes& unitType,
		char                            initValBuff[NApi::BUFF_SIZE]);


    virtual bool starting(NApiCore::IApiManager_1_0& apiManager, int numThreads);

private:
    /** particle manager. */
    NApiCore::IParticleManagerApi_1_0* m_particleMngr;
    double          m_maxcycles;
	vector<string>  m_particleNames;
	vector<double>  m_expansionType;
	vector<double>  m_expansionFinalScale;
	vector<double>  m_expansionStartTime;
    // declaring the new vectors  : edit:abhilash
    vector<double>  m_DwellTime;
};

#endif // CREPLACEPBF_H
