#ifndef IPLUGINPARTICLEBODYFORCEV3_1_0_H
#define IPLUGINPARTICLEBODYFORCEV3_1_0_H

/***************************************************************************/
/* This header file contains the V3.1.0 plugin particle body force API     */
/* definition.  Include this header and PluginParticleBodyForceCore.h into */
/* your plugin project then implement the methods from                     */
/* PluginParticleBodyForceCore.h and create a new class derived from       */
/* IPluginParticleBodyForceV3_1_0 that implements your desired             */
/* functionality.                                                          */
/***************************************************************************/

// Include ALL required headers.  Do not use forward declarations, this
// makes things easier on the end user
#include "ApiTypes.h"
#include "IApiManager_1_0.h"
#include "ICustomPropertyDataApi_1_0.h"
#include "IParticleManagerApi_1_0.h"
#include "IPluginParticleBodyForce.h"
#include "NExternalForceTypesV3_0_0.h"
#include "PluginConstants.h"

namespace NApiPbf
{
    /**
     * This interface contains all of the methods required to create a
     * particle body force plugin.  A new class should be created that
     * derives from this interface and implements all of its methods.
     * Additionally the methods from the PluginParticleBodyForceCore.h file
     * need to be implemented.
     *
     * NAME:              Particle Body Force Plugin API
     * VERSION:           3.1.0
     * CUSTOM PROPERTIES: Contact, Geometry, Particle, Simulation
     *
     * The methods required by the interface are:
     *
     * calculateForce(...)
     *     does the actual force calculation
     *
     * The following additional methods can be over-ridden to change
     * the default implementations:
     *
     * getPreferenceFileName(...)
     *     Returns the name of the config file used by this model (if any)
     *     Default: Indicates that no preference file is used
     *
     * isThreadSafe(...)
     *     Indicates that the plugin's externalForce(...) method is thread
     *     safe (returns true) or not (returns false).  A thread safe
     *     externalForce(...) can be called in parallel by threads running
     *     on multiple processors. This can speed up the simulation
     *     processing.
     *     Default: Indicates the plugin is not threadsafe
     *
     * usesCustomProperties(...)
     *     Indicates that the plugin wishes to register and/or receive custom
     *     property data.
     *     Default: Indicates the plugin does not use custom properties
     *
     * setFilePath(...)
     *     Provides the simulation file path to the plugin.
     *     Called within the starting(...) method.
     *
     * getGpuFileName(...)
     *     Provides the file name to the gpu plugin.
     *     Default:: empty
     *
     * setup(...)
     *     Performs any one-off setup and initialization just after the
     *     plugin is loaded.
     *     Default: Performs no work but returns true to indicate plugin
     *              loaded cleanly.
     *
     * starting(...)
     *     Called just before simulation starts to allow API handles to
     *     be acquired, temporary storage allocated, files opened etc.
     *     Default: Performs no work but returns true to indicate plugin
     *              is ready to start processing
     *
     * stopping()
     *     Called just after simulation stops to allow API handles to
     *     be released, temporary storage freed, files closed etc.
     *     Default: Does nothing
     *
     * getNumberOfRequiredProperties(...)
     *     Called by EDEM to query the number of custom properties required
     *     by this plugin to run.
     *     Default: Returns 0 to indicate no properties are required for
     *              any category
     *
     * getDetailsForProperty(...)
     *     Returns the details for the custom properties required for this plugin
     *     to run.
     *     Default: Returns false.
     *
     * If you need per plugin instance data simply add entries to your
     * plugin's class definition as you would with any other C++ class
     * definition.
     */
    class IPluginParticleBodyForceV3_1_0 : public IPluginParticleBodyForce
    {
    public:
        /**
         * Constructor, does nothing
         */
        IPluginParticleBodyForceV3_1_0() {}

        /**
         * Destructor, does nothing
         */
        virtual ~IPluginParticleBodyForceV3_1_0() {}

        /**
         * Retrieves the name of the config file used by the plugin.
         *
         * If the plugin does not need a config file then prefFileName
         * should be set to the empty string.
         *
         * @param prefFileName (RETURN VALUE)
         *                     A character array to be populated with the
         *                     config file name. This path is relative to
         *                     the directory the plugin is stored in.
         *                     EDEM will prepend the full directory the plugin
         *                     is stored in and pass it back to the setup method.
         */
        virtual void getPreferenceFileName(char prefFileName[NApi::FILE_PATH_MAX_LENGTH]) {prefFileName[0] = '\0';}

        /**
         * If the plugin implementations calculateForce() method is thread
         * safe then this method should return true.
         *
         * When a plug-in is thread-safe, edem allows multiple threads to
         * call the plug-in at the same time.  This can speed-up
         * calculations substantially on multi-processor machines.
         *
         * Thread safe programming requires a number of conventions and
         * restrictions to be followed.  If in doubt set this method to return
         * false.
         *
         * @return Bool to indicate if the calculateForce() method is thread
         *         safe
         */
        virtual bool isThreadSafe() {return false;}

        /**
         * Indicates whether the plugin wishes to register or receive custom
         * property data.
         *
         * @return Bool to indicate if custom properties are to be registered
         *         or should be supplied to the calculateForce(...) method.
         */
        virtual bool usesCustomProperties() {return false;}

        /**
         * Initializes the plugin by giving it the path to the simulation files.
         *
         * This method is called once, shortly after the plugin is first loaded,
         * within the call to function starting(...).
         *
         * IMPORTANT: Plugins should not cache API handles in this
         * method.  See the starting(...) and stopping(...) methods.
         *
         * @param simFile Full path to simulation file or empty
         *                 string if none
         * @return void
         */
        virtual void setFilePath(const char simFile[]) {;}

        /**
         * Initializes the plugin by giving it the name of the cl file for GPU solver.
         *
         * If empty, the model will not be supported on the GPU solver.
         * @param nameFile the name of the cl file. Do not include the extention .cl
         */
        virtual void getGpuFileName(char nameFile[NApi::FILE_PATH_MAX_LENGTH]) {;}

        /**
         * Initializes the plugin by giving it a chance to read any config
         * files, open temporary files, generate data structures or any other
         * one-off setup work.
         *
         * This method is called once, shortly after the plugin is first loaded.
         * If this method returns false EDEM will immediately delete the plugin
         * and an error message will be reported.
         *
         * IMPORTANT: Plugins should not cache API handles in this
         * method.  See the starting(...) and stopping(...) methods.
         *
         * @param   apiManager  The api manager for use by plugin models.
         * @param   prefFile    Full path to optional preferences file or empty
         *                      string if none.
         * @param   customMsg   (RETURN VALUE)
         *                      Character buffer to pass a custom error message to EDEM.
         * @return  bool        To say if setup was a success.
         */
        virtual bool setup(NApiCore::IApiManager_1_0& apiManager,
                           const char                 prefFile[],
                           char                       customMsg[NApi::ERROR_MSG_MAX_LENGTH]) {return true;}

        /**
         * Called to indicate processing is about to begin and the
         * model should allocate any temporary storage and retrieve any
         * file/api/socket handles it may need
         *
         * If the method returns false then processing will not start.
         *
         * IMPORTANT: Plugins should only retrieve API handles in this
         * method. API handles may change between one processing
         * run and another. Attempting to keep and re-use handles
         * will cause system instability.
         *
         * @param apiManager The api manager for use by plugin models
         * @param numThreads The number of threads this will be run with.
         * @return true if model is ready to start, else false
         */
        virtual bool starting(NApiCore::IApiManager_1_0& apiManager, int numThreads) {return true;}

        /**
         * Called to indicate processing is finished and that
         * the model should free any temporary storage and close/release
         * file/api/socket handles.
         *
         * The implementation must be able to handle this method being
         * called multiple times in a row without intervening calls
         * to starting.  This can occur when one or more loaded models
         * abort processing.
         *
         * IMPORTANT: Plugins must release all API handles in this
         * method. API handles may change between one processing
         * run and another. Attempting to keep and re-use handles
         * will cause system instability.
         */
        virtual void stopping(NApiCore::IApiManager_1_0& apiManager) {}

        /**
         * Use externalForce to add particle body forces (such as
         * electromagnetic or drag forces) to particles. This function
         * is called every single time step for every single particle.
         *
         * @param threadID                      The ID of the thread, in which this method is running,
         *                                      if multi-threaded (isThreadSafe must return true).
         * @param timeStepData                  Stores this time step's current time and the length of this time step.
         * @param particle                      Particle element, see NExternalForceTypesV3_0_0::SParticle for further details.
         * @param particleCustomProperties      Versioned interface providing access to the particle's
         *                                      custom property data and corresponding changeset.
         * @param simulationCustomProperties    Versioned interface providing access to custom
         *                                      property data and corresponding changeset for the simulation.
         * @param results                       (RETURN VALUE), see NExternalForceTypesV3_0_0::SResults for further details.
         * @return enum                         Value to indicate function result.
         */
        virtual NApi::ECalculateResult externalForce(
                                           int threadID,
                                           const NExternalForceTypesV3_0_0::STimeStepData& timeStepData,
                                           const NExternalForceTypesV3_0_0::SParticle& particle,
                                           NApiCore::ICustomPropertyDataApi_1_0* particleCustomProperties,
                                           NApiCore::ICustomPropertyDataApi_1_0* simulationCustomProperties,
                                           NExternalForceTypesV3_0_0::SResults& results) = 0;

        /**
         * Returns the number of custom properties this plugin wants to
         * register with the system for the supplied category.
         *
         * This version of the API supports the following property
         * categories:
         *     Contact Properties
         *     Geometry Properties
         *     Particle Properties
         *     Simulation Properties
         *
         * The method will be called once for each category at load time.
         * The implementation should return how many properties of that
         * category it wishes to register.
         *
         * If the plugin does not use custom properties this method should
         * return 0 for all categories.
         *
         * @param category The category of the custom property.
         * @return The number of custom properties the plugin wishes to
         *         register.
         */
        virtual unsigned int getNumberOfRequiredProperties(
                                 const NApi::EPluginPropertyCategory category) {return 0;}

        /**
         * Retrieves details for a given property.  This method will be
         * called for each category for propertyIndex values
         * 0...(getNumberOfRequiredProperties(category) - 1) to retrieve
         * the details for that property from the plugin.  These properties
         * will then be registered with the system if they do not clash
         * with any existing properties.
         *
         * This version of the API supports the following property
         * categories:
         *     Contact Properties
         *     Geometry Properties
         *     Particle Properties
         *     Simulation Properties
         *
         * If the plugin does not use custom properties this method should
         * always return false.
         *
         * @param propertyIndex    The index of the property to retrieve data
         *                         for
         * @param category         The category of the custom property to return
         *                         details for.
         * @param name             (RETURN VALUE)
         *                         A CUSTOM_PROP_MAX_NAME_LENGTH char array
         *                         is supplied to be populated with the name
         *                         of the property
         * @param dataType         (RETURN VALUE)
         *                         The data type of the property should always
         *                         be set to eDouble
         * @param numberOfElements (RETURN VALUE)
         *                         The number of elements (min 1)
         * @param unitType         (RETURN VALUE)
         *                         The unit type of the property
         * @param initValBuff      (RETURN VALUE)
         *                         Delimited string with details of the initial values
         *                         for each of the properties elements
         *
         * @return bool to say if data exists for the property
         */
        virtual bool getDetailsForProperty(
                         unsigned int                    propertyIndex,
                         NApi::EPluginPropertyCategory   category,
                         char                            name[NApi::CUSTOM_PROP_MAX_NAME_LENGTH],
                         NApi::EPluginPropertyDataTypes& dataType,
                         unsigned int&                   numberOfElements,
                         NApi::EPluginPropertyUnitTypes& unitType,
                         char                            initValBuff[NApi::BUFF_SIZE]) {return false;}

        /**
         * Enables custom functionality to be carried out on a per timestep basis.
         * The function is called at the start of each timestep to enable per timestep
         * configuration to be carried out
         *
         * If the plugin does not use custom properties this method should
         * return 0 for all categories.
         *
         * @param simData Details of the associated simulation property details
         * @param particleManager Particle property manager
         * @param time The current simulation time
         */
        virtual void configForTimeStep(
                        NApiCore::ICustomPropertyDataApi_1_0* simData,
                        NApiCore::IParticleManagerApi_1_0* particleManager,
                        double time) {}
        
        /**
         * Gets the particle parameter data in a buffer format
         * @param particleType The name/type of the particle
         * @param parameterData A preallocated buffer of size PARAM_MAX_SIZE containing parameter values for the particle type
         * @return The actual size of the parameter data
         */
        virtual unsigned int getParticleParameterData(const char particleType[], void* parameterData) {return 0;}

        /**
         * Gets the simulation parameter data in a buffer format
         * @param parameterData A preallocated buffer of size PARAM_MAX_SIZE containing simulation parameter values
         * @return The actual size of the parameter data
         */
        virtual unsigned int getSimulationParameterData(void* parameterData) {return 0;}

        /**
         * Process particles that were marked for additional processing in externalForce call.
         * Particle modifications are not allowed.
         * @param threadID The ID of the thread, in which this method is running,
                           if multi-threaded (isThreadSafe must return true).
         * @param particleOfInterestId One of the particle IDs that requires additional processing.
         */
        virtual void processParticleOfInterest(int threadID, int particleOfInterestId) {}

    };
};

#endif 

