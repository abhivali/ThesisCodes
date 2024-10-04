#if !defined(CBonded_h)
#define CBonded_h

#include <string>
#include <fstream>
#include <sstream>

#include "IPluginContactModelV3_3_0.h"
#include "CBondParametersList.h"

/**
 * This class provides an implementation of IPluginContactModelV3_3_0
 */

 /******
  * That adds performs hertz mindlin contact force calculation
  * with added bonded particle physics
  *
  * This plugin reads a config file from the same directory as
  * it is stored in.  The file format is a single line containing
  *     The time at which bonds should be formed
  *
  * Followed by any number of lines defining bond types:
  *     type1:type_2 normal_stiffness
  *                  tangential_stiffness
  *                  normal_strength
  *                  shear_strength
  *                  bond_disk_radius;
  *
  * NOTE: This version of the hertz-midnlin & bonded particle code is
  * provided only as an example of how to create a contact model
  * plugin and is not updated on a regular basis.
  ******/


class CBonded : public NApiCm::IPluginContactModelV3_3_0
{
public:
    /**
     * Name of the preferences file to load bond information
     * from
     * Declare a string for the the preferences file name, custom property name
     */
    static const std::string GPU_FILE;
    static const std::string PREFS_FILE;

    // Declare a string for the the preferences file name, custom property name
    static const std::string BOND_STATUS;
    static const std::string BOND_PREFS;
    static const std::string BOND_NORMAL_FORCE;
    static const std::string BOND_TANGENTIAL_FORCE;
    static const std::string BOND_NORMAL_TORQUE;
    static const std::string BOND_TANGENTIAL_TORQUE;
    static const std::string BOND_SEPERATION;

    unsigned int iBOND_STATUS;
    unsigned int iBOND_PREFS;
    unsigned int iBOND_NORMAL_FORCE;
    unsigned int iBOND_TANGENTIAL_FORCE;
    unsigned int iBOND_NORMAL_TORQUE;
    unsigned int iBOND_TANGENTIAL_TORQUE;
    unsigned int iBOND_SEPERATION;

    /**
     * Constructor, does nothing
     */
    CBonded();

    /**
     * Destructor, does nothing
     */
    ~CBonded() override;

    /**
     * Sets prefFileName to PREFS_FILE
     *
     * Implementation of method from IPluginContactModelV3_3_0
     */
    void getPreferenceFileName(char prefFileName[NApi::FILE_PATH_MAX_LENGTH]) override;

    /**
     * Returns true to indicate plugin is thread safe
     *
     * Implementation of method from IPluginContactModelV3_3_0
     */
    bool isThreadSafe() override;

    /**
     * Returns true to indicate the plugin uses custom properties
     *
     * Implementation of method from IPluginContactModelV3_3_0
     */
    bool usesCustomProperties() override;

    /**
     * How the model should interact with chaining
     */
    NApi::EPluginModelType getModelType() override;
    NApi::EPluginExecutionChainPosition getExecutionChainPosition() override;

    /**
    * Returns void
    *
    * Give full path to simulation file or empty
    * string if none
    *
    * Implementation of method from IPluginContactModelV3_3_0
    */
    void setFilePath(const char simFile[]) override;

    /**
     * Gets the name of the OpenCL file
     *
     * Implementation of method from IPluginContactModelV3_3_0
     */
    void getGpuFileName(char nameFile[NApi::FILE_PATH_MAX_LENGTH]) override;

    /**
     * Setup function called when the plugin is loaded
     * Used to read preference file
     *
     * Implementation of method from IPluginContactModelV3_3_0
     */
    bool setup(NApiCore::IApiManager_1_0& apiManager,
        const char                 prefFile[],
        char                       customMsg[NApi::ERROR_MSG_MAX_LENGTH]) override;

    /**
     * Called at the beginning of the processing
     *
     * Implementation of method from IPluginContactModelV3_3_0
     */
    bool starting(NApiCore::IApiManager_1_0& apiManager, int numThreads) override;


    /**
     * Called at the end of the processing
     *
     * Implementation of method from IPluginContactModelV3_3_0
     */
    void stopping(NApiCore::IApiManager_1_0& apiManager) override;

    /**
     * Implementation of method from IPluginContactModelV3_3_0
     */
    NApi::ECalculateResult calculateForce(int                                            threadID,
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
        NCalcForceTypesV3_0_0::SContactResult& contactResults) override;


    /**
     * Indicates to EDEM the number and type of custom properties
     *
     * Implementation of method from IPluginContactModelV3_3_0
     */
    unsigned int getNumberOfRequiredProperties(const NApi::EPluginPropertyCategory category) override;
    /**
     * Define details of each of the custom properties
     *
     * Implementation of method from IPluginContactModelV3_3_0
     */
    bool getDetailsForProperty(unsigned int                    propertyIndex,
        NApi::EPluginPropertyCategory   category,
        char                            name[NApi::CUSTOM_PROP_MAX_NAME_LENGTH],
        NApi::EPluginPropertyDataTypes& dataType,
        unsigned int& numberOfElements,
        NApi::EPluginPropertyUnitTypes& unitType,
        char                            initValBuff[NApi::BUFF_SIZE]) override;

    unsigned int getPartPartContactParameterData(const char elem1Type[], const char elem2Type[], void* parameterData) override;

    unsigned int getSimulationParameterData(void* parameterData) override;

private:
    // Geometry API manager. Used to access geometry custome property with the configForTimeStep function
    //NApiCore::IGeometryManagerApi_1_0* m_geomMgr;

    double m_requestedBondTime;
    double m_bondingTimestep;
    int torque_feedback;
    int break_compression;
    double damping_coeff;
    double rotation_coeff;
    CBondParametersList m_bondParameters;
    double m_contactRadiusScale;

};

#endif
