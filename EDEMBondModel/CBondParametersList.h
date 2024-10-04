#if !defined(cbondparameterslist_h)
#define cbondparameterslist_h

#include <map>
#include <string>

typedef struct
{
    double      m_nElasticModulus;
    double      m_nPlasticModulus;
    double      m_nStiffnessRatio;
    double      m_nNormalStrength;
    double      m_nShearStrength;
    double      m_nBondDiskScale;
    double      m_nYieldStress;
}SBondParameters;

class CBondParametersList : public std::map<std::string, SBondParameters>
{
public:
    static const std::string JOIN_STRING;

    void addBondParameters(const std::string& key1,
        const std::string& key2,
        SBondParameters bondParameters);

    SBondParameters getBondParameters(const std::string& key1,
        const std::string& key2);
};

#endif