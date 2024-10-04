#if !defined(pluginexternalforcecore_h)
#define pluginexternalforcecore_h

/***************************************************************************/
/* This header file contains methods that must be implemented by           */
/* all V2.0.0 and later particle body force plugins.                       */
/***************************************************************************/

#include "PluginConstants.h"

namespace NApiPbf
{
    class IPluginParticleBodyForce;
};

/**
 * This method should return a newly allocated instance of the plugin
 * each time it is called.
 *
 * TO BE IMPLEMENTED BY END USER IN EACH PARTICLE BODY FORCE PLUGIN
 *
 * @return A newly allocated plugin instance
 */
EXPORT_MACRO NApiPbf::IPluginParticleBodyForce* GETPBFINSTANCE();

/**
 * This method should de-allocate the supplied instance.
 *
 * A plugin will only ever be supplied instances to
 * de-allocate that it created via it's GETPBFINSTANCE() method.
 *
 * TO BE IMPLEMENTED BY END USER IN EACH PARTICLE BODY FORCE PLUGIN
 *
 * @param instance The instance to be released
 */
EXPORT_MACRO void RELEASEPBFINSTANCE(NApiPbf::IPluginParticleBodyForce* instance);

/**
 * Returns the version of the plugin interface implemented
 * by this plugin file.
 *
 * TO BE IMPLEMENTED BY END USER
 *
 * The version number must be packed in to a single
 * 32 bit value.  Consider the following example code:
 *
 *  static const int INTERFACE_VERSION_MAJOR = 0x02;
 *  static const int INTERFACE_VERSION_MINOR = 0x00;
 *  static const int INTERFACE_VERSION_PATCH = 0x00;
 *
 *  return (INTERFACE_VERSION_MAJOR << 16 |
 *          INTERFACE_VERSION_MINOR << 8 |
 *          INTERFACE_VERSION_PATCH);
 *
 * @return Interface version packed in to a single 32bit value
 */
EXPORT_MACRO int GETEFINTERFACEVERSION();

#endif
