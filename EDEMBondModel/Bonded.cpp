#include "PluginContactModelCore.h"
#include "CBonded.h"

using namespace NApiCm;

/******
 * NOTE: In this empty project example, it will ask EDEM 
 * to get a plugin instance called CContactModel. 
 * If you want to change the name of this class, 
 * search and replace "CContactModel" for all 
 * files in the project
 ******/

EXPORT_MACRO IPluginContactModel* GETCMINSTANCE()
{
	/*** Get an instance for CContactModel ***/
    return new CBonded();
}

EXPORT_MACRO void RELEASECMINSTANCE(IPluginContactModel* plugin)
{
    //Delete this instance, do not modify
    if (nullptr != plugin)
    {
        delete plugin;
    }
}

EXPORT_MACRO int GETCMINTERFACEVERSION()
{
	//Set the version number of the plugin, do not modify
    static const int INTERFACE_VERSION_MAJOR = 0x03;
    static const int INTERFACE_VERSION_MINOR = 0x03;
    static const int INTERFACE_VERSION_PATCH = 0x00;

    return (INTERFACE_VERSION_MAJOR << 16 |
            INTERFACE_VERSION_MINOR << 8 |
            INTERFACE_VERSION_PATCH);
}

