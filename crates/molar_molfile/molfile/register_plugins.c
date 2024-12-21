// Static registration of all needed plugins
#include "molfile_plugin.h"
#include "vmdplugin.h"
#include <stddef.h>
#include "bindings.h" // Brings plugin pointer variables

#define DECLARE_PLUGIN(name)\
    VMDPLUGIN_EXTERN int name##plugin_init();\
    VMDPLUGIN_EXTERN int name##plugin_register(void *v, vmdplugin_register_cb cb);\
    VMDPLUGIN_EXTERN int name##plugin_fini();\
    static molfile_plugin_t* name##_plugin_ptr = NULL;\
    static int name##_register_cb(void *v, vmdplugin_t *p){\
        name##_plugin_ptr = (molfile_plugin_t *)p;\
        return VMDPLUGIN_SUCCESS;\
    }


// Declare plugins
DECLARE_PLUGIN(pdb)
DECLARE_PLUGIN(xyz)
DECLARE_PLUGIN(dcd)

#define PLUGIN_GETTER(name)\
molfile_plugin_t* name##_get_plugin_ptr(){\
    if(!name##_plugin_ptr){\
        name##plugin_init();\
        name##plugin_register(NULL, &name##_register_cb);\
    }\
    return name##_plugin_ptr;\
}

PLUGIN_GETTER(pdb)
PLUGIN_GETTER(xyz)
PLUGIN_GETTER(dcd)
