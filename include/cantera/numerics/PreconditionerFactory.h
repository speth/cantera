//! @file PreconditionerFactory.h

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#ifndef PRECONDITIONER_FACTORY_H
#define PRECONDITIONER_FACTORY_H

#include "cantera/numerics/PreconditionerBase.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

//! PreconditionerFactory created for use with the python interface to
//! create preconditioners based on a string type
class PreconditionerFactory : public Factory<PreconditionerBase>
{
public:
    static PreconditionerFactory* factory() {
        std::unique_lock<std::mutex> lock(precon_mutex);
        if (!s_factory) {
            s_factory = new PreconditionerFactory;
        }
        return s_factory;
    };

    //! Delete preconditioner factory
    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(precon_mutex);
        delete s_factory;
        s_factory = 0;
    };

    //! Create a new preconditioner by type name.
    /*!
     * @param preconType the type to be created.
     */
    virtual PreconditionerBase* newPreconditioner(const std::string& preconType);

private:
    static PreconditionerFactory* s_factory;
    static std::mutex precon_mutex;
    PreconditionerFactory();
};

//! Create a Preconditioner object of the specified type
inline PreconditionerBase* newPreconditioner(const std::string& precon)
{
    return PreconditionerFactory::factory()->newPreconditioner(precon);
};

}

#endif
