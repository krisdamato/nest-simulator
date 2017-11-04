/*
 *  sammodule.cpp
 *
 *  This file is part of SAM, an extension of NEST.
 *
 *  Copyright (C) 2017 D'Amato
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "sammodule.h"

// Generated includes:
#include "config.h"

// include headers with your own stuff
#include "sam_names.h"
#include "srm_pecevski_alpha.h"
#include "stdp_pecevski_connection.h"

// Includes from nestkernel:
#include "connection_manager_impl.h"
#include "connector_model_impl.h"
#include "dynamicloader.h"
#include "exceptions.h"
#include "genericmodel.h"
#include "genericmodel_impl.h"
#include "kernel_manager.h"
#include "model.h"
#include "model_manager_impl.h"
#include "nestmodule.h"
#include "target_identifier.h"

// Includes from sli:
#include "booldatum.h"
#include "integerdatum.h"
#include "sliexceptions.h"
#include "tokenarray.h"

// -- Interface to dynamic module loader ---------------------------------------
#if defined( LTX_MODULE ) | defined( LINKED_MODULE )
sam::SamModule sammodule_LTX_mod;
#endif
// -- DynModule functions ------------------------------------------------------

sam::SamModule::SamModule()
{
#ifdef LINKED_MODULE
  // register this module at the dynamic loader
  // this is needed to allow for linking in this module at compile time
  // all registered modules will be initialized by the main app's dynamic loader
  nest::DynamicLoaderModule::registerLinkedModule(this);
#endif
}

sam::SamModule::~SamModule()
{

}

const std::string sam::SamModule::name(void) const
{
  return std::string("SAM module - Pecevski et al. 2016 models");
}

const std::string sam::SamModule::commandstring(void) const
{
  // Instruct the interpreter to load sammodule-init.sli
  return std::string( "(sammodule-init) run" );
}

//-------------------------------------------------------------------------------------

void sam::SamModule::init(SLIInterpreter* i)
{
  nest::kernel().model_manager.register_node_model<SrmPecevskiAlpha>("srm_pecevski_alpha");

  nest::kernel().model_manager.register_connection_model<StdpPecevskiConnection<nest::TargetIdentifierPtrRport> >("stdp_pecevski_synapse");
} // SamModule::init()
