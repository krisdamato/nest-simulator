/*
*  sam_names.h
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

#ifndef SAM_NAMES_H
#define SAM_NAMES_H

#include "name.h"

namespace sam
{
    /**
     * This namespace contains global Name objects. These can be used in
     * Node::get_status and Node::set_status to make data exchange
     * more efficient and consistent. Creating a Name from a std::string
     * is in O(log n), for n the number of Names already created. Using
     * predefined names should make data exchange much more efficient.
     */
    namespace names
    {
        extern const Name adaptive_threshold;
        extern const Name b_baseline;
        extern const Name e_0_exc;
        extern const Name e_0_inh;
        extern const Name eta_bias;
        extern const Name eta_0;
        extern const Name eta_final;
        extern const Name input_conductance;
        extern const Name learning_time;
        extern const Name max_bias;
        extern const Name min_bias;
        extern const Name r_m;
        extern const Name rect_exc;
        extern const Name rect_inh;
        extern const Name t;
        extern const Name tau_bias;
        extern const Name tau_exc;
        extern const Name tau_inh;
        extern const Name w_baseline;
    }
}

#endif //SAM_NAMES_H
