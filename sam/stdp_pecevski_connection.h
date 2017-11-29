/*
*  stdp_pecevski_connection.h
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
*  File based on stdp_connection.h in NEST.
*
*/

#ifndef STDP_PECEVSKI_CONNECTION_H
#define STDP_PECEVSKI_CONNECTION_H

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"
#include "numerics.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"

namespace sam
{
    using namespace nest;

    /**
	* @brief Rectangular-profile STDP based on Pecevski et al 2016[1].
	*
	* StdpPecevskiConnection is a variant of STDP with a simple rectangular
    * profile, where potentiation occurs if a pre-synaptic spike arrives at
    * the post-synaptic neuron within tau of the post-synaptic spike, and
    * depression occurs otherwise.
    *
    * Furthermore, a linearly-changing learning rate can be set, starting
    * at eta_0 and linearly changing to eta_final over learning_time ms.
    * If learning_time is set to -1, learning carries on forever at eta_0.
	*
	* <b>Parameters</b>
	*
	* The following parameters can be set in the status dictionary.
	* (default values and constraints are given in parentheses):
	*
	* <table>
	* <tr><th>name</th>                       <th>type</th>   <th>comment</th></tr>
	* <tr><td>\a tau</td>                     <td>double</td> <td>Potentiation window of STDP (20.0) [ms]</td></tr>
    * <tr><td>\a max_weight</td>              <td>double</td> <td>Maximum connection weight (4.0)</td></tr>
    * <tr><td>\a min_weight</td>              <td>double</td> <td>Minimum connection weight (0.0)</td></tr>
    * <tr><td>\a eta_0</td>                   <td>double</td> <td>Initial learning rate (0.1)</td></tr>
    * <tr><td>\a eta_final</td>               <td>double</td> <td>Final learning rate (0.0)</td></tr>
    * <tr><td>\a weight</td>                  <td>double</td> <td>Initial connection weight (1.0)</td></tr>
    * <tr><td>\a w_baseline</td>              <td>double</td> <td>Baseline weight offset (2.5 ln 0.2)</td></tr>
    * <tr><td>\a learning_time</td>           <td>long</td> <td>Learning duration in ms; -1 = forever (-1) [ms]</td></tr>
    * <tr><td>\a T</td>                       <td>long</td> <td>Scaling parameter in weight update.</td></tr>
    * <tr><td>\a tau_multiplier</td>                       <td>long</td> <td>Multiplier of tau that limits depressing effects.</td></tr>
    * </table>
	*
	* <i>Sends:</i> SpikeEvent
	*
	* <b>References</b>
	*
	* [1] Learning Probabilistic Inference through Spike-Timing-Dependent
	* Plasticity (2016) D. Peceveski and W. Maass, eNeuro
	*
	* @author  D'Amato
	*/
    template <typename targetidentifierT>
    class StdpPecevskiConnection : public Connection<targetidentifierT>
    {

    public:
        typedef CommonSynapseProperties CommonPropertiesType;
        typedef Connection<targetidentifierT> ConnectionBase;

        StdpPecevskiConnection();
        StdpPecevskiConnection( const StdpPecevskiConnection& );

        // Explicitly declare all methods inherited from the dependent base
        // ConnectionBase. This avoids explicit name prefixes in all places these
        // functions are used. Since ConnectionBase depends on the template parameter,
        // they are not automatically found in the base class.
        using ConnectionBase::get_delay_steps;
        using ConnectionBase::get_delay;
        using ConnectionBase::get_rport;
        using ConnectionBase::get_target;

        /**
         * Get all properties of this connection and put them into a dictionary.
         */
        void get_status(DictionaryDatum& d) const;

        /**
         * Set properties of this connection from the values given in dictionary.
         */
        void set_status(const DictionaryDatum& d, ConnectorModel& cm);

        /**
         * Send an event to the receiver of this connection.
         * \param e The event to send
         * \param t_lastspike Point in time of last spike sent.
         * \param cp common properties of all synapses (empty).
         */
        void send(Event& e, thread t, double t_lastspike, const CommonSynapseProperties& cp);

        class ConnTestDummyNode : public ConnTestDummyNodeBase
        {
        public:
            // Ensure proper overriding of overloaded virtual functions.
            // Return values from functions are ignored.
            using ConnTestDummyNodeBase::handles_test_event;
            port handles_test_event(SpikeEvent&, rport)
            {
                return invalid_port_;
            }
        };

        void check_connection(Node& s, Node& t, rport receptor_type, double t_lastspike, const CommonPropertiesType&)
        {
            ConnTestDummyNode dummy_target;
            ConnectionBase::check_connection_(dummy_target, s, t, receptor_type);

            t.register_stdp_connection(t_lastspike - get_delay());
        }

        void set_weight(double w)
        {
            weight_ = w;
        }

    private:
        double facilitate(double current_weight, double t) const
        {
            double delta_weight = numerics::expm1(-t_ * (w_baseline_ + current_weight));
            return std::max(std::min(current_weight + current_eta(t) * delta_weight, Wmax_), Wmin_);
        }

        double depress(double current_weight, double t) const
        {
            return std::max(std::min(current_weight - current_eta(t), Wmax_), Wmin_);
        }

        double current_eta(double t) const
        {
            if (learning_time_ == -1)
            {
                // If synapse learns forever, return the initial learning rate.
                return eta_0_;
            }

            if (learning_time_ > 0 && t <= learning_time_)
            {
                // If learning time is finite, return a weighted value of the initial and final
                // learning rates.
                return (1 - t / learning_time_) * eta_0_ + (t / learning_time_) * eta_final_;
            }

            // Otherwise, the synapse does not learn.
            return 0;
        }

        // Data members of each connection.
        double eta_0_;
        double eta_final_;
        double weight_;
        double tau_;
        double t_;
        double Wmax_;
        double Wmin_;
        double w_baseline_;
        double depress_tau_multiplier_;
        long learning_time_;
    };


/**
* Send an event to the receiver of this connection.
* \param e The event to send
* \param t The thread on which this connection is stored.
* \param t_lastspike Time point of last spike emitted
* \param cp Common properties object, containing the stdp parameters.
*/
    template < typename targetidentifierT >
    inline void StdpPecevskiConnection<targetidentifierT >::send(Event& e,
                                                         thread t,
                                                         double t_lastspike,
                                                         const CommonSynapseProperties& )
    {
        // Get this spike's timestamp.
        double t_spike = e.get_stamp().get_ms();

        // Use accessor functions (inherited from Connection) to obtain delay and target
        Node* target = get_target(t);
        double dendritic_delay = get_delay();

        // Get spike history in relevant range (t1, t2] from post-synaptic neuron
        std::deque<histentry>::iterator start;
        std::deque<histentry>::iterator finish;

        // For a new synapse, t_lastspike contains the point in time of the last
        // spike. So we initially read the
        // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
        // which increases the access counter for these entries.
        // At registration, all entries' access counters of
        // history[0, ..., t_last_spike - dendritic_delay] have been
        // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
        // details.
        target->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay, &start, &finish);

        // Plasticity due to post-synaptic spikes since last pre-synaptic spike.
        while (start != finish)
        {
            double dt = start->t_ + dendritic_delay - t_lastspike;

            if (dt >= 0 && dt <= tau_)
            {
                // Delta is within the potentiation window.
                weight_ = facilitate(weight_, start->t_ + dendritic_delay);
            }
            else if (dt <= depress_tau_multiplier_ * tau_)
            {
                // Delta is outside potentiation window.
                weight_ = depress(weight_, start->t_ + dendritic_delay);
            }

            ++start;
        }

        // NOTE: There is no depression from pre-synaptic spikes, since the rule explicitly states that updates
        // occur only on post-synaptic spikes.

        e.set_receiver(*target);
        e.set_weight(weight_);
        e.set_delay(get_delay_steps());
        e.set_rport(get_rport());
        e();
    }


    template <typename targetidentifierT>
    StdpPecevskiConnection<targetidentifierT>::StdpPecevskiConnection()
            : ConnectionBase(),
              eta_0_(0.1),
              eta_final_(0),
              weight_(1.0),
              tau_(20.0),
              t_(0.58),
              Wmax_(4.0),
              Wmin_(0.0),
              w_baseline_(-4.02359), // 2.5 * ln(0.2) (see [1]).
              depress_tau_multiplier_(5.0),
              learning_time_(-1) // i.e. learn forever.
    {

    }

    template <typename targetidentifierT>
    StdpPecevskiConnection<targetidentifierT>::StdpPecevskiConnection(const StdpPecevskiConnection<targetidentifierT>& rhs)
            : ConnectionBase(rhs),
              eta_0_(rhs.eta_0_),
              eta_final_(rhs.eta_final_),
              weight_(rhs.weight_),
              tau_(rhs.tau_),
              t_(rhs.t_),
              Wmax_(rhs.Wmax_),
              Wmin_(rhs.Wmin_),
              w_baseline_(rhs.w_baseline_),
              depress_tau_multiplier_(rhs.depress_tau_multiplier_),
              learning_time_(rhs.learning_time_)
    {

    }

    template <typename targetidentifierT>
    void StdpPecevskiConnection<targetidentifierT>::get_status(DictionaryDatum& d) const
    {
        ConnectionBase::get_status(d);
        def<double>(d, nest::names::weight, weight_);
        def<double>(d, nest::names::tau, tau_);
        def<double>(d, sam::names::t, t_);
        def<double>(d, sam::names::max_weight, Wmax_);
        def<double>(d, sam::names::min_weight, Wmin_);
        def<double>(d, sam::names::eta_0, eta_0_);
        def<double>(d, sam::names::eta_final, eta_final_);
        def<double>(d, sam::names::w_baseline, w_baseline_);
        def<long>(d, sam::names::learning_time, learning_time_);
        def<double>(d, sam::names::depress_multiplier, depress_tau_multiplier_);
        def<long>(d, nest::names::size_of, sizeof( *this ) );
    }

    template <typename targetidentifierT>
    void StdpPecevskiConnection<targetidentifierT>::set_status(const DictionaryDatum& d, ConnectorModel& cm)
    {
        ConnectionBase::set_status(d, cm);
        updateValue<double>(d, nest::names::weight, weight_);
        updateValue<double>(d, nest::names::tau, tau_);
        updateValue<double>(d, sam::names::t, t_);
        updateValue<double>(d, sam::names::max_weight, Wmax_);
        updateValue<double>(d, sam::names::min_weight, Wmin_);
        updateValue<double>(d, sam::names::eta_0, eta_0_);
        updateValue<double>(d, sam::names::eta_final, eta_final_);
        updateValue<double>(d, sam::names::w_baseline, w_baseline_);
        updateValue<long>(d, sam::names::learning_time, learning_time_);
        updateValue<double>(d, sam::names::depress_multiplier, depress_tau_multiplier_);

        // check if weight_ and Wmax_ has the same sign
        if (not(((weight_ >= 0 ) - (weight_ < 0)) == ((Wmax_ >= 0 ) - (Wmax_ < 0))))
        {
            throw BadProperty("Weight and Wmax must have same sign.");
        }

        if (!(std::abs(Wmin_) <= std::abs(Wmax_)))
        {
            throw BadProperty("Min weight must be smaller (absolute value) than max weight.");
        }

        if (eta_final_ < 0 || eta_0_ < 0)
        {
            throw BadProperty("Learning rates must be positive.");
        }

        if (tau_ <= 0)
        {
            throw BadProperty("Tau must be positive.");
        }

        if (depress_tau_multiplier_ < 1.0)
        {
            throw BadProperty("Tau multiplier for depression must be greater than or equal to 1.0.");
        }
    }
}

#endif // of #ifndef STDP_PECEVSKI_CONNECTION_H
