/*
*  srm_pecevski_alpha.cpp
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
*  File based on poisson_dbl_exp_neuron.h in NEST::SPORE extension.
*
*/

#include "srm_pecevski_alpha.h"

#include "compose.hpp"
#include "doubledatum.h"
#include "dict.h"
#include "dictutils.h"
#include "exceptions.h"
#include "integerdatum.h"
#include "numerics.h"
#include "sam_defines.h"
#include "sam_names.h"
#include "universal_data_logger_impl.h"

#include "../nestkernel/logging_manager.h"


namespace nest
{
	/*
	* Recordables map of SrmPecevskiAlpha.
	*/
	template<>
	void RecordablesMap<sam::SrmPecevskiAlpha>::create()
	{
		// use standard names wherever you can for consistency!
		insert_(nest::names::V_m, &sam::SrmPecevskiAlpha::get_V_m_);
		insert_(nest::names::E_sfa, &sam::SrmPecevskiAlpha::get_E_sfa_);
	}
}

namespace sam
{
	using nest::BadProperty;

	/*
	* Recordables map instance.
	*/
	nest::RecordablesMap<SrmPecevskiAlpha> SrmPecevskiAlpha::recordablesMap_;


	//
	// SrmPecevskiAlpha::Parameters_ implementation.
	//
	SrmPecevskiAlpha::Parameters_::Parameters_():
    use_rect_psp_exc_(false),
    use_rect_psp_inh_(false),
    resistance_(1),
	epsilon_0_exc_(2.8),			// mv
	epsilon_0_inh_(2.8),			// mv
	tau_alpha_exc_(8.5),			// ms
	tau_alpha_inh_(8.5),			// ms
    tau_membrane_(10.0),            // ms
	dead_time_(1.0),				// ms
	dead_time_random_(false),		// ms
	dead_time_shape_(1l),
	with_reset_(true),
	c_1_(0),
	c_2_(1.238),
	c_3_(0.25),
	I_e_(0),						// pA
	t_ref_remaining_(0),			// ms
	bias_baseline_(-1.0),
    eta_bias_(0.1),
    tau_bias_(8.5),
    max_bias_(5.0),
    min_bias_(-30.0),
    t_(0.58)
	{

	}

	void SrmPecevskiAlpha::Parameters_::get(DictionaryDatum& d) const
	{
        def<bool>(d, sam::names::rect_exc, use_rect_psp_exc_);
        def<bool>(d, sam::names::rect_inh, use_rect_psp_inh_);
        def<double>(d, sam::names::r_m, resistance_);
		def<double>(d, nest::names::dead_time, dead_time_);
		def<double>(d, nest::names::dead_time_random, dead_time_random_);
		def<long>(d, nest::names::dead_time_shape, dead_time_shape_);
		def<double>(d, sam::names::e_0_exc, epsilon_0_exc_);
		def<double>(d, sam::names::e_0_inh, epsilon_0_inh_);
		def<double>(d, sam::names::tau_exc, tau_alpha_exc_);
        def<double>(d, sam::names::tau_inh, tau_alpha_inh_);
        def<double>(d, nest::names::tau_m, tau_membrane_);
		def<bool>(d, nest::names::with_reset, with_reset_);
		def<double>(d, nest::names::c_1, c_1_);
		def<double>(d, nest::names::c_2, c_2_);
		def<double>(d, nest::names::c_3, c_3_);
		def<double>(d, nest::names::I_e, I_e_);
		def<double>(d, nest::names::t_ref_remaining, t_ref_remaining_);
        def<double>(d, sam::names::b_baseline, bias_baseline_);
        def<double>(d, sam::names::eta_bias, eta_bias_);
        def<double>(d, sam::names::tau_bias, tau_bias_);
        def<double>(d, sam::names::max_bias, max_bias_);
        def<double>(d, sam::names::min_bias, min_bias_);
        def<double>(d, sam::names::t, t_);
	}

	void SrmPecevskiAlpha::Parameters_::set(const DictionaryDatum& d)
	{
        updateValue<bool>(d, sam::names::rect_exc, use_rect_psp_exc_);
        updateValue<bool>(d, sam::names::rect_inh, use_rect_psp_inh_);
        updateValue<double>(d, sam::names::r_m, resistance_);
		updateValue<double>(d, nest::names::dead_time, dead_time_);
		updateValue<double>(d, nest::names::dead_time_random, dead_time_random_);
		updateValue<long>(d, nest::names::dead_time_shape, dead_time_shape_);
		updateValue<double>(d, sam::names::e_0_exc, epsilon_0_exc_);
		updateValue<double>(d, sam::names::e_0_inh, epsilon_0_inh_);
		updateValue<double>(d, sam::names::tau_exc, tau_alpha_exc_);
		updateValue<double>(d, sam::names::tau_inh, tau_alpha_inh_);
        updateValue<double>(d, nest::names::tau_m, tau_membrane_);
		updateValue<bool>(d, nest::names::with_reset, with_reset_);
		updateValue<double>(d, nest::names::c_1, c_1_);
		updateValue<double>(d, nest::names::c_2, c_2_);
		updateValue<double>(d, nest::names::c_3, c_3_);
		updateValue<double>(d, nest::names::I_e, I_e_);
		updateValue<double>(d, nest::names::t_ref_remaining, t_ref_remaining_);
		updateValue<double>(d, sam::names::b_baseline, bias_baseline_);
        updateValue<double>(d, sam::names::eta_bias, eta_bias_);
        updateValue<double>(d, sam::names::tau_bias, tau_bias_);
        updateValue<double>(d, sam::names::max_bias, max_bias_);
        updateValue<double>(d, sam::names::min_bias, min_bias_);
        updateValue<double>(d, sam::names::t, t_);

		if (dead_time_ < 0.0)
		{
			throw BadProperty("Dead time must be >= 0.");
		}

		if (dead_time_shape_ < 1)
		{
			throw BadProperty("Dead time shape must be >= 1.");
		}

		if (tau_alpha_exc_ <= 0.0 || tau_alpha_inh_ <= 0.0)
		{
			throw BadProperty("All decay constants must be greater than 0.");
		}

        if (tau_membrane_ < 0.0)
        {
            throw BadProperty("Membrane time constant must be 0 or greater.");
        }

		if (epsilon_0_exc_ <= 0.0 || epsilon_0_inh_ <= 0.0)
		{
			throw BadProperty("All PSP absolute amplitudes bust be greater than 0.");
		}

        if (resistance_ <= 0)
        {
            throw BadProperty("Resistance must be greater than 0.");
        }

		if (c_3_< 0.0)
		{
			throw BadProperty("c_3 must be >= 0.");
		}

		if (t_ref_remaining_ < 0.0)
		{
			throw BadProperty("t_ref_remaining must be >= 0.");
		}

        if (eta_bias_ < 0)
        {
            throw BadProperty("eta_bias must be >= 0");
        }

        if (tau_bias_ <= 0)
        {
            throw BadProperty("tau_bias must be > 0");
        }

        if (max_bias_ < min_bias_)
        {
            throw BadProperty("max_bias must be greather than min_bias");
        }
	}

	//
	// SrmPecevskiAlpha::State_ implementation.
	//


	/**
	* Default constructor.
	*/
	SrmPecevskiAlpha::State_::State_():
		u_membrane_(0.0),
		input_current_(0.0),
		adaptive_threshold_(0.0),
        u_i_(0.0),
		r_(0)
	{

	}

	/**
	* State getter function.
	*/
	void SrmPecevskiAlpha::State_::get(DictionaryDatum& d, const Parameters_&) const
	{
		def<double>(d, nest::names::V_m, u_membrane_); // Membrane potential
		def<double>(d, sam::names::adaptive_threshold, adaptive_threshold_);
	}

	/**
	* Sate setter function.
	*/
	void SrmPecevskiAlpha::State_::set(const DictionaryDatum& d, const Parameters_&)
	{
		updateValue<double>(d, nest::names::V_m, u_membrane_);
		updateValue<double>(d, names::adaptive_threshold, adaptive_threshold_);
	}

	//
	// SrmPecevskiAlpha::Buffers_ implementation.
	//

	/**
	* Constructor.
	*/
	SrmPecevskiAlpha::Buffers_::Buffers_(SrmPecevskiAlpha& n)
		: logger_(n)
	{
	}

	/**
	* Constructor.
	*/
	SrmPecevskiAlpha::Buffers_::Buffers_(const Buffers_&, SrmPecevskiAlpha& n)
		: logger_(n)
	{
	}

	//
	// SrmPecevskiAlpha implementation.
	//

	/**
	* Default Constructor.
	*/
	SrmPecevskiAlpha::SrmPecevskiAlpha()
		: Archiving_Node(),
		P_(),
		S_(),
		B_(*this)
	{
		recordablesMap_.create();
	}

	/**
	* Copy Constructor.
	*/
	SrmPecevskiAlpha::SrmPecevskiAlpha(const SrmPecevskiAlpha& n)
		: Archiving_Node(n),
		P_(n.P_),
		S_(n.S_),
		B_(n.B_, *this)
	{
	}

	/**
	* Node state initialization.
	*/
	void SrmPecevskiAlpha::init_state_(const nest::Node& proto)
	{
		const SrmPecevskiAlpha& pr = downcast<SrmPecevskiAlpha>(proto);
		S_ = pr.S_;
		S_.r_ = nest::Time(nest::Time::ms(P_.t_ref_remaining_)).get_steps();
	}

	/**
	* Initialize the node's spike and current buffers.
	*/
	void SrmPecevskiAlpha::init_buffers_()
	{
		B_.exc_queue_.Clear();
		B_.inh_queue_.Clear();
		B_.currents_.clear(); //!< includes resize
		B_.logger_.reset(); //!< includes resize
	}

	/**
	* Calibrate the node.
	*/
	void SrmPecevskiAlpha::calibrate()
	{
		B_.logger_.init();

		V_.h_ = nest::Time::get_resolution().get_ms();
        V_.propagator_ = std::exp(-V_.h_ / P_.tau_membrane_);
		V_.rng_ = nest::kernel().rng_manager.get_rng(get_thread());

		V_.t_1 = 0.23196095298653444; // First solution of t exp(1-t) = 0.5

		if (P_.dead_time_ != 0 && P_.dead_time_ < V_.h_)
			P_.dead_time_ = V_.h_;

		// TauR specifies the length of the absolute refractory period as
		// a double in ms. The grid based iaf_psp_delta can only handle refractory
		// periods that are integer multiples of the computation step size (h).
		// To ensure consistency with the overall simulation scheme such conversion
		// should be carried out via objects of class nest::Time. The conversion
		// requires 2 steps:
		//
		//     1. A time object r is constructed defining the representation of
		//        TauR in tics. This representation is then converted to computation time
		//        steps again by a strategy defined by class nest::Time.
		//     2. The refractory time in units of steps is read out by get_steps(), a member
		//        function of class nest::Time.
		//
		// The definition of the refractory period of the SrmPecevskiAlpha is consistent
		// with the one of iaf_neuron_ps.
		//
		// Choosing a TauR that is not an integer multiple of the computation time
		// step h will lead to accurate (up to the resolution h) and self-consistent
		// results. However, a neuron model capable of operating with real valued spike
		// time may exhibit a different effective refractory time.
		if (P_.dead_time_random_)
		{
			// Choose dead time rate parameter such that mean equals dead_time
			V_.dt_rate_ = P_.dead_time_shape_ / P_.dead_time_;
			V_.gamma_dev_.set_order(P_.dead_time_shape_);
		}
		else
		{
			V_.DeadTimeCounts_ = nest::Time(nest::Time::ms(P_.dead_time_)).get_steps();
			assert(V_.DeadTimeCounts_ >= 0); // Since t_ref_ >= 0, this can only fail in error
		}
	}

	/*
	* Spike response kernel: half-alpha shape or rectangular (depending on user parameters).
	*/
	double SrmPecevskiAlpha::kernel(const double time_since_spike, const bool use_exc_kernel = true) const
	{
		const double& t = time_since_spike;
        const double& tau_alpha = use_exc_kernel ? P_.tau_alpha_exc_ : P_.tau_alpha_inh_;
        const double& epsilon = use_exc_kernel ? P_.epsilon_0_exc_ : P_.epsilon_0_inh_;
        bool use_rect_psp = use_exc_kernel ? P_.use_rect_psp_exc_ : P_.use_rect_psp_inh_;

        // Use the rectangular PSP if requested.
        if (use_rect_psp)
        {
            return (t >= 0 && t <= tau_alpha) ? epsilon / 2 : 0.0;
        }

        // Otherwise return the half-alpha PSP.
        return epsilon * ((t / tau_alpha + V_.t_1) * (std::exp(1 - (t / tau_alpha + V_.t_1))) - 0.5);
	}

	/*
	 * Sums up the PSPs from excitatory or inhibitory spikes.
	 */
	double SrmPecevskiAlpha::get_psp_sum(const nest::Time& now, const bool use_exc_psp)
	{
		double psp = 0;
        SpikeQueue& queue = use_exc_psp ? B_.exc_queue_ : B_.inh_queue_;

        for (SpikeQueue::IteratorType it = queue.Begin(); it != queue.End();)
        {
            // Get spike time and value.
            long spike_time = it->first;
            double amplitude = it->second;

            double delta = (now - nest::Time::step(spike_time)).get_ms();
            double this_psp = amplitude * kernel(delta, use_exc_psp);
            if (this_psp <= sam::effective_zero && delta > 0) // We'll remove spikes when they aren't affecting membrane voltage.
            {
                // Erase the spike, because we won't need it anymore.
                it = queue.EraseItemAt(it);
            }
            else
            {
                ++it;
            }

            psp += std::max(0., this_psp);
        }

		return use_exc_psp? psp : -psp;
	}

	/**
	* Update the node to the given time point.
	*/
	void SrmPecevskiAlpha::update(nest::Time const& origin, const long from, const long to)
	{
		assert(from < to);

		for (long lag = from; lag < to; ++lag)
		{
			nest::Time now = nest::Time::step(origin.get_steps() + lag);

            // Calculate PSP responses.
			double psp_exc = get_psp_sum(now, true);
			double psp_inh = get_psp_sum(now, false);

            // Update current-induced voltage.
            if (V_.propagator_ >= 0.05)
            {
                // If the membrane time constant is large enough, use low-pass filtering.
                S_.u_i_ = V_.propagator_ * S_.u_i_ + (1 - V_.propagator_) * P_.resistance_ * (S_.input_current_ + P_.I_e_);
            }
            else
            {
                // Otherwise, use a straightforward jump.
                S_.u_i_ = P_.resistance_ * (S_.input_current_ + P_.I_e_);
            }

            // Update total potential.
			S_.u_membrane_ = psp_exc + psp_inh + S_.u_i_;

            // Update intrinsic bias and clip.
			S_.adaptive_threshold_ -= P_.eta_bias_ * V_.h_ * 1e-3; // The 1e-3 is necessary since V_.h_ is in ms.
            S_.adaptive_threshold_ = std::max(std::min(S_.adaptive_threshold_, P_.max_bias_), P_.min_bias_);

			if (S_.r_ == 0)
			{
				// Neuron is not refractory.

				// Calculate instantaneous rate from transfer function:
				//     rate = c1 * u' + c2 * exp(c3 * u')
				double V_eff = S_.u_membrane_ + S_.adaptive_threshold_;

				double rate = (P_.c_1_ * V_eff + P_.c_2_ * std::exp(P_.c_3_ * V_eff));
				double spike_probability = -numerics::expm1(-rate * V_.h_ * 1e-3);
				long n_spikes = 0;

				if (rate > 0.0)
				{
					if (P_.dead_time_ > 0.0)
					{
						// Draw random number and compare to probability to have a spike
						if (V_.rng_->drand() <= spike_probability)
							n_spikes = 1;
					}
					else
					{
						// Draw Poisson random number of spikes
						V_.poisson_dev_.set_lambda(rate);
						n_spikes = V_.poisson_dev_.ldev(V_.rng_);
					}

					if (n_spikes > 0) // Is there a spike? Then set the new dead time.
					{
						// Set dead time interval according to parameters.
						if (P_.dead_time_random_)
						{
							S_.r_ = nest::Time(nest::Time::ms(V_.gamma_dev_(V_.rng_) / V_.dt_rate_)).get_steps();
						}
						else
							S_.r_ = V_.DeadTimeCounts_;

						// Set spike time.
                        set_spiketime(nest::Time::step(origin.get_steps() + lag + 1));

						// And send the spike event.
						nest::SpikeEvent se;
						se.set_multiplicity(n_spikes);
						nest::kernel().event_delivery_manager.send(*this, se, lag);

						// Reset the potential if applicable.
						if (P_.with_reset_)
						{
							B_.exc_queue_.Clear();
							B_.inh_queue_.Clear();

							S_.u_membrane_ = 0.0;
                            S_.u_i_ = 0.0;
						}

                        // Update intrinsic bias and clip.
						S_.adaptive_threshold_ += P_.eta_bias_ * P_.tau_bias_ *
                                std::exp(-P_.t_ * (S_.adaptive_threshold_ + P_.bias_baseline_));
                        S_.adaptive_threshold_ = std::max(std::min(S_.adaptive_threshold_, P_.max_bias_), P_.min_bias_);
                    }
				}
			}
			else // Neuron is within dead time
			{
				--S_.r_;
			}

			// Set new input current
			S_.input_current_ = B_.currents_.get_value(lag);

			// Voltage logging
			B_.logger_.record_data(origin.get_steps() + lag);
		}
	}

	/**
	* SpikeEvent handling.
	* @param e the event.
	*/
	void SrmPecevskiAlpha::handle(nest::SpikeEvent& e)
	{
		assert(e.get_delay() > 0);

		if (e.get_rport() == 0)
		{
			// Add spike to the queue.
			// Note: we need to compute the absolute number of steps since the beginning of simulation time.
			B_.exc_queue_.AddSpike(e.get_rel_delivery_steps(nest::Time()), e.get_weight() * e.get_multiplicity());
		}
		else if (e.get_rport() == 1)
		{
			// Add spike to the queue.
			// Note: we need to compute the absolute number of steps since the beginning of simulation time.
			B_.inh_queue_.AddSpike(e.get_rel_delivery_steps(nest::Time()), e.get_weight() * e.get_multiplicity());
		}
		else
		{
			std::ostringstream msg;
			msg << "Unexpected rport id: " << e.get_rport();
			throw nest::BadProperty(msg.str());
		}
	}

	/**
	* CurrentEvent handling.
	* @param e the event.
	*/
	void SrmPecevskiAlpha::handle(nest::CurrentEvent& e)
	{
		assert(e.get_delay() > 0);

		const double c = e.get_current();
		const double w = e.get_weight();

		B_.currents_.add_value(e.get_rel_delivery_steps(nest::kernel().simulation_manager.get_slice_origin()), w * c);
	}

	/**
	* DataLoggingRequest handling.
	* @param e the event.
	*/
	void SrmPecevskiAlpha::handle(nest::DataLoggingRequest& e)
	{
		B_.logger_.handle(e);
	}
}
