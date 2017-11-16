/*
*  srm_pecevski_alpha.h
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

#ifndef SRM_PECEVSKI_ALPHA_H
#define SRM_PECEVSKI_ALPHA_H

#include "event.h"
#include "nest.h"
#include "clipped_randomdev.h"
#include "gamma_randomdev.h"
#include "normal_randomdev.h"
#include "poisson_randomdev.h"
#include "ring_buffer.h"
#include "spike_queue.h"
#include "universal_data_logger.h"

#include "archiving_node.h"

namespace sam
{
	/**
	* @brief Point process neuron with alpha-shaped PSPs based on Pecevski et al 2016[1].
	*
	* SrmPecevskiAlpha is a variant of the spike response model model with
	* alpha-shaped PSP modified from vanilla alpha so that only the top half of the
	* alpha kernel is used as a response.
	*
	* \f$\epsilon(t) \;=\; \epsilon_o \left( e^1 (\frac{t}{\tau_\alpha} + t_1)
	* e^{-(\frac{t}{\tau_\alpha} + t_1)} - \frac{1}{2}\right)\f$
	*
	* Spikes are generated randomly according to the current value of the
	* transfer function which operates on the membrane potential. Spike
	* generation is followed by an optional dead time (refractory state).
	* Setting \a with_reset to true will reset the membrane potential after
	* each spike.
	*
	* The transfer function can be chosen to be linear, exponential or a sum of
	* both by adjusting three parameters:
	*
	*     \f[
	*        rate = Rect[ c1 * V' + c2 * exp(c3 * V') ],
	*     \f]
	* where the effective potential \f$ V' = V_m + E_{sfa}\f$ and \f$ E_{sfa} \f$
	* is the the intrinsic bias. By setting \a c3 = 0, \a c2 can be used as an
	* offset spike rate for an otherwise linear rate model.
	*
	* The dead time enables refractoriness. If dead time is 0, the
	* number of spikes in one time step might exceed one and is drawn from the
	* Poisson distribution accordingly. Otherwise, the probability for a spike
	* is given by \f$ 1 - exp(-rate*h) \f$, where h is the simulation time step.
	* If dead_time is smaller than the simulation resolution (time step), it is
	* internally set to the time step. Note that, even if non-refractory neurons
	* are to be modeled, a small value of dead_time, like \a dead_time=1e-8, might
	* be the value of choice since it uses faster uniform random numbers than
	* \a dead_time=0, which draws Poisson numbers. Only for very large spike rates
	* (> 1 spike/h) this will cause errors.
	*
	* The model implements intrinsic plasticity. If the neuron spikes, the
	* bias increases and the membrane potential will take less time to reach it.
	* If the neuron does not spike, the threshold linearly decays over time,
	* decrease the firing probability of the neuron.
	*
	* Bias can be generated randomly from a Gaussian distribution clipped by [min_bias,
	* max_bias], mean mu_bias, standard deviation sigma_bias. In that case initial bias
	* parameter is ignored.
	*
	* This model has been adapted from poisson_dbl_exp_neuron. The default parameters
	* are set to the mean values in [2], which have were matched to spike-train
	* recordings.
	*
	* <b>Parameters</b>
	*
	* The following parameters can be set in the status dictionary.
	* (default values and constraints are given in parentheses):
	*
	* <table>
	* <tr><th>name</th>                       <th>type</th>   <th>comment</th></tr>
	* <tr><td>\a V_m</td>                     <td>double</td> <td>Membrane potential [mV]</td></tr>
	* <tr><td>\a rect_exc</td>                <td>bool</td>   <td>Use a rectangular EPSP (false)</td></tr>
	* <tr><td>\a rect_inh</td>                <td>bool</td>   <td>Use a rectangular IPSP (false)</td></tr>
	* <tr><td>\a r_m</td>                     <td>double</td> <td>Membrane resistance (1.0) [GOhm] </td></tr>
	* <tr><td>\a e_0_exc</td>				  <td>double</td> <td>Amplitude factor of EPSPs [mV]</td></tr>
	* <tr><td>\a e_0_inh</td>				  <td>double</td> <td>Amplitude factor of IPSPs [mV]</td></tr>
	* <tr><td>\a tau_exc</td>				  <td>double</td> <td>Alpha EPSP time constant [ms]</td></tr>
	* <tr><td>\a tau_inh</td>				  <td>double</td> <td>Alpha IPSP time constant [ms]</td></tr>
    * <tr><td>\a tau_m</td>		              <td>double</td> <td>Impulse current response time constant [ms]</td></tr>
	* <tr><td>\a dead_time</td>               <td>double</td> <td>Duration of the dead time (1.0, &ge 0.0) [ms]</td></tr>
	* <tr><td>\a dead_time_random</td>        <td>bool</td>   <td>Should a random dead time be drawn after each
	*                                                             spike? (false) </td></tr>
	* <tr><td>\a dead_time_shape</td>         <td>int</td>    <td>Shape parameter of dead time gamma distribution
	*                                                             (1, &ge 1) </td></tr>
	* <tr><td>\a t_ref_remaining</td>         <td>double</td> <td>Remaining dead time at simulation start (>0.0) [ms]
	*                                                             (0.0, &ge 0.0) </td></tr>
	* <tr><td>\a with_reset</td>              <td>bool</td>   <td>Should the membrane potential be reset after a
	*                                                             spike?</td></tr>
	* <tr><td>\a I_e</td>                     <td>double</td> <td>Constant input current (0.0) [pA]</td></tr>
	* <tr><td>\a c_1</td>                     <td>double</td> <td>Slope of linear part of transfer function
	*                                                             (0.0) [Hz/mV]</td></tr>
	* <tr><td>\a c_2</td>                     <td>double</td> <td>Prefactor of exponential part of transfer function
	*                                                             (1.238) [Hz/mV]</td></tr>
	* <tr><td>\a c_3</td>                     <td>double</td> <td>Coefficient of exponential non-linearity of transfer
	*                                                             function (0.25, &ge 0.0) [1/mV]</td></tr>
	* <tr><td>\a b_baseline</td>              <td>double</td> <td>Intrinsic excitability baseline (-1.0)</td></tr>
	* <tr><td>\a eta_bias</td>                <td>double</td> <td>Intrinsic excitability learning rate (0.1)</td></tr>
	* <tr><td>\a tau_bias</td>                <td>double</td> <td>Coefficient of bias updates (15.0) [ms]</td></tr>
	* <tr><td>\a max_bias</td>                <td>double</td> <td>Maximum bias value (5.0)</td></tr>
	* <tr><td>\a min_bias</td>                <td>double</td> <td>Minimum bias value (-30.0)</td></tr>
	* <tr><td>\a T</td>                       <td>double</td> <td>Bias update scaling parameter (0.58)</td></tr>
	* <tr><td>\a bias</td>                    <td>double</td> <td>Initial bias (5.0)</td></tr>
	* </table>
	*
	* <i>Sends:</i> SpikeEvent
	*
	* <i>Receives:</i> SpikeEvent, CurrentEvent, DataLoggingRequest
	*
	* <b>References</b>
	*
	* [1] Learning Probabilistic Inference through Spike-Timing-Dependent
	* Plasticity (2016) D. Peceveski and W. Maass, eNeuro
	*
	* [2] Predicting spike timing of neocortical pyramidal neurons by simple
	* threshold models (2006) Jolivet R, Rauch A, Luescher H-R, Gerstner W.
	* J Comput Neurosci 21:35-49.
	*
	* [3] Pozzorini C, Naud R, Mensi S, Gerstner W (2013) Temporal whitening by
	* power-law adaptation in neocortical neurons. Nat Neurosci 16: 942-948.
	* (uses a similar model of multi-timescale adaptation).
	*
	* [4] Grytskyy D, Tetzlaff T, Diesmann M and Helias M (2013) A unified view
	* on weakly correlated recurrent networks. Front. Comput. Neurosci. 7:131.
	*
	* [5] Deger M, Schwalger T, Naud R, Gerstner W (2014) Fluctuations and
	* information filtering in coupled populations of spiking neurons with
	* adaptation. Physical Review E 90:6, 062704.
	*
	* @author  D'Amato; (of poisson_dbl_exp_neuron) Kappel, Hsieh; (of pp_psc_delta) July 2009, Deger, Helias; January 2011, Zaytsev; May 2014, Setareh
	* @see TracingNode
	*/
	class SrmPecevskiAlpha : public nest::Archiving_Node
	{
	public:

		SrmPecevskiAlpha();
		SrmPecevskiAlpha(const SrmPecevskiAlpha&);

		/**
		* Import sets of overloaded virtual functions.
		* @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
		*/
		using nest::Node::handle;
		using nest::Node::handles_test_event;

		nest::port send_test_event(nest::Node&, nest::rport, nest::synindex, bool);

		void handle(nest::SpikeEvent &);
		void handle(nest::CurrentEvent &);
		void handle(nest::DataLoggingRequest &);

		nest::port handles_test_event(nest::SpikeEvent&, nest::rport);
		nest::port handles_test_event(nest::CurrentEvent&, nest::rport);
		nest::port handles_test_event(nest::DataLoggingRequest&, nest::rport);

		void get_status(DictionaryDatum &) const;
		void set_status(const DictionaryDatum &);

	private:

		void init_state_(const nest::Node& proto);
		void init_buffers_();
		void calibrate();
		double kernel(const double time_since_spike, const bool use_exc_kernel) const;
		double get_psp_sum(const nest::Time& now, const bool use_exc_psp);

		void update(nest::Time const &, const long, const long);

		// The next two classes need to be friends to access the State_ class/member
		friend class nest::RecordablesMap<SrmPecevskiAlpha>;
		friend class nest::UniversalDataLogger<SrmPecevskiAlpha>;

		// ----------------------------------------------------------------

		/**
		* Independent parameters of the model.
		*/
		struct Parameters_
		{
            /** If true, uses a rectangular PSP rather than the half-alpha shape. */
            bool use_rect_psp_exc_;

            /** If true, uses a rectangular PSP rather than the half-alpha shape. */
            bool use_rect_psp_inh_;

            /** Prefactor that converts input current into voltages. */
            double resistance_;

			/** Amplitude of excitatory alpha PSP before truncation (if not using rect PSPs),
			 * otherwise it is twice rect amplitude. */
			double epsilon_0_exc_;

            /** Amplitude of inhibitory alpha PSP before truncation (if not using rect PSPs),
             * otherwise it is twice rect amplitude. */
			double epsilon_0_inh_;

			/*** Excitatory alpha PSP time constant, or rect PSP duration. */
			double tau_alpha_exc_;

			/*** Inhibitory alpha PSP time constant, or rect PSP duration. */
			double tau_alpha_inh_;

            /** Membrane time constant (for current effects on voltage). */
            double tau_membrane_;

			/** Dead time in ms. */
			double dead_time_;

			/** Do we use random dead time? */
			bool dead_time_random_;

			/** Shape parameter of random dead time gamma distribution. */
			long dead_time_shape_;

			/** Do we reset the membrane potential after each spike? */
			bool with_reset_;

			/** Slope of the linear part of transfer function. */
			double c_1_;

			/** Prefactor of exponential part of transfer function. */
			double c_2_;

			/** Coefficient of exponential non-linearity of transfer function. */
			double c_3_;

			/** External DC current. */
			double I_e_;

			/** Dead time from simulation start. */
			double t_ref_remaining_;

            /** Intrinsic plasticity baseline parameter, see Peceveski et al. 2016. */
            double bias_baseline_;

            /** Learning rate of intrinsic plasticity. */
            double eta_bias_;

            /** Coefficient of bias updates on spiking. */
            double tau_bias_;

            /** Bias value range. */
            double max_bias_;
            double min_bias_;

            /** Bias update scaling factor, for use in learning rule. */
            double t_;

			Parameters_(); //!< Sets default parameter values

			void get(DictionaryDatum&) const; //!< Store current values in dictionary
			void set(const DictionaryDatum&); //!< Set values from dictionary
		};

		/**
		* State variables of the model.
		*/
		struct State_
		{
			double u_membrane_; //!< The membrane potential
			double input_current_; //!< The piecewise linear input currents
			double bias_; //!< Intrinsic plasticity bias
            double u_i_; //!< Current-induced voltage
			int r_; //!< Number of refractory steps remaining

			State_(); //!< Default initialization

			void get(DictionaryDatum&, const Parameters_&) const;
			void set(const DictionaryDatum&, const Parameters_&);
		};

		/**
		* Buffers of the model.
		*/
		struct Buffers_
		{
			Buffers_(SrmPecevskiAlpha &);
			Buffers_(const Buffers_ &, SrmPecevskiAlpha &);

			/** Buffers and sums up incoming currents */
			nest::RingBuffer currents_;

			/** Queues time-amplitude pairs for SRM kernel calculations */
			SpikeQueue exc_queue_;
			SpikeQueue inh_queue_;

			//! Logger for all analog data
			nest::UniversalDataLogger<SrmPecevskiAlpha> logger_;
		};

		/**
		* Internal variables of the model.
		*/
		struct Variables_
		{
			double t_1; //!< First time at half-amplitude of generic alpha kernel
			double h_; //!< simulation time step in ms
			double dt_rate_; //!< rate parameter of dead time distribution
            double propagator_;

			librandom::RngPtr rng_; // random number generator of my own thread
			librandom::PoissonRandomDev poisson_dev_; // random deviate generator
			librandom::GammaRandomDev gamma_dev_; // random deviate generator

			int DeadTimeCounts_;

		};

		// Access functions for UniversalDataLogger

		//! Read out the real membrane potential

		double get_V_m_() const
		{
			return S_.u_membrane_;
		}

		//! Read out the adaptive threshold potential

		double get_E_sfa_() const
		{
			return S_.bias_;
		}

		/**
		* Instances of private data structures for the different types
		* of data pertaining to the model.
		* @note The order of definitions is important for speed.
		* @{
		*/
		Parameters_ P_;
		State_ S_;
		Variables_ V_;
		Buffers_ B_;
		/** @} */

		//! Mapping of recordables names to access functions
		static nest::RecordablesMap<SrmPecevskiAlpha> recordablesMap_;
	};

	/**
	* SrmPecevskiAlpha test event.
	*/
	inline
		nest::port SrmPecevskiAlpha::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
	{
		nest::SpikeEvent e;
		e.set_sender(*this);

		return target.handles_test_event(e, receptor_type);
	}

	/**
	* PoissonDblExpNeuron test event.
	*/
	inline
		nest::port SrmPecevskiAlpha::handles_test_event(nest::SpikeEvent&, nest::rport receptor_type)
	{
		if ((receptor_type != 0) && (receptor_type != 1))
		{
			throw nest::UnknownReceptorType(receptor_type, get_name());
		}

		return receptor_type;
	}

	/**
	* SrmPecevskiAlpha test event.
	*/
	inline
		nest::port SrmPecevskiAlpha::handles_test_event(nest::CurrentEvent&, nest::rport receptor_type)
	{
		if (receptor_type != 0)
		{
			throw nest::UnknownReceptorType(receptor_type, get_name());
		}

		return 0;
	}

	/**
	* SrmPecevskiAlpha test event.
	*/
	inline
		nest::port SrmPecevskiAlpha::handles_test_event(nest::DataLoggingRequest &dlr,
			nest::rport receptor_type)
	{
		if (receptor_type != 0)
		{
			throw nest::UnknownReceptorType(receptor_type, get_name());
		}

		return B_.logger_.connect_logging_device(dlr, recordablesMap_);
	}

	/**
	* Status getter function.
	*/
	inline
		void SrmPecevskiAlpha::get_status(DictionaryDatum &d) const
	{
		P_.get(d);
		S_.get(d, P_);

		(*d)[nest::names::recordables] = recordablesMap_.get_list();
	}

	/**
	* Status setter function.
	*/
	inline
		void SrmPecevskiAlpha::set_status(const DictionaryDatum &d)
	{
		Parameters_ ptmp = P_; // temporary copy in case of errors
		ptmp.set(d); // throws if BadProperty
		State_ stmp = S_; // temporary copy in case of errors
		stmp.set(d, ptmp); // throws if BadProperty

						   // We now know that (ptmp, stmp) are consistent. We do not
						   // write them back to (P_, S_) before we are also sure that
						   // the properties to be set in the parent class are internally
						   // consistent.
						   //Archiving_Node::set_status(d);

						   // if we get here, temporaries contain consistent set of properties
		P_ = ptmp;
		S_ = stmp;
	}

}

#endif