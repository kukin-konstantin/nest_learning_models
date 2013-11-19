/*
 *  iaf_neuron_const.h
 *
 *
 */

#ifndef IAF_NEURON_CONST_H
#define IAF_NEURON_CONST_H

#include <map>
#include "nest.h"
#include "event.h"
#include "mt19937.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"

namespace nest
{
  class Network;

  /* BeginDocumentation
Name: iaf_neuron_const - Leaky integrate-and-fire neuron model.

Description:

  iaf_neuron_const is an implementation of a leaky integrate-and-fire model
  with alpha-function shaped synaptic currents. 
  Inspired by [1].
  
  Solving system of differential equations:
  
  dV_m/dt = (E_L - V_m)/tau_m + I_e(t)/C_m + I_syn(t)/C_m,

 where I_syn(t) = sum(Wi),i - number of synapse connected with given neuron)
 

 if V_m(t)>=V_th => V_m=E_L and neuron fire spike

 

  


Parameters: 

  The following parameters can be set in the status dictionary.

  V_m        double - Membrane potential in mV 
  E_L        double - Resting membrane potential in mV. 
  C_m        double - Capacity of the membrane in pF
  tau_m      double - Membrane time constant in ms.
  t_ref      double - Duration of refractory period in ms. 
  V_th       double - Spike threshold in mV.
  V_reset    double - Reset potential of the membrane in mV.
  I_e        double - Constant external input current in pA.

  
 


Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest
      
*/

  /**
   * Leaky integrate-and-fire neuron.
   */
  class iaf_neuron_const: public Archiving_Node
  {
    
  public:
    
    iaf_neuron_const();
    iaf_neuron_const(const iaf_neuron_const&);
    //~iaf_neuron_const();
    

    /**
     * Import sets of overloaded virtual functions.
     * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
     */

    using Node::connect_sender;
    using Node::handle;

    port check_connection(Connection&, port);
    
    void handle(SpikeEvent &);
    void handle(CurrentEvent &);
    void handle(DataLoggingRequest &);
    
    port connect_sender(SpikeEvent&, port);
    port connect_sender(CurrentEvent&, port);
    port connect_sender(DataLoggingRequest&, port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

	
  private:
	// The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<iaf_neuron_const>;
    friend class UniversalDataLogger<iaf_neuron_const>;

    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();

    void update(Time const &, const long_t, const long_t); 

    
     
    // ---------------------------------------------------------------- 

    /** 
     * Independent parameters of the model. 
     */
    struct Parameters_ {
      /** Membrane capacitance in pF. */
      double_t C_;
    
      /** Membrane time constant in ms. */
      double_t Tau_; 
      
      /** Refractory period in ms. */
      double_t TauR_;

      /** Resting potential in mV. */
      double_t E_L_;
      
      /** Reset_value */
      double_t V_reset_;
		
	  /** Threshold in mV.  */
      double_t V_th_;

      /** External current in pA */
      double_t I_e_;
      
      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&);  //!< Set values from dicitonary
    };

    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     */
    struct State_ {
      double_t     v_; // membrane potential
      double_t     I_; // input current
      
      int_t    r_;  //!< number of refractory steps remaining

      State_();  //!< Default initialization
      
      void get(DictionaryDatum&, const Parameters_&) const;
      void set(const DictionaryDatum&, const Parameters_&);
    };
    
    // ---------------------------------------------------------------- 

    /**
     * Buffers of the model.
     */
     struct Buffers_ {
      Buffers_(iaf_neuron_const &);
      Buffers_(const Buffers_ &, iaf_neuron_const &);

     /** buffers and summs up incoming spikes/currents */
      RingBuffer spikes_;
      RingBuffer currents_;

      //! Logger for all analog data
      UniversalDataLogger<iaf_neuron_const> logger_;
    };
    

    // Access functions for UniversalDataLogger -------------------------------
    struct Variables_ { 
     int_t    RefractoryCounts_;  //!< refractory time in steps
   };
    
    //! Read out the real membrane potential
    double_t get_V_m_() const { return S_.v_; }

   // ---------------------------------------------------------------- 
  
   Parameters_ P_;
   State_      S_;
   Variables_  V_;
   Buffers_    B_;

   //! Mapping of recordables names to access functions
   static RecordablesMap<iaf_neuron_const> recordablesMap_;

  };


  
inline
port iaf_neuron_const::check_connection(Connection& c, port receptor_type)
{
  SpikeEvent e;
  e.set_sender(*this);
  c.check_event(e);
  return c.get_target()->connect_sender(e, receptor_type);
}

inline
port iaf_neuron_const::connect_sender(SpikeEvent&, port receptor_type)
{
  if (receptor_type != 0)
    throw UnknownReceptorType(receptor_type, get_name());
  return 0;
}
 
inline
port iaf_neuron_const::connect_sender(CurrentEvent&, port receptor_type)
{
  if (receptor_type != 0)
    throw UnknownReceptorType(receptor_type, get_name());
  return 0;
}
 
inline
port iaf_neuron_const::connect_sender(DataLoggingRequest &dlr, 
		      port receptor_type)
{
  if (receptor_type != 0)
    throw UnknownReceptorType(receptor_type, get_name());
  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

inline
void iaf_neuron_const::get_status(DictionaryDatum &d) const
{
  P_.get(d);
  S_.get(d, P_);
  Archiving_Node::get_status(d);
  (*d)[names::recordables] = recordablesMap_.get_list(); //danger
}

inline
void iaf_neuron_const::set_status(const DictionaryDatum &d)
{
  Parameters_ ptmp = P_;  // temporary copy in case of errors
  ptmp.set(d);                       // throws if BadProperty
  State_      stmp = S_;  // temporary copy in case of errors
  stmp.set(d, ptmp);                 // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not 
  // write them back to (P_, S_) before we are also sure that 
  // the properties to be set in the parent class are internally 
  // consistent.
  Archiving_Node::set_status(d);

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}




} // namespace


#endif /* #ifndef iaf_neuron_const_H */
