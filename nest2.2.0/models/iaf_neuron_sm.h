/*
 *  iaf_neuron_sm.h
 *
 *
 */

#ifndef IAF_NEURON_SM_H
#define IAF_NEURON_SM_H

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
Name: iaf_neuron_sm - Leaky integrate-and-fire neuron model.

Description:

  iaf_neuron_sm is an implementation of a leaky integrate-and-fire model
  with alpha-function shaped synaptic currents. 
  Inspired by [1].
  
  Solving system of differential equations:
  
  dV_m/dt = (E_L - V_m)/tau_m + I_e(t)/C_m + I_syn(t)/C_m,

 where I_syn(t) = sum(g_i(t)*(E_rev_i - V_m)),i - number of synapse connected with given neuron), g_i(t) = g_max_i*w_i(t)*sum(f_i(t-t_j),t_j<t),
 
 t_j - time of spike,
 H(t) - Heaviside function.

 f(t) = (t/tau_a)*exp(-t/tau_a)*H(t)

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
  tau_a		double - time constant in ms

  Following variables obtained from synapse:

  g_max		double - Maximal conductivity of synapse in nS
  w(t)		double - synaptic weight
  E_rev		double - synaptic reversal potentials in mV
  
 


Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest
      
*/

  /**
   * Leaky integrate-and-fire neuron.
   */
  class iaf_neuron_sm: public Archiving_Node
  {
    
  public:
    
    iaf_neuron_sm();
    iaf_neuron_sm(const iaf_neuron_sm&);
    ~iaf_neuron_sm();
    

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

    void set_learning_variables(DictionaryDatum &d, double t,std::vector<double > &currents);//prepare before learning
    void clear_spike_history();
    double get_sum_I_learn(double_t t);
    double get_sum_I(double_t t); // get summary current from all incoming synapse
    double heaviside (double t);
    double alpha_function(double t, double tau_a);
    double gauss(double t,double mu,double sigma);
    //bool learn_procedure(double_t time,double_t h,double_t v);
    bool emission_procedure(double_t potential);
    double normalization(double_t v);
	
  private:
	// The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<iaf_neuron_sm>;
    friend class UniversalDataLogger<iaf_neuron_sm>;

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
		
	  /** Threshold in mV.  */
      double_t V_th_;

      /** External current in pA */
      double_t I_e_;
      double_t tau_a_;//time constant in ms
      
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
      Buffers_(iaf_neuron_sm &);
      Buffers_(const Buffers_ &, iaf_neuron_sm &);

     /** buffers and summs up incoming spikes/currents */
      RingBuffer spikes_;
      RingBuffer currents_;

      //! Logger for all analog data
      UniversalDataLogger<iaf_neuron_sm> logger_;
    };
    

    // Access functions for UniversalDataLogger -------------------------------
    struct Variables_ { 
     int_t    RefractoryCounts_;  //!< refractory time in steps
     librandom::RngPtr rng_; // random number generator of my own thread
   };
    
    //! Read out the real membrane potential
    double_t get_V_m_() const { return S_.v_; }

   // ---------------------------------------------------------------- 
  
   Parameters_ P_;
   State_      S_;
   Variables_  V_;
   Buffers_    B_;

   //! Mapping of recordables names to access functions
   static RecordablesMap<iaf_neuron_sm> recordablesMap_;

   static int id_gl;
   //static int sh_gl;
   Iaf_neuron_sm_behavior* behav_loc; //nunost?
   Gradient_descent_behavior* behav_learning_loc;
  };


  
inline
port iaf_neuron_sm::check_connection(Connection& c, port receptor_type)
{
  SpikeEvent e;
  e.set_sender(*this);
  c.check_event(e);
  return c.get_target()->connect_sender(e, receptor_type);
}

inline
port iaf_neuron_sm::connect_sender(SpikeEvent&, port receptor_type)
{
  if (receptor_type != 0)
    throw UnknownReceptorType(receptor_type, get_name());
  return 0;
}
 
inline
port iaf_neuron_sm::connect_sender(CurrentEvent&, port receptor_type)
{
  if (receptor_type != 0)
    throw UnknownReceptorType(receptor_type, get_name());
  return 0;
}
 
inline
port iaf_neuron_sm::connect_sender(DataLoggingRequest &dlr, 
		      port receptor_type)
{
  if (receptor_type != 0)
    throw UnknownReceptorType(receptor_type, get_name());
  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

inline
void iaf_neuron_sm::get_status(DictionaryDatum &d) const
{
  P_.get(d);
  S_.get(d, P_);
  behav_learning_loc->get(d);//KAK
  Archiving_Node::get_status(d);
  (*d)[names::recordables] = recordablesMap_.get_list(); //danger
}

inline
void iaf_neuron_sm::set_status(const DictionaryDatum &d)
{
  Parameters_ ptmp = P_;  // temporary copy in case of errors
  ptmp.set(d);                       // throws if BadProperty
  State_      stmp = S_;  // temporary copy in case of errors
  stmp.set(d, ptmp);                 // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not 
  // write them back to (P_, S_) before we are also sure that 
  // the properties to be set in the parent class are internally 
  // consistent.
  behav_learning_loc->set(d); //KAK
  Archiving_Node::set_status(d);

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

inline
double iaf_neuron_sm::heaviside (double t)
{
	if (t>0)
	{
		return 1.0;
	}
	else
	{
		return 0.0;
	}
}

inline
double iaf_neuron_sm::alpha_function(double t, double tau_a)
{
	return (t/tau_a)*exp(-t/tau_a)*heaviside(t);
}

inline 
double iaf_neuron_sm::gauss(double t,double mu,double sigma)
{
  return exp(-((t-mu)*(t-mu))/(2*sigma*sigma))*heaviside(t);
}

inline
void iaf_neuron_sm::clear_spike_history()
{
  std::map<std::pair<index,index>,struct_iaf_neuron_sm >::iterator it=behav_loc->val.begin();
  while (it!=behav_loc->val.end())
  {
    (*it).second.base_val.v_time_spike_times.clear();
    it++;
  }
}

inline
double iaf_neuron_sm::get_sum_I(double_t t)
{
	/*заменить внизу:*/
	std::map<std::pair<index,index>,struct_iaf_neuron_sm >::iterator it=behav_loc->val.begin();
	double sum=0;
	while (it!=behav_loc->val.end())
	{
	    double sum_inter=0;
	    double t_spike;
	    //std::cout<<"g_max "<<(*it).second.g_max<<"\t";
	    //std::cout<<"E_rev "<<(*it).second.E_rev<<"\t";
	    //std::cout<<"weight "<<(*it).second.base_val.weight<<"\n";
	    for (int i=0;i!=(*it).second.base_val.v_time_spike_times.size();i++)
	    {
		t_spike=(*it).second.base_val.v_time_spike_times[i];
		//std::cout<<"t "<<t<<" t_spike "<<t_spike<<"\n";
		sum_inter+=alpha_function(t-t_spike,P_.tau_a_);
	    }
	    sum+=(*it).second.g_max*((*it).second.E_rev-S_.v_)*(*it).second.base_val.weight*sum_inter;
	    it++;
	}
	return sum;
}


inline
double iaf_neuron_sm::get_sum_I_learn(double_t t)
{
	//double learn_time=t-P_.bias_;
	double learn_time=t-behav_learning_loc->get_bias();
	int number_prim=behav_learning_loc->get_number_prim();
	//std::cout<<"learn time "<<learn_time<<"\n";
	std::map<std::pair<index,index>,struct_iaf_neuron_sm >::iterator it=behav_loc->val.begin();
	double sum=0;
	int sh=0;
	while (it!=behav_loc->val.end())
	{
	    double sum_inter=0;
	    double t_spike;
	    //std::cout<<"g_max "<<(*it).second.g_max<<"\t";
	    //std::cout<<"E_rev "<<(*it).second.E_rev<<"\t";
	    //std::cout<<"weight "<<(*it).second.base_val.weight<<"\n";
	    for (int i=0;i!=behav_learning_loc->v_learn_set_input[sh][number_prim].size();i++)
	    {
		t_spike=behav_learning_loc->v_learn_set_input[sh][number_prim][i];
		//std::cout<<"t "<<t<<" t_spike "<<t_spike<<"\n";
		if (t_spike>0.0)
		{
		  sum_inter+=alpha_function(learn_time-t_spike,P_.tau_a_);
		}
	    }
	    sum+=(*it).second.g_max*((*it).second.E_rev-S_.v_)*(*it).second.base_val.weight*sum_inter;
	    it++;
	    sh++;
	}
	return sum;
}



} // namespace


#endif /* #ifndef iaf_neuron_sm_H */
