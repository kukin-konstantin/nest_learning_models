/*
 *  stdp_tsodyk.h
 *
 *
 */

#ifndef STDP_TSODYK_CONNECTION_H
#define STDP_TSODYK_CONNECTION_H

#include "connection_het_wd.h"
#include "archiving_node.h"
#include "generic_connector.h"
#include <cmath>

namespace nest
{
  class STDP_Tsodyk_Connection : public ConnectionHetWD
  {

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDP_Tsodyk_Connection();

  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDP_Tsodyk_Connection(const STDP_Tsodyk_Connection &);

  /**
   * Default Destructor.
   */
  ~STDP_Tsodyk_Connection() {}

  void check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike); // проверить еще

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status(DictionaryDatum & d) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status(const DictionaryDatum & d, ConnectorModel &cm);

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void set_status(const DictionaryDatum & d, index p, ConnectorModel &cm);

  /**
   * Create new empty arrays for the properties of this connection in the given
   * dictionary. It is assumed that they are not existing before.
   */
  void initialize_property_arrays(DictionaryDatum & d) const;

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void append_properties(DictionaryDatum & d) const;

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses (empty).
   */
  void send(Event& e, double_t t_lastspike, const CommonSynapseProperties &cp);

  // overloaded for all supported event types
  using Connection::check_event;
  void check_event(SpikeEvent&) {}

 private:

  double_t facilitate_(double_t w, double_t kplus);
  double_t depress_(double_t w, double_t kminus);

  // data members of each connection
  double_t tau_psc_;   //!< [ms] time constant of postsyn current
  double_t tau_fac_;   //!< [ms] time constant for fascilitation
  double_t tau_rec_;   //!< [ms] time constant for recovery
  double_t U_;         //!< asymptotic value of probability of release
  double_t x_;         //!< amount of resources in recovered state
  double_t y_;         //!< amount of resources in active state
  double_t u_;         //!< actual probability of release
  double_t tau_plus_;
  double_t lambda_;
  double_t alpha_;
  double_t mu_plus_;
  double_t mu_minus_;  
  double_t Wmax_;
  double_t Kplus_;

  };


inline
double_t STDP_Tsodyk_Connection::facilitate_(double_t w, double_t kplus)
{
  double_t norm_w = (w / Wmax_) + (lambda_ * std::pow(1.0 - (w/Wmax_), mu_plus_) * kplus);
  //std::cout<<"fac norm_w="<<norm_w<<"\t"<<"w="<<w<<"\tWmax="<<Wmax_<<"\t"<<"new w="<<norm_w * Wmax_<<"\n";
  return norm_w < 1.0 ? norm_w * Wmax_ : Wmax_;
}

inline 
double_t STDP_Tsodyk_Connection::depress_(double_t w, double_t kminus)
{
  double_t norm_w = (w / Wmax_) - (alpha_ * lambda_ * std::pow(w/Wmax_, mu_minus_) * kminus);
  //std::cout<<"depr norm_w="<<norm_w<<"\t"<<"w="<<w<<"\tWmax="<<Wmax_<<"\t"<<"new w="<<norm_w * Wmax_<<"\n";;
  return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
}


inline 
void STDP_Tsodyk_Connection::check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike)
{
  ConnectionHetWD::check_connection(s, r, receptor_type, t_lastspike);

  // For a new synapse, t_lastspike contains the point in time of the last spike.
  // So we initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of history[0, ..., t_last_spike - dendritic_delay] will be 
  // incremented by the following call to Archiving_Node::register_stdp_connection().
  // See bug #218 for details.
  r.register_stdp_connection(t_lastspike - Time(Time::step(delay_)).get_ms());
}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void STDP_Tsodyk_Connection::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{
  // synapse STDP depressing/facilitation dynamics
  double_t t_tmp=e.get_stamp().get_ms();
  double_t h =t_tmp - t_lastspike;
  double_t t_spike = t_tmp;
  // t_lastspike_ = 0 initially
  // this has no influence on the dynamics, IF y = z = 0 initially
  // !!! x != 1.0 -> z != 0.0 -> t_lastspike_=0 has influence on dynamics
  double_t dendritic_delay = Time(Time::step(delay_)).get_ms();
   
  // propagator
  // TODO: use expm1 here instead, where applicable
  double_t Puu = (tau_fac_ == 0.0) ? 0.0 : std::exp(-h/tau_fac_);
  double_t Pyy = std::exp(-h/tau_psc_);
  double_t Pzz = std::exp(-h/tau_rec_);
    
  //double_t Pzy = (Pyy - Pzz) * tau_rec_ / (tau_psc_ - tau_rec_);
  double_t Pxy = ((Pzz - 1.0)*tau_rec_ - (Pyy - 1.0)*tau_psc_) / (tau_psc_ - tau_rec_);
  double_t Pxz = 1.0 - Pzz;

  double_t z = 1.0 - x_ - y_;

  // propagation t_lastspike -> t_spike
  // don't change the order !

  u_ *= Puu;
  x_ += Pxy * y_ + Pxz * z;
  //z = Pzz * z_ + Pzy * y_;
  y_ *= Pyy;

  // delta function u
  u_ += U_*(1.0-u_);

  // postsynaptic current step caused by incoming spike
  double_t delta_y_tsp = u_*x_;

  // delta function x, y
  x_ -= delta_y_tsp;
  y_ += delta_y_tsp;

  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;

  // For a new synapse, t_lastspike contains the point in time of the last spike.
  // So we initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of history[0, ..., t_last_spike - dendritic_delay] have been 
  // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for details.
  target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,
                         &start, &finish);
  //facilitation due to post-synaptic spikes since last pre-synaptic spike
  double_t minus_dt;
  //std::cout<<"weight before "<<weight_<<"\n";
  while (start != finish)
  {
    minus_dt = t_lastspike - (start->t_ + dendritic_delay);
    start++;
    if (minus_dt == 0)
      continue;
    weight_ = facilitate_(weight_, Kplus_ * std::exp(minus_dt / tau_plus_));
  }

  //depression due to new pre-synaptic spike
//   std::cout<<"weight middle "<<weight_<<"\n";
  weight_ = depress_(weight_, target_->get_K_value(t_spike - dendritic_delay));
  
//   std::cout<<"weight after "<<weight_<<"\n";
  e.set_receiver(*target_);
  e.set_weight(weight_ * delta_y_tsp);
  e.set_delay(delay_);
  e.set_rport(rport_);
  e();
 
  Kplus_ = Kplus_ * std::exp((t_lastspike - t_spike) / tau_plus_) + 1.0;
}

} // of namespace nest

#endif // of #ifndef STDP_TSODYK_CONNECTION_H
