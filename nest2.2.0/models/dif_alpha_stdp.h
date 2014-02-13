/*
 *  dif_alpha_stdp.h
 */

#ifndef DIF_ALPHA_STDP_H
#define DIF_ALPHA_STDP_H

#include "connection_het_wd.h"
#include "archiving_node.h"
#include "generic_connector.h"
#include "dictutils.h"
#include "nest_names.h"
#include <cmath>

/* BeginDocumentation
  Name: dif_alpha_stdp - Synapse type with simple plasticity.


   Parameters: 
     The following parameters can be set in the status dictionary:
	 tau_a : double time for post synaptic current

  Transmits: SpikeEvent
       
*/



namespace nest {

class Dif_alpha_stdp : public ConnectionHetWD
{
 public:

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  Dif_alpha_stdp();
  
  Dif_alpha_stdp(const Dif_alpha_stdp &);

  /**
   * Default Destructor.
   */
  ~Dif_alpha_stdp() {/*std::cout<<"sm destr "<<"\n";*/}
  
  void check_connection(Node & s, Node & r, rport receptor_type, double_t t_lastspike);

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
   * \param cp Common properties to all synapses (empty).
   */
  void send(Event& e, double_t t_lastspike, const CommonSynapseProperties &cp);

  // overloaded for all supported event types
  using Connection::check_event;
  void check_event(SpikeEvent&) {}
  
 private:
   
  double_t facilitate_(double_t w, double_t kplus);
  double_t depress_(double_t w, double_t kminus);
  
  //tsodyks param
  double_t U_;
  double_t D_;
  double_t F_;
  double_t u_;
  double_t R;
  //trosdyks param
  double_t tau_plus_;
  double_t tau_minus_;
  double_t lambda_;
  double_t alpha_;
  double_t mu_plus_;
  double_t mu_minus_;  
  double_t Wmax_;
   
  double_t tau_a_;		 
  int t_count;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */

inline
double_t Dif_alpha_stdp::facilitate_(double_t w, double_t kplus)
{
  double_t norm_w = (w / Wmax_) + (lambda_ * std::pow(1.0 - (w/Wmax_), mu_plus_) * kplus);
  return norm_w < 1.0 ? norm_w * Wmax_ : Wmax_;
}

inline 
double_t Dif_alpha_stdp::depress_(double_t w, double_t kminus)
{
  //double_t norm_w = (w / Wmax_) - (alpha_ * lambda_ * std::pow(w/Wmax_, mu_minus_) * kminus);
  //return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
  double t_w_new=w-alpha_ * lambda_*kminus;
  if (0.0>=t_w_new)
  {
    return 0.0;
  }
  else
  {
    return t_w_new;
  }
}

inline 
void Dif_alpha_stdp::check_connection(Node & s, Node & r, rport receptor_type, double_t t_lastspike)
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

inline
void Dif_alpha_stdp::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{
  if (t_count==0)
  {
    //target_->behav->set_variables(id_target,id_source,tau_a_,weight_,weight_); //modif
    DictionaryDatum *d = new DictionaryDatum(new Dictionary);
    def<double>((*d), "weight", weight_);
    def<double>((*d), "tau_a", tau_a_);
    def<double>((*d), "lambda", lambda_);
    def<double>((*d), "tau_plus", tau_plus_);
    def<double>((*d), "Wmax", Wmax_);
    target_->behav->set(id_target,id_source,(*d));
    t_count++;
  }
  double_t t_weight;
  DictionaryDatum *d3 = new DictionaryDatum(new Dictionary);
  target_->behav->get_parametrs(id_target,id_source,t_weight,(*d3));
  delete d3;
  weight_ =t_weight;
  
  double_t t_spike = e.get_stamp().get_ms();
  // t_lastspike_ = 0 initially
  double_t dendritic_delay = Time(Time::step(delay_)).get_ms();
  
  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;
  
  target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,
                         &start, &finish);
  double_t minus_dt;
  std::deque<histentry>::iterator post_first_spike=start;
  std::deque<histentry>::iterator post_last_spike=finish;
  while (start != finish)
  {
    post_last_spike=start;
    minus_dt = t_lastspike - (start->t_ + dendritic_delay);
    start++;
    if (minus_dt == 0)
      continue;
    //weight_ = facilitate_(weight_,  std::exp(minus_dt / tau_plus_));
  }
  if (post_last_spike!=finish)
  {
    minus_dt = -t_spike + (post_last_spike->t_ + dendritic_delay);
    //if ((minus_dt != 0)&&(std::abs(minus_dt)>2.0)&&(std::abs(minus_dt)<3*20.0))
    if ((minus_dt != 0)&&(std::abs(minus_dt)>2.0))
      //weight_ = depress_(weight_, std::exp((minus_dt) / 20.0));
      weight_ = depress_(weight_, std::exp(minus_dt/ tau_minus_));
      //weight_ = depress_(weight_, std::exp((minus_dt) / 20.0));
  }
  // synapse STDP depressing/facilitation dynamics
  //tsodyk
  
  double_t u_old=u_;
  double_t R_old=R;
  u_=U_+u_old*(1.0-U_)*std::exp(-(t_spike-t_lastspike)/F_);
  R=1.0+(R_old-u_old*R_old-1.0)*std::exp(-(t_spike-t_lastspike)/D_);
 
  //tsodyk
 
  //std::cout<<"target_adress behav "<<target_->behav<<"\n";
  //std::cout<<"target_adress behav_loc "<<target_->behav_loc<<"\n";
  //std::cout<<"weight in synapse "<<e.get_weight()<<"\n";
 
  
  
  //depression due to new pre-synaptic spike
  double_t dynamic_weight=weight_*u_*R;
  //std::cout<<"id_source="<<id_source<<" dynamic_weight="<<dynamic_weight<<"\n";
  target_->behav->add_spike_time(id_target,id_source,0,e.get_stamp().get_ms());
  DictionaryDatum *d2 = new DictionaryDatum(new Dictionary);
  def<double>((*d2), "weight",weight_);
  def<double>((*d2), "tau_a",tau_a_); 
  def<double>((*d2),"dynamic_weight",dynamic_weight);
  target_->behav->set_variables(id_target,id_source,(*d2));
  delete d2;
  e.set_receiver(*target_);
  e.set_weight(weight_*u_*R);
  //e.set_weight(weight_ );
  e.set_delay( delay_ );
  e.set_rport( rport_ );
  e();

}
 
} // namespace

#endif
