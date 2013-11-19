/*
 *  sm_connection.h
 */

#ifndef SM_CONNECTION_H
#define SM_CONNECTION_H

#include "connection_het_wd.h"

/* BeginDocumentation
  Name: sm_synapse - Synapse type with simple plasticity.


   Parameters: 
     The following parameters can be set in the status dictionary:
	 g        double - Conductance in nS;
	 E		  double - Reversal potential;

  Transmits: SpikeEvent
       
*/



namespace nest {

class SMConnection : public ConnectionHetWD
{
 public:

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  SMConnection();
  
  SMConnection(const SMConnection &);

  /**
   * Default Destructor.
   */
  ~SMConnection() {/*std::cout<<"sm destr "<<"\n";*/}

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
  double_t g_;       //!< Leak Conductance in nS
  double_t E_;		 //!< Reversal potential in mV
  int t_count;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void SMConnection::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{  
  //std::cout<<"target_adress behav "<<target_->behav<<"\n";
  //std::cout<<"target_adress behav_loc "<<target_->behav_loc<<"\n";
  //std::cout<<"weight in synapse "<<e.get_weight()<<"\n";
  if (t_count==0)
  {
    target_->behav->set_variables(id_target,id_source,g_,weight_,E_); //modif
    //std::cout<<"YYYYYYYYYYY"<<"\n";
    t_count++;
  }
  
  target_->behav->add_spike_time(id_target,id_source,e.get_stamp().get_ms());
  
  e.set_receiver(*target_);
  e.set_weight( g_*weight_ );
  e.set_delay( delay_ );
  e.set_rport( rport_ );
  std::cout<<"spike send "<<e.get_stamp().get_ms()<<"\n";
  e();
}
 
} // namespace

#endif
