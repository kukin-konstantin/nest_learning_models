/*
 *  dif_alpha_connection.cpp
 *
 */

#include "dif_alpha_connection.h"
#include "network.h"
#include "connector_model.h"

namespace nest
{

  Dif_alpha_connection::Dif_alpha_connection() :
    ConnectionHetWD(),
    tau_a_(10.0),
    t_count(0.0),
    local_id_synapse(0.0)
  {
	  //std::cout<<"sm born "<<"\n";
 	  //target_->behav;
	  //target_->behav->set_variables(id_target,id_source,g_,weight_,E_);
	  //target_->behav->set_variables(id_target,id_source,0,0,0);
	  //target_->set_reversal_potential(значение потенциала,ссылка на синапс(или ид синапса))
	  //target_->set_conductance(значение проводимость,ссылка на синапс(или ид синапса))
	  //сделать у таргета ссылку на connection
	  //здесь надо один раз назначать неменяемые параметры на таргете -
	  //реверсальный потенциал, проводимость если она неизменна  
  }

  Dif_alpha_connection::Dif_alpha_connection(const Dif_alpha_connection &rhs):
  ConnectionHetWD(rhs)
  {
    tau_a_=rhs.tau_a_;
    t_count=rhs.t_count;
    local_id_synapse=rhs.local_id_synapse;
    //std::cout<<"sm born by copy"<<"\n";
  }
  
  
  void Dif_alpha_connection::get_status(DictionaryDatum & d) const
  {
    //std::cout<<"get_status "<<"example"<<"\n";
    double_t t_weight;
    DictionaryDatum *d2 = new DictionaryDatum(new Dictionary); 
    def<double>((*d2),"local_id_synapse",local_id_synapse);
    target_->behav->get_parametrs(id_target,id_source,t_weight,(*d2));
    delete d2;
    //weight_=t_weight;
    //ConnectionHetWD::get_status(d);
    def<double_t>(d, names::weight, t_weight);
    def<double_t>(d, names::delay, Time(Time::step(delay_)).get_ms());
    def<double_t>(d, names::tau_a, tau_a_);
    def<double_t>(d, names::local_id_synapse, local_id_synapse);
    
  }
  
  void Dif_alpha_connection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
    //std::cout<<"set_status "<<"example 1"<<"\n";
    updateValue<double_t>(d, names::tau_a, tau_a_);
    updateValue<double_t>(d,names::local_id_synapse,local_id_synapse);
  }

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void Dif_alpha_connection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    //std::cout<<"set_status "<<"example 2"<<"\n";
    set_property<double_t>(d, "tau_as", p, tau_a_);
    set_property<double_t>(d,"local_id_synapses",p,local_id_synapse);
  }

  void Dif_alpha_connection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "tau_as");
    initialize_property_array(d,"local_id_synapses");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void Dif_alpha_connection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, "tau_as", tau_a_);
    append_property<double_t>(d,"local_id_synapses",local_id_synapse);
  }

} // of namespace nest
