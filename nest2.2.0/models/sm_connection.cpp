/*
 *  sm_connection.cpp
 *
 */

#include "sm_connection.h"
#include "network.h"
#include "connector_model.h"

namespace nest
{

  SMConnection::SMConnection() :
    ConnectionHetWD(),
    g_(1.0),
    E_(50.0),
    t_count(0.0)
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

  SMConnection::SMConnection(const SMConnection &rhs):
  ConnectionHetWD(rhs)
  {
    g_=rhs.g_;
    E_=rhs.E_;
    t_count=rhs.t_count;
    //std::cout<<"sm born by copy"<<"\n";
  }
  
  
  void SMConnection::get_status(DictionaryDatum & d) const
  {
    //std::cout<<"get_status "<<"example"<<"\n";
    double_t t_weight;
    target_->behav->get_parametrs(id_target,id_source,t_weight,(*d));
    //weight_=t_weight;
    //ConnectionHetWD::get_status(d);
    def<double_t>(d, names::weight, t_weight);
    def<double_t>(d, names::delay, Time(Time::step(delay_)).get_ms());
    def<double_t>(d, names::g, g_);
    def<double_t>(d, names::E, E_);
    
  }
  
  void SMConnection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
    std::cout<<"set_status "<<"example 1"<<"\n";
    updateValue<double_t>(d, names::g, g_);
    updateValue<double_t>(d, names::E, E_);
  }

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void SMConnection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    std::cout<<"set_status "<<"example 2"<<"\n";
    set_property<double_t>(d, "gs", p, g_);
    set_property<double_t>(d, "Es", p, E_);
  }

  void SMConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "gs");
    initialize_property_array(d, "Es");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void SMConnection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, "gs", g_);
    append_property<double_t>(d, "Es", E_);
  }

} // of namespace nest
