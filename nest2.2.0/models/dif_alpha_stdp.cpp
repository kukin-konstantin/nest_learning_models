/*
 *  dif_alpha_stdp.cpp
 *
 */

#include "dif_alpha_stdp.h"
#include "network.h"
#include "connector_model.h"

namespace nest
{

  Dif_alpha_stdp::Dif_alpha_stdp() :
    ConnectionHetWD(),
    U_(0.5),
    D_(1100.0),
    F_(50.0),
    u_(0.5),
    R(1.0),
    tau_plus_(20.0),
    lambda_(0.01),
    alpha_(1.0),
    mu_plus_(1.0),
    mu_minus_(1.0),
    Wmax_(100.0),
    tau_a_(10.0),
    t_count(0)
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

  Dif_alpha_stdp::Dif_alpha_stdp(const Dif_alpha_stdp &rhs):
  ConnectionHetWD(rhs)
  {
    U_= rhs.U_;
    D_= rhs.D_;
    F_= rhs.F_;
    u_= rhs.u_;
    R= rhs.R;
    tau_plus_ = rhs.tau_plus_;
    lambda_ = rhs.lambda_;
    alpha_ = rhs.alpha_;
    mu_plus_ = rhs.mu_plus_;
    mu_minus_ = rhs.mu_minus_;
    Wmax_ = rhs.Wmax_;
    tau_a_=rhs.tau_a_;
    t_count=rhs.t_count;
    //std::cout<<"sm born by copy"<<"\n";
  }
  
  
  void Dif_alpha_stdp::get_status(DictionaryDatum & d) const
  {
    //std::cout<<"get_status "<<"example"<<"\n";
    double_t t_weight;
    target_->behav->get_parametrs(id_target,id_source,t_weight);
    //weight_=t_weight;
    //ConnectionHetWD::get_status(d);
    def<double_t>(d, names::weight, t_weight);
    def<double_t>(d, names::delay, Time(Time::step(delay_)).get_ms());
    def<double_t>(d, names::tau_a, tau_a_);
    def<double_t>(d, "U", U_);
    def<double_t>(d, "D", D_);
    def<double_t>(d, "F", F_);
    def<double_t>(d, "u", u_);
    def<double_t>(d, "R", R);
    def<double_t>(d, "tau_plus", tau_plus_);
    def<double_t>(d, "lambda", lambda_);
    def<double_t>(d, "alpha", alpha_);
    def<double_t>(d, "mu_plus", mu_plus_);
    def<double_t>(d, "mu_minus", mu_minus_);
    def<double_t>(d, "Wmax", Wmax_);
  }
  
  void Dif_alpha_stdp::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
    //std::cout<<"set_status "<<"example 1"<<"\n";
    updateValue<double_t>(d, names::tau_a, tau_a_);
    updateValue<double_t>(d, "U", U_);
    updateValue<double_t>(d, "D", D_);
    updateValue<double_t>(d, "F", F_);
    updateValue<double_t>(d, "u", u_);
    updateValue<double_t>(d, "R", R);
    updateValue<double_t>(d, "tau_plus", tau_plus_);
    updateValue<double_t>(d, "lambda", lambda_);
    updateValue<double_t>(d, "alpha", alpha_);
    updateValue<double_t>(d, "mu_plus", mu_plus_);
    updateValue<double_t>(d, "mu_minus", mu_minus_);
    updateValue<double_t>(d, "Wmax", Wmax_);
  }

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void Dif_alpha_stdp::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    //std::cout<<"set_status "<<"example 2"<<"\n";
    set_property<double_t>(d, "tau_as", p, tau_a_);
    set_property<double_t>(d, "Us", p, U_);
    set_property<double_t>(d, "Ds", p, D_);
    set_property<double_t>(d, "Fs", p, F_);
    set_property<double_t>(d, "us", p, u_);
    set_property<double_t>(d, "Rs", p, R);
    set_property<double_t>(d, "tau_pluss", p, tau_plus_);
    set_property<double_t>(d, "lambdas", p, lambda_);
    set_property<double_t>(d, "alphas", p, alpha_);
    set_property<double_t>(d, "mu_pluss", p, mu_plus_);
    set_property<double_t>(d, "mu_minuss", p, mu_minus_);
    set_property<double_t>(d, "Wmaxs", p, Wmax_);
  }

  void Dif_alpha_stdp::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "tau_as");
    
    initialize_property_array(d, "Us");
    initialize_property_array(d, "Ds");
    initialize_property_array(d, "Fs");
    initialize_property_array(d, "us");
    initialize_property_array(d, "Rs");
    initialize_property_array(d, "tau_pluss"); 
    initialize_property_array(d, "lambdas"); 
    initialize_property_array(d, "alphas"); 
    initialize_property_array(d, "mu_pluss"); 
    initialize_property_array(d, "mu_minuss");
    initialize_property_array(d, "Wmaxs");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void Dif_alpha_stdp::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, "tau_as", tau_a_);
    
    append_property<double_t>(d, "Us", U_);
    append_property<double_t>(d, "Ds", D_);
    append_property<double_t>(d, "Fs", F_);
    append_property<double_t>(d, "us", u_);
    append_property<double_t>(d, "Rs", R);
    append_property<double_t>(d, "tau_pluss", tau_plus_); 
    append_property<double_t>(d, "lambdas", lambda_); 
    append_property<double_t>(d, "alphas", alpha_); 
    append_property<double_t>(d, "mu_pluss", mu_plus_); 
    append_property<double_t>(d, "mu_minuss", mu_minus_);
    append_property<double_t>(d, "Wmaxs", Wmax_);
  }

} // of namespace nest
