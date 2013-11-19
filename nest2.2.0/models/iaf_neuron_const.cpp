/*
 *  iaf_neuron_const.cpp
 *
 *
 */

#include "exceptions.h"
#include "iaf_neuron_const.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"

#include <limits>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<nest::iaf_neuron_const> nest::iaf_neuron_const::recordablesMap_;


namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<iaf_neuron_const>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, &iaf_neuron_const::get_V_m_);
  }
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::iaf_neuron_const::Parameters_::Parameters_()
  : C_      (200.0    ),  // pF
    Tau_    ( 10.0    ),  // ms
    TauR_   (  2.0    ),  // ms
    E_L_     (-70.0    ),  // mV
    V_reset_(-70.0),
    V_th_  (-54.0),  // mV
    I_e_    (  0.0    )   // pA
    //n_in_(0)
    
{}

nest::iaf_neuron_const::State_::State_()
  : 
  v_(-70.0),
    I_(0.0),  
    r_ (0)
{}




/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::iaf_neuron_const::Parameters_::get(DictionaryDatum &d) const
{
  
  def<double>(d, names::E_L, E_L_);   // Resting potential
  def<double>(d, names::V_reset, V_reset_);   // Resting potential
  def<double>(d, names::I_e, I_e_);
  def<double>(d, names::V_th, V_th_); // threshold value
  def<double>(d, names::C_m, C_);
  def<double>(d, names::tau_m, Tau_);
  def<double>(d, names::t_ref, TauR_);
}

void nest::iaf_neuron_const::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, names::E_L, E_L_);
  updateValue<double>(d, names::V_reset, V_reset_);
  updateValue<double>(d, names::V_th, V_th_); 
  updateValue<double>(d, names::I_e, I_e_);
  updateValue<double>(d, names::C_m, C_);
  updateValue<double>(d, names::tau_m, Tau_);
  updateValue<double>(d, names::t_ref, TauR_);
  if ( C_ <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  if ( Tau_ <= 0 || TauR_ <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
}

void nest::iaf_neuron_const::State_::get(DictionaryDatum &d, const Parameters_& p) const
{
  def<double>(d, names::V_m, v_); // Membrane potential
}

void nest::iaf_neuron_const::State_::set(const DictionaryDatum& d, const Parameters_& p)
{
  updateValue<double>(d, names::V_m, v_);
}


nest::iaf_neuron_const::Buffers_::Buffers_(iaf_neuron_const &n)
  : logger_(n)
{}

nest::iaf_neuron_const::Buffers_::Buffers_(const Buffers_ &, iaf_neuron_const &n)
  : logger_(n)
{}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::iaf_neuron_const::iaf_neuron_const()
  : Archiving_Node(), 
    P_(), 
    S_(),
    B_(*this)
{
  recordablesMap_.create();
}

nest::iaf_neuron_const::iaf_neuron_const(const iaf_neuron_const& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{}



/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::iaf_neuron_const::init_state_(const Node& proto)
{
  const iaf_neuron_const& pr = downcast<iaf_neuron_const>(proto);
  S_ = pr.S_;
}



void nest::iaf_neuron_const::init_buffers_()
{
  B_.spikes_.clear();    // includes resize
  B_.currents_.clear();  // include resize
  B_.logger_.reset(); // includes resize
  Archiving_Node::clear_history();
}



void nest::iaf_neuron_const::calibrate()
{
  B_.logger_.init();
  V_.RefractoryCounts_ = Time(Time::ms(P_.TauR_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0); 
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */



void nest::iaf_neuron_const::update(Time const & origin, const long_t from, const long_t to)
{
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);
  const double_t h = Time::get_resolution().get_ms();
  
  double_t v_old;
  for ( long_t lag = from ; lag < to ; ++lag )
  {
    if ( S_.r_ == 0 ) // neuron not refractory
    {
      v_old = S_.v_;
      S_.v_ +=h*( ((P_.E_L_-v_old)/P_.Tau_) + (P_.I_e_/P_.C_) + (S_.I_/P_.C_) )+B_.spikes_.get_value(lag);
      //S_.v_ +=h*( ((P_.E_L_-v_old)/P_.Tau_) + (P_.I_e_/P_.C_) + (S_.I_/P_.C_) );
    }
    else // neuron is absolute refractory
      --S_.r_;
  
    if (S_.v_ >= P_.V_th_)
    {
	S_.r_ = V_.RefractoryCounts_;
	S_.v_= P_.V_reset_;   
	set_spiketime(Time::step(origin.get_steps()+lag+1)); // change -1
	//S_.I_=0;
	//std::cout<<"fun time spike from neuron current"<<S_.I_<<"\n";// change -1
	SpikeEvent se;
	network()->send(*this, se, lag);
    }
    //B_.spikes_.get_value(lag);// Need for clean buffer. It's necessary for working log-device wrapper.
    //std::cout<<"neuron id="<<get_gid()<<"\told s.I="<<S_.I_<<"\t";
    S_.I_=B_.currents_.get_value(lag);
    //std::cout<<"neuron id="<<get_gid()<<"\tnew s.I="<<S_.I_<<"\n";
    //S_.I_ =S_.I_+get_sum_I(t+h*lag); 
    // voltage logging
    B_.logger_.record_data(origin.get_steps()+lag);
  }
}                           
                     
void nest::iaf_neuron_const::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  B_.spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
		       e.get_weight() * e.get_multiplicity() );
  //std::cout<<"spike handle "<<e.get_stamp().get_ms()<<"\n";
}

void nest::iaf_neuron_const::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			 w * c);
}

void nest::iaf_neuron_const::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}
