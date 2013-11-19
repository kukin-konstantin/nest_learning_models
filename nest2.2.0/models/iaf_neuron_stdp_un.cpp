/*
 *  iaf_neuron_stdp_un.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
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
 */

#include "exceptions.h"
#include "iaf_neuron_stdp_un.h"
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

int nest::iaf_neuron_stdp_un::id_gl=-1; //nebodimo dlya togochtobi numeracia neuronov nacinalas s edinici
nest::RecordablesMap<nest::iaf_neuron_stdp_un> nest::iaf_neuron_stdp_un::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<iaf_neuron_stdp_un>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, &iaf_neuron_stdp_un::get_V_m_);
  }
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::iaf_neuron_stdp_un::Parameters_::Parameters_()
  : C_      (250.0    ),  // pF
    Tau_    ( 10.0    ),  // ms
    tau_syn_(  2.0    ),  // ms
    TauR_   (  2.0    ),  // ms
    U0_     (-70.0    ),  // mV
    V_reset_(-70.0-U0_),  // mV, rel to U0_
    Theta_  (-55.0-U0_),  // mV, rel to U0_
    I_e_    (  0.0    )   // pA
{}

nest::iaf_neuron_stdp_un::State_::State_()
  : y0_(0.0),
    y1_(0.0),
    y2_(0.0),
    y3_(0.0),  
    r_ (0)
{}

nest::iaf_neuron_stdp_un::~iaf_neuron_stdp_un()
{
  delete behav_loc;
  delete behav_learning_loc;
  std::cout<<"iaf_neuron_stdp_un destruct counter="<<iaf_neuron_stdp_un::id_gl<<"\n";
  iaf_neuron_stdp_un::id_gl--;
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::iaf_neuron_stdp_un::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, names::E_L, U0_);   // Resting potential
  def<double>(d, names::I_e, I_e_);
  def<double>(d, names::V_th, Theta_+U0_); // threshold value
  def<double>(d, names::V_reset, V_reset_+U0_);
  def<double>(d, names::C_m, C_);
  def<double>(d, names::tau_m, Tau_);
  def<double>(d, names::tau_syn, tau_syn_);
  def<double>(d, names::t_ref, TauR_);
}

double nest::iaf_neuron_stdp_un::Parameters_::set(const DictionaryDatum& d)
{
  // if U0_ is changed, we need to adjust all variables defined relative to U0_
  const double ELold = U0_;
  updateValue<double>(d, names::E_L, U0_);
  const double delta_EL = U0_ - ELold;

  if(updateValue<double>(d, names::V_reset, V_reset_))
    V_reset_ -= U0_;   // here we use the new U0_, no need for adjustments
  else
    V_reset_ -= delta_EL;  // express relative to new U0_

  if (updateValue<double>(d, names::V_th, Theta_)) 
    Theta_ -= U0_;
  else
    Theta_ -= delta_EL;  // express relative to new U0_

  updateValue<double>(d, names::I_e, I_e_);
  updateValue<double>(d, names::C_m, C_);
  updateValue<double>(d, names::tau_m, Tau_);
  updateValue<double>(d, names::tau_syn, tau_syn_);
  updateValue<double>(d, names::t_ref, TauR_);

  if ( V_reset_ >= Theta_ )
    throw BadProperty("Reset potential must be smaller than threshold.");
    
  if ( C_ <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  if ( Tau_ <= 0 || tau_syn_ <= 0 || TauR_ <= 0 )
    throw BadProperty("All time constants must be strictly positive.");

  if ( Tau_ == tau_syn_ )
    throw BadProperty("Membrane and synapse time constant(s) must differ."
		      "See note in documentation.");

  return delta_EL;
}

void nest::iaf_neuron_stdp_un::State_::get(DictionaryDatum &d, const Parameters_& p) const
{
  def<double>(d, names::V_m, y3_ + p.U0_); // Membrane potential
}

void nest::iaf_neuron_stdp_un::State_::set(const DictionaryDatum& d, const Parameters_& p, double delta_EL)
{
  if ( updateValue<double>(d, names::V_m, y3_) )
    y3_ -= p.U0_;
  else
    y3_ -= delta_EL;
}

nest::iaf_neuron_stdp_un::Buffers_::Buffers_(iaf_neuron_stdp_un &n)
  : logger_(n)
{}

nest::iaf_neuron_stdp_un::Buffers_::Buffers_(const Buffers_ &, iaf_neuron_stdp_un &n)
  : logger_(n)
{}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::iaf_neuron_stdp_un::iaf_neuron_stdp_un()
  : Archiving_Node(), 
    P_(), 
    S_(),
    B_(*this)
{
  recordablesMap_.create();
  behav=new Iaf_neuron_sm_behavior;
  //
  DictionaryDatum *d = new DictionaryDatum(new Dictionary);
  behav_learning=new STDP_learning_behavior((*d));
  std::cout<<"non-copy iaf_neuron_stdp_un::id_gl "<<iaf_neuron_stdp_un::id_gl<<"\n"; 
  behav_loc= dynamic_cast<Iaf_neuron_sm_behavior *>(behav);
  behav_learning_loc=dynamic_cast<STDP_learning_behavior *>(behav_learning);
  delete d;
  //vnesti suda Behavior_learning
  //nest::iaf_neuron_sm::sh_gl=0;
}

nest::iaf_neuron_stdp_un::iaf_neuron_stdp_un(const iaf_neuron_stdp_un& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
  behav=new Iaf_neuron_sm_behavior;
  DictionaryDatum *d = new DictionaryDatum(new Dictionary);
  behav_learning=new STDP_learning_behavior((*d));
  //std::cout<<"iaf born by copy "<<"behav too"<<"\n"; 
  behav_loc= dynamic_cast<Iaf_neuron_sm_behavior *>(behav);
  behav_learning_loc=dynamic_cast<STDP_learning_behavior *>(behav_learning);
  /*learn function*/
  //vnesti suda Behavior_learning
  std::stringstream r;
  iaf_neuron_stdp_un::id_gl++;
  id_neuron=iaf_neuron_stdp_un::id_gl;
  std::cout<<"copy iaf_neuron_stdp_un::id_gl "<<iaf_neuron_stdp_un::id_gl<<"\n";
  r << iaf_neuron_stdp_un::id_gl;
  //r<<++nest::iaf_neuron_sm::sh_gl;
  //id_gl=nest::iaf_neuron_sm::sh_gl;
  std::string s_directory_learn_set_desired="/home/neuron/animat_work/learn_set/desired/neuron_"+r.str()+".txt";
  std::string s_directory_learn_set_input="/home/neuron/animat_work/learn_set/input/neuron_";
  //std::cout<<"s_tmp "<<s_tmp<<"\n";
  behav_learning_loc->read_learn_set_desired_from_file(s_directory_learn_set_desired);//perepisat vizov
  //behav_learning_loc->read_learn_set_input_from_files(s_directory_learn_set_input,iaf_neuron_stdp_un::id_gl); //ano
  /*learn function*/
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::iaf_neuron_stdp_un::init_state_(const Node& proto)
{
  const iaf_neuron_stdp_un& pr = downcast<iaf_neuron_stdp_un>(proto);
  S_ = pr.S_;
}

void nest::iaf_neuron_stdp_un::init_buffers_()
{
  B_.spikes_.clear();    // includes resize
  B_.currents_.clear();  // include resize
  B_.logger_.reset(); // includes resize
  Archiving_Node::clear_history();
}

void nest::iaf_neuron_stdp_un::calibrate()
{
  B_.logger_.init();

  const double h = Time::get_resolution().get_ms(); 

  // these P are independent
  V_.P11_ = V_.P22_ = std::exp(-h/P_.tau_syn_);
  V_.P33_ = std::exp(-h/P_.Tau_);
  V_.P21_ = h * V_.P11_;
  
  // these depend on the above. Please do not change the order.
  V_.P30_ = 1/P_.C_*(1-V_.P33_)*P_.Tau_;
  V_.P31_ = 1/P_.C_ * ((V_.P11_-V_.P33_)/(-1/P_.tau_syn_- -1/P_.Tau_)- h*V_.P11_)
    /(-1/P_.Tau_ - -1/P_.tau_syn_);
  V_.P32_ = 1/P_.C_*(V_.P33_-V_.P11_)/(-1/P_.Tau_ - -1/P_.tau_syn_);
  V_.PSCInitialValue_=1.0 * numerics::e/P_.tau_syn_;


  // TauR specifies the length of the absolute refractory period as 
  // a double_t in ms. The grid based iaf_neuron_stdp_un can only handle refractory
  // periods that are integer multiples of the computation step size (h).
  // To ensure consistency with the overall simulation scheme such conversion
  // should be carried out via objects of class nest::Time. The conversion 
  // requires 2 steps:
  //     1. A time object is constructed defining representation of 
  //        TauR in tics. This representation is then converted to computation time
  //        steps again by a strategy defined by class nest::Time.
  //     2. The refractory time in units of steps is read out get_steps(), a member
  //        function of class nest::Time.
  //
  // The definition of the refractory period of the iaf_neuron_stdp_un is consistent 
  // the one of iaf_psc_alpha_ps.
  //
  // Choosing a TauR that is not an integer multiple of the computation time 
  // step h will lead to accurate (up to the resolution h) and self-consistent
  // results. However, a neuron model capable of operating with real valued spike
  // time may exhibit a different effective refractory time.

  V_.RefractoryCounts_ = Time(Time::ms(P_.TauR_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::iaf_neuron_stdp_un::update(Time const & origin, const long_t from, const long_t to)
{
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);
  const double_t h = Time::get_resolution().get_ms();
  /*add KAK*/
  DictionaryDatum *d = new DictionaryDatum(new Dictionary);
  //behav_learning_loc->get((*d));
  get_status_user((*d)); //KAK
  //get_status((*d));
  double_t t=origin.get_ms();
  bool trig_seal_;
  bool trig_learn_oper_;
  updateValue<bool>((*d), names::trig_seal,trig_seal_);// d->trig_seal_
  updateValue<bool>((*d), names::trig_learn_oper,trig_learn_oper_);
  if (trig_seal_)
  {
    //std::cout<<"trig_seal\n";
    //std::cout<<"grad "<<"\n";
    trig_seal_=0;
    //std::cout<<id_neuron<<"behav_learning_loc "<<behav_learning_loc<<"\n";
    def<bool>((*d), names::trig_seal, trig_seal_);;
    behav_learning_loc->reset_iterator_desired();
    behav_learning_loc->set((*d));
  }
  /*add KAK*/
  for ( long_t lag = from ; lag < to ; ++lag )
    {
      if ( S_.r_ == 0 )
	{
	  // neuron not refractory
	  S_.y3_ = V_.P30_*(S_.y0_ + P_.I_e_) + V_.P31_*S_.y1_ + V_.P32_*S_.y2_ + V_.P33_*S_.y3_;
	}
      else // neuron is absolute refractory
	--S_.r_;

      // alpha shape PSCs
      S_.y2_  = V_.P21_*S_.y1_ + V_.P22_ * S_.y2_;
      S_.y1_ *= V_.P11_;
    
      // Apply spikes delivered in this step: The spikes arriving at T+1 have an 
      // immediate effect on the state of the neuron
      S_.y1_ += V_.PSCInitialValue_* B_.spikes_.get_value(lag);   
    
      def<double>((*d), names::time,t+h*lag);
      //std::cout<<id_neuron<<" time="<<t+h*lag<<"\n";
      // threshold crossing;
      if (trig_learn_oper_) //KAK 
      {
	std::vector<double > currents;
	//std::cout<<id_neuron<<"behav_learning_loc "<<behav_learning_loc<<"\n";
	//std::cout<<id_neuron<<" behav_learning_loc->learn(behav_loc,(*d),currents)"<<"\n";
	if (behav_learning_loc->learn(behav_loc,(*d),currents))
	{
	  S_.r_ = V_.RefractoryCounts_;
	  S_.y3_= P_.V_reset_; ;
	  // A supra-threshold membrane potential should never be observable.
	  // The reset at the time of threshold crossing enables accurate integration
	  // independent of the computation step size, see [2,3] for details.
  	  //std::cout<<"origin "<<origin.get_ms()<<"\t"<<"lag*h "<<Time::get_resolution().get_ms()*lag<<"\t";//add tmp
	  //std::cout<<"spike time "<<origin.get_ms()+Time::get_resolution().get_ms()*lag<<"\n";//add tmp
	  //std::cout<<"learning id_neuron="<<id_neuron<<" time="<<t+h*lag<<"\n";
	  set_spiketime(Time::step(origin.get_steps()+lag+1));
	  SpikeEvent se;
	  network()->send(*this, se, lag);
	}
	if (S_.y3_ >= P_.Theta_)
	{
	  S_.r_ = V_.RefractoryCounts_;
	  S_.y3_= P_.V_reset_;
	  set_spiketime(Time::step(origin.get_steps()+lag+1));
	  SpikeEvent se;
	  network()->send(*this, se, lag);
	}
      }
      else
      {
	if (S_.y3_ >= P_.Theta_)
	{
	  //std::cout<<"living id_neuron="<<id_neuron<<" time="<<t+h*lag<<"\n";
	  S_.r_ = V_.RefractoryCounts_;
	  S_.y3_= P_.V_reset_; 
	  // A supra-threshold membrane potential should never be observable.
	  // The reset at the time of threshold crossing enables accurate integration
	  // independent of the computation step size, see [2,3] for details.
  	  //std::cout<<"origin "<<origin.get_ms()<<"\t"<<"lag*h "<<Time::get_resolution().get_ms()*lag<<"\t";//add tmp
	  //std::cout<<"spike time "<<origin.get_ms()+Time::get_resolution().get_ms()*lag<<"\n";//add tmp
	  set_spiketime(Time::step(origin.get_steps()+lag+1));
	  SpikeEvent se;
	  network()->send(*this, se, lag);
	}
      }

      // set new input current
      S_.y0_ = B_.currents_.get_value(lag);

      // voltage logging
      B_.logger_.record_data(origin.get_steps()+lag);
    }
    delete d;
}                           
                     
void nest::iaf_neuron_stdp_un::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  B_.spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
		       e.get_weight() * e.get_multiplicity() );
}

void nest::iaf_neuron_stdp_un::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			 w * c);
}

void nest::iaf_neuron_stdp_un::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}
