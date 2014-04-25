/*
 *  iaf_neuron_dif_alpha_stdp.cpp
 *
 *
 */

#include "exceptions.h"
#include "iaf_neuron_dif_alpha_stdp.h"
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

int nest::iaf_neuron_dif_alpha_stdp::id_gl=-1; //nebodimo dlya togochtobi numeracia neuronov nacinalas s edinici
nest::RecordablesMap<nest::iaf_neuron_dif_alpha_stdp> nest::iaf_neuron_dif_alpha_stdp::recordablesMap_;


namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<iaf_neuron_dif_alpha_stdp>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, &iaf_neuron_dif_alpha_stdp::get_V_m_);
  }
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::iaf_neuron_dif_alpha_stdp::Parameters_::Parameters_()
  : C_      (200.0    ),  // pF
    Tau_    ( 10.0    ),  // ms
    TauR_   (  2.0    ),  // ms
    E_L_     (-70.0    ),  // mV
    V_reset_ (0.0),
    V_th_  (-54.0),  // mV
    I_e_    (  0.0    ),   // pA
    switch_work_mode_(false)
    //n_in_(0)
    
{}

nest::iaf_neuron_dif_alpha_stdp::State_::State_()
  : 
  v_(-70.0),
    I_(0.0),  
    r_ (0)
{
  time_list_.clear();
}

nest::iaf_neuron_dif_alpha_stdp::~iaf_neuron_dif_alpha_stdp()
{
  /*std::cout<<"destruct learn_set_input "<<"\n";
   for (int i=0;i!=behav_learning_loc->v_learn_set_input.size();i++)
   {
     for (int j=0;j!=behav_learning_loc->v_learn_set_input[i].size();j++)
     {
       for (int k=0;k!=behav_learning_loc->v_learn_set_input[i][j].size();k++)
       {
	  std::cout<<behav_learning_loc->v_learn_set_input[i][j][k]<<"\t";
       }
        std::cout<<"\n";
     }
     std::cout<<"\n\n";
   }*/
  delete behav_loc;
  std::cout<<"iaf_neuron_dif_alpha_stdp destruct counter="<<iaf_neuron_dif_alpha_stdp::id_gl<<"\n";
  iaf_neuron_dif_alpha_stdp::id_gl--;
}



/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::iaf_neuron_dif_alpha_stdp::Parameters_::get(DictionaryDatum &d) const
{
  
  def<double>(d, names::E_L, E_L_);   // Resting potential
  def<double>(d, names::I_e, I_e_);
  def<double>(d, names::V_th, V_th_); // threshold value
  def<double>(d, names::V_reset, V_reset_);
  def<double>(d, names::C_m, C_);
  def<double>(d, names::tau_m, Tau_);
  def<double>(d, names::t_ref, TauR_);
  def<double>(d,names::switch_work_mode,switch_work_mode_);
  //dobavit peremennie obuchenia
}

void nest::iaf_neuron_dif_alpha_stdp::Parameters_::set(const DictionaryDatum& d)
{
  
  updateValue<double>(d, names::E_L, E_L_);
  updateValue<double>(d, names::V_th, V_th_);
  updateValue<double>(d, names::V_reset, V_reset_);
  updateValue<double>(d, names::I_e, I_e_);
  updateValue<double>(d, names::C_m, C_);
  updateValue<double>(d, names::tau_m, Tau_);
  updateValue<double>(d, names::t_ref, TauR_);
  updateValue<double>(d, names::switch_work_mode, switch_work_mode_);
  //dobavit peremennie obuchenia
  if ( C_ <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  if ( Tau_ <= 0 || TauR_ <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
}

void nest::iaf_neuron_dif_alpha_stdp::State_::get(DictionaryDatum &d, const Parameters_& p) const
{
  def<double>(d, names::V_m, v_); // Membrane potential
  
  ArrayDatum time_list_ad(time_list_);
  def<ArrayDatum>(d,"time_list", time_list_ad);
}

void nest::iaf_neuron_dif_alpha_stdp::State_::set(const DictionaryDatum& d, const Parameters_& p)
{
  updateValue<double>(d, names::V_m, v_);
  
  std::vector<double> tau_tmp;
  updateValue<std::vector<double> >(d, "time_list", time_list_);
}


nest::iaf_neuron_dif_alpha_stdp::Buffers_::Buffers_(iaf_neuron_dif_alpha_stdp &n)
  : logger_(n)
{}

nest::iaf_neuron_dif_alpha_stdp::Buffers_::Buffers_(const Buffers_ &, iaf_neuron_dif_alpha_stdp &n)
  : logger_(n)
{}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::iaf_neuron_dif_alpha_stdp::iaf_neuron_dif_alpha_stdp()
  : Archiving_Node(), 
    P_(), 
    S_(),
    B_(*this)
{
  recordablesMap_.create();
  behav=new Iaf_neuron_dif_alpha_stdp_behavior;
  //
  DictionaryDatum *d = new DictionaryDatum(new Dictionary);
  std::cout<<"non-copy iaf_neuron_dif_alpha_stdp::id_gl "<<iaf_neuron_dif_alpha_stdp::id_gl<<"\n"; 
  behav_loc= dynamic_cast<Iaf_neuron_dif_alpha_stdp_behavior *>(behav);
  //vnesti suda Behavior_learning
  //nest::iaf_neuron_dif_alpha_stdp::sh_gl=0;
}

nest::iaf_neuron_dif_alpha_stdp::iaf_neuron_dif_alpha_stdp(const iaf_neuron_dif_alpha_stdp& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
  behav=new Iaf_neuron_dif_alpha_stdp_behavior;
  DictionaryDatum *d = new DictionaryDatum(new Dictionary);
  //std::cout<<"iaf born by copy "<<"behav too"<<"\n"; 
  behav_loc= dynamic_cast<Iaf_neuron_dif_alpha_stdp_behavior *>(behav);
  /*learn function*/
  //vnesti suda Behavior_learning
  std::stringstream r;
  iaf_neuron_dif_alpha_stdp::id_gl++;
  id_neuron=iaf_neuron_dif_alpha_stdp::id_gl;
  /*learn function*/
}




/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::iaf_neuron_dif_alpha_stdp::init_state_(const Node& proto)
{
  const iaf_neuron_dif_alpha_stdp& pr = downcast<iaf_neuron_dif_alpha_stdp>(proto);
  S_ = pr.S_;
}



void nest::iaf_neuron_dif_alpha_stdp::init_buffers_()
{
  B_.spikes_.clear();    // includes resize
  B_.currents_.clear();  // include resize
  B_.logger_.reset(); // includes resize
  Archiving_Node::clear_history();
}



void nest::iaf_neuron_dif_alpha_stdp::calibrate()
{
  B_.logger_.init();
  V_.rng_ = net_->get_rng(get_thread());
  V_.RefractoryCounts_ = Time(Time::ms(P_.TauR_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);
  const double_t h = Time::get_resolution().get_ms();
  V_.C1=std::exp(-h/P_.Tau_);
  V_.C2=(P_.Tau_/P_.C_)*(1-V_.C1);
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

double nest::iaf_neuron_dif_alpha_stdp::normalization(double_t v)
{
  double k=0.1;//coff learning
  double tresh=1.0;
  double u=(v-P_.E_L_)/(P_.V_th_);
  if (u<0)
    u=0;
  double z=(u-tresh)/k;
  return z;
}

bool nest::iaf_neuron_dif_alpha_stdp::emission_procedure(double_t v)
{
  double Th=1.0;
  double k=0.1;
  const double_t h = Time::get_resolution().get_ms();
  double rand=V_.rng_->drand();
  double lambd_val=std::exp((v-Th)/k);
  //double z=normalization(v);
  double z=v;
  //double sigm=1.0/(1.0+std::exp(-z));
  //double sigm=z;
  double sigm=1.0-std::exp(-lambd_val*h);
  if (rand<=sigm)
  {
    return true;
  }
  else
  {
    return false;
  }
}

double nest::iaf_neuron_dif_alpha_stdp::facilitate(double w, double kplus,double lambda,double Wmax)
{
  double t_w_new=w+lambda*kplus;
  if (Wmax<t_w_new)
  {
    return Wmax;
  }
  else
  {
    return t_w_new;
  }
}

void nest::iaf_neuron_dif_alpha_stdp::stdp_procedure(double time_post_spike)
{
  std::map<std::pair<index,index>,struct_iaf_neuron_dif_alpha_stdp>::iterator it=behav_loc->val.begin();
        //std::cout<<"time_post_spike="<<time_post_spike<<"\n";//rubbish
	while (it!=behav_loc->val.end())
	{
	    double t_spike;
	    double t_Wmax=(*it).second.Wmax;
	    double t_lambda=(*it).second.lambda;
	    double t_tau_plus=(*it).second.tau_plus;
	    /*if (!((*it).second.base_val.v_time_spike_times.empty()))//rubbish
	    {
	      std::cout<<"Wmax="<<t_Wmax<<"\t";//rubbish
	    }*/
	    for (int i=0;i!=(*it).second.base_val.v_time_spike_times.size();i++)
	    {
		t_spike=(*it).second.base_val.v_time_spike_times[i];
		double minus_dt=t_spike-time_post_spike;
		double t_weight=(*it).second.base_val.weight;
		if (t_weight>=0.0)
		{
		  //if ((std::abs(minus_dt)>2.0))
		  //std::cout<<"minus_dt="<<minus_dt<<"\t";//rubbish
		  if ((std::abs(minus_dt)>2.0)&&(std::abs(minus_dt)<3*20.0))
		  {
		  (*it).second.base_val.weight=facilitate(t_weight,std::exp(minus_dt / t_tau_plus),t_lambda,t_Wmax);
		  }
		}
	    }
	    /*if (!((*it).second.base_val.v_time_spike_times.empty()))//rubbish
	    {
	      std::cout<<"\n\n";//rubbish
	    }*/
	    it++;
	}
}

void nest::iaf_neuron_dif_alpha_stdp::update(Time const & origin, const long_t from, const long_t to)
{
  //bool ok=true;
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);
  double h = Time::get_resolution().get_ms();
  double_t v_old;
  double t=origin.get_ms();
  //std::cout<<"t="<<t<<"\t";
  //std::cout<<S_.v_<<"\t";
  //std::cout<<std::fixed<<std::setprecision(2);
  for ( long_t lag = from ; lag < to ; ++lag )
   {
     std::cout<<"t="<<t<<"\t";
     std::cout<<"h*lag="<<h*lag<<"\t";
     double z=t+h*lag;
     std::cout<<"z-t="<<z-t<<"\t";
     //std::cout<<"before S_.v="<<S_.v_<<"\t";
     //std::cout<<"before S_.I="<<S_.I_<<"\t";
     if (S_.r_!=0) // neuron is absolute refractory
     {
       --S_.r_;
     }
     else if ( S_.r_ == 0 ) // neuron not refractory
     {
	  v_old = S_.v_;
	  bool b_spiking=false;
	  if (P_.switch_work_mode_) //true - learning
	  {
	    b_spiking=true;
	  }
	  else
	  {
	    if (S_.v_ >= P_.V_th_)
	    {
	      b_spiking=true;
	    }
	  }
	  if (b_spiking) //change there
	  {
		stdp_procedure(t+h*lag);
		S_.r_ = V_.RefractoryCounts_;
		S_.v_= P_.V_reset_;   
		set_spiketime(Time::step(origin.get_steps()+lag)); // change -1
		SpikeEvent se;
		network()->send(*this, se, lag);
 		clear_spike_history(t+h*lag);
	  }
	  else
	  {
	    //S_.v_ +=h*( ((P_.E_L_-v_old)/P_.Tau_) + (P_.I_e_/P_.C_) + (S_.I_/P_.C_)) ;
	    S_.v_=V_.C1*v_old+V_.C2*(P_.I_e_+S_.I_);
	    //std::cout<<"after S_.v="<<S_.v_<<"\t";
	    
	  }
     }
     S_.I_=B_.currents_.get_value(lag);
     S_.I_ =S_.I_+get_sum_I(t+h*lag);
     //std::cout<<"after I="<<S_.I_<<"\t";
     B_.logger_.record_data(origin.get_steps()+lag); 
     
	//std::cout<<"\n";
      }
      //S_.I_=0;//change
      //std::cout<<"\n"; // delete after use
}                           
                     
void nest::iaf_neuron_dif_alpha_stdp::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  B_.spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
		       e.get_weight() * e.get_multiplicity() );
  //std::cout<<"spike handle "<<e.get_stamp().get_ms()<<"\n";
}

void nest::iaf_neuron_dif_alpha_stdp::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			 w * c);
}

void nest::iaf_neuron_dif_alpha_stdp::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}
