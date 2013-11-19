
/*
 *  stdp_tsodyk.cpp
 *
 *
 */
#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "stdp_tsodyk.h"
#include "event.h"


namespace nest
{

  STDP_Tsodyk_Connection::STDP_Tsodyk_Connection() :
    ConnectionHetWD(),
	tau_psc_(3.0),
    tau_fac_(0.0),
    tau_rec_(800.0),
    U_(0.5),
    x_(1.0),
    y_(0.0),
    u_(0.0),
    tau_plus_(20.0),
    lambda_(0.01),
    alpha_(1.0),
    mu_plus_(1.0),
    mu_minus_(1.0),
    Wmax_(100.0),
    Kplus_(0.0)
  { }


  STDP_Tsodyk_Connection::STDP_Tsodyk_Connection(const STDP_Tsodyk_Connection &rhs) :
    ConnectionHetWD(rhs)
  {
    tau_psc_= rhs.tau_psc_;
    tau_fac_= rhs.tau_fac_;
    tau_rec_= rhs.tau_rec_;
    U_= rhs.U_;
    x_= rhs.x_;
    y_= rhs.y_;
    u_= rhs.u_;
    tau_plus_ = rhs.tau_plus_;
    lambda_ = rhs.lambda_;
    alpha_ = rhs.alpha_;
    mu_plus_ = rhs.mu_plus_;
    mu_minus_ = rhs.mu_minus_;
    Wmax_ = rhs.Wmax_;
    Kplus_ = rhs.Kplus_;
  }

  void STDP_Tsodyk_Connection::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);

	def<double_t>(d, "U", U_);
    def<double_t>(d, "tau_psc", tau_psc_);
    def<double_t>(d, "tau_rec", tau_rec_);
    def<double_t>(d, "tau_fac", tau_fac_);
    def<double_t>(d, "x", x_);
    def<double_t>(d, "y", y_);
    def<double_t>(d, "u", u_);
    def<double_t>(d, "tau_plus", tau_plus_);
    def<double_t>(d, "lambda", lambda_);
    def<double_t>(d, "alpha", alpha_);
    def<double_t>(d, "mu_plus", mu_plus_);
    def<double_t>(d, "mu_minus", mu_minus_);
    def<double_t>(d, "Wmax", Wmax_);
  }

  void STDP_Tsodyk_Connection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
	
	updateValue<double_t>(d, "U", U_);
    updateValue<double_t>(d, "tau_psc", tau_psc_);
    updateValue<double_t>(d, "tau_rec", tau_rec_);
    updateValue<double_t>(d, "tau_fac", tau_fac_);

    double_t x = x_;
    double_t y = y_;
    updateValue<double_t>(d, "x", x);
    updateValue<double_t>(d, "y", y);

    if (x + y > 1.0)
    {
      cm.network().message(SLIInterpreter::M_ERROR,
			   "TsodyksConnection::set_status()", "x + y must be <= 1.0.");
      throw BadProperty();
    }
    else
    {
      x_ = x;
      y_ = y;
    }

    updateValue<double_t>(d, "u", u_);

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
  void STDP_Tsodyk_Connection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
	
	set_property<double_t>(d, "Us", p, U_);
    set_property<double_t>(d, "tau_pscs", p, tau_psc_);
    set_property<double_t>(d, "tau_facs", p, tau_fac_);
    set_property<double_t>(d, "tau_recs", p, tau_rec_);


    double_t x = x_;
    double_t y = y_;
    set_property<double_t>(d, "xs", p, x);
    set_property<double_t>(d, "ys", p, y);

    if (x + y > 1.0)
    {
      cm.network().message(SLIInterpreter::M_ERROR, 
			   "TsodyksConnection::set_status()", 
			   "x + y must be <= 1.0.");
      throw BadProperty();
    }
    else
    {
      x_ = x;
      y_ = y;
    }

    set_property<double_t>(d, "us", p, u_);

    set_property<double_t>(d, "tau_pluss", p, tau_plus_);
    set_property<double_t>(d, "lambdas", p, lambda_);
    set_property<double_t>(d, "alphas", p, alpha_);
    set_property<double_t>(d, "mu_pluss", p, mu_plus_);
    set_property<double_t>(d, "mu_minuss", p, mu_minus_);
    set_property<double_t>(d, "Wmaxs", p, Wmax_);
  }

  void STDP_Tsodyk_Connection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
	
	initialize_property_array(d, "Us"); 
    initialize_property_array(d, "tau_pscs");
    initialize_property_array(d, "tau_facs");
    initialize_property_array(d, "tau_recs");  
    initialize_property_array(d, "xs"); 
    initialize_property_array(d, "ys");
    initialize_property_array(d, "us");

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
  void STDP_Tsodyk_Connection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
	
	append_property<double_t>(d, "Us", U_); 
    append_property<double_t>(d, "tau_pscs", tau_psc_);
    append_property<double_t>(d, "tau_facs", tau_fac_);
    append_property<double_t>(d, "tau_recs", tau_rec_);  
    append_property<double_t>(d, "xs", x_); 
    append_property<double_t>(d, "ys", y_);
    append_property<double_t>(d, "us", u_);

    append_property<double_t>(d, "tau_pluss", tau_plus_); 
    append_property<double_t>(d, "lambdas", lambda_); 
    append_property<double_t>(d, "alphas", alpha_); 
    append_property<double_t>(d, "mu_pluss", mu_plus_); 
    append_property<double_t>(d, "mu_minuss", mu_minus_);
    append_property<double_t>(d, "Wmaxs", Wmax_);
  }

} // of namespace nest
