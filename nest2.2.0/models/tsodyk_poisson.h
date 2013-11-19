/*
 *  tsodyk_poisson.h
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

#ifndef TSODYK_POISSON_H
#define TSODYK_POISSON_H


#include "connection_het_wd.h"
#include "archiving_node.h"
#include "generic_connector.h"
#include <cmath>

namespace nest
{
  class TsodykPoisson : public ConnectionHetWD
  {

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  TsodykPoisson();

  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  TsodykPoisson(const TsodykPoisson &);

  /**
   * Default Destructor.
   */
  ~TsodykPoisson() {}


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

  //tsodyks param
  double_t U_;
  double_t D_;
  double_t F_;
  double_t u_;
  double_t R;
  //trosdyks para
  // data members of each connection

  };





/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void TsodykPoisson::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{

  double_t t_spike = e.get_stamp().get_ms();
  double_t dendritic_delay = Time(Time::step(delay_)).get_ms();
  
  double_t u_old=u_;
  double_t R_old=R;
  u_=U_+u_old*(1.0-U_)*std::exp(-(t_spike-t_lastspike)/F_);
  R=1.0+(R_old-u_old*R_old-1.0)*std::exp(-(t_spike-t_lastspike)/D_);
  
  //tsodyk
  e.set_receiver(*target_);
  e.set_weight(weight_*u_*R);
  //e.set_weight(weight_);
  e.set_delay(delay_);
  e.set_rport(rport_);
  e();
}

} // of namespace nest

#endif // of #ifndef STDP_CONNECTION_WITHOUT_MIN_H
