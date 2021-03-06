/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSLOCALIZATION_H_
#define OOPS_INTERFACE_OBSLOCALIZATION_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/ObsLocalizationBase.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Encapsulates the observation-space localization
/// Note: to see methods that need to be implemented in the ObsLocalization implementation,
/// see ObsLocalizationBase class
template <typename MODEL, typename OBS>
class ObsLocalization : public util::Printable,
                        private boost::noncopyable {
  typedef ObsLocalizationBase<MODEL, OBS>  ObsLocBase_;
  typedef GeometryIterator<MODEL>  GeometryIterator_;
  typedef ObsSpace<OBS>            ObsSpace_;
  typedef ObsDataVector<OBS, int>  ObsDataVector_;
  typedef ObsVector<OBS>           ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsLocalization";}

  ObsLocalization(const eckit::Configuration &, const ObsSpace_ &);
  ~ObsLocalization();

  /// compute obs-space localization: fill \p obsvector with observation-space
  /// localization values between observations and \p point in model-space, and
  /// fill \p outside with flags on whether obs is local or not (1: outside of
  /// localization, 0: inside of localization, local)
  void computeLocalization(const GeometryIterator_ & point,
                           ObsDataVector_ & outside, ObsVector_ & obsvector) const override;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<ObsLocBase_> obsloc_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
ObsLocalization<MODEL, OBS>::ObsLocalization(const eckit::Configuration & conf,
                                             const ObsSpace_ & obspace)
  : obsloc_()
{
  Log::trace() << "ObsLocalization<MODEL, OBS>::ObsLocalization starting" << std::endl;
  util::Timer timer(classname(), "ObsLocalization");
  obsloc_.reset(ObsLocalizationFactory<MODEL, OBS>::create(conf, obspace));
  Log::trace() << "ObsLocalization<MODEL, OBS>::ObsLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
ObsLocalization<MODEL, OBS>::~ObsLocalization() {
  Log::trace() << "ObsLocalization<MODEL, OBS>::~ObsLocalization starting" << std::endl;
  util::Timer timer(classname(), "~ObsLocalization");
  obsloc_.reset();
  Log::trace() << "ObsLocalization<MODEL, OBS>::~ObsLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void ObsLocalization<MODEL, OBS>::computeLocalization(const GeometryIterator_ & p,
                                  ObsDataVector_ & local, ObsVector_ & obsvector) const {
  Log::trace() << "ObsLocalization<MODEL, OBS>:: computeLocalization starting" << std::endl;
  util::Timer timer(classname(), "computeLocalization");
  obsloc_->computeLocalization(p, local, obsvector);
  Log::trace() << "ObsLocalization<MODEL, OBS>:: computeLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void ObsLocalization<MODEL, OBS>::print(std::ostream & os) const {
  os << *obsloc_;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSLOCALIZATION_H_
