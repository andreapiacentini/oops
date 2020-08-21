/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/DateTime.h"

#include <stdint.h>

#include <cmath>
#include <iomanip>
#include <istream>
#include <limits>
#include <sstream>
#include <string>

#include "oops/util/abor1_cpp.h"
#include "oops/util/dateFunctions.h"
#include "oops/util/Duration.h"

namespace df = util::datefunctions;

namespace util {

// -----------------------------------------------------------------------------

DateTime::DateTime() : date_(0), time_(0) {
}

// -----------------------------------------------------------------------------
DateTime::~DateTime() {
}

// -----------------------------------------------------------------------------

DateTime::DateTime(const std::string & str) {
  this->set(str);
}

// -----------------------------------------------------------------------------

DateTime::DateTime(const int &YYYY, const int &MM, const int &DD,
                   const int &hh, const int &mm, const int &ss)
    : date_(df::dateToJulian(YYYY, MM, DD)),
      time_(df::hmsToSeconds(hh, mm, ss)) {
}

// -----------------------------------------------------------------------------

void DateTime::set(const std::string & str) {
  int year;
  int month;
  int day;
  int hour;
  int minute;
  int second;
  df::stringToYYYYMMDDhhmmss(str, year, month, day, hour, minute, second);
  date_ = df::dateToJulian(year, month, day);
  time_ = df::hmsToSeconds(hour, minute, second);
}


// -----------------------------------------------------------------------------

std::string DateTime::toString() const {
  failIfUnset();

  int year;
  int month;
  int day;
  int hour;
  int minute;
  int second;

  df::julianToDate(date_, year, month, day);
  df::secondToHms(time_, hour, minute, second);

  std::ostringstream os;
  os << std::setfill('0');
  os << std::setw(4) << year;
  os.put('-');
  os << std::setw(2) << month;
  os.put('-');
  os << std::setw(2) << day;
  os.put('T');
  os << std::setw(2) << hour;
  os.put(':');
  os << std::setw(2) << minute;
  os.put(':');
  os << std::setw(2) << second;
  os.put('Z');
  return os.str();
}

// -----------------------------------------------------------------------------

void DateTime::toYYYYMMDDhhmmss(int & year, int & month, int & day,
                                int & hour, int & minute, int & second) const {
  df::julianToDate(date_, year, month, day);
  df::secondToHms(time_, hour, minute, second);
}

// -----------------------------------------------------------------------------

void DateTime::addSeconds(const int64_t &seconds) {
  failIfUnset();

  int64_t days = 0;
  int64_t new_time = 0;
  int64_t secondsPerDay = 86400;

  new_time = seconds;
  days = new_time / secondsPerDay;

  date_ += days;
  new_time %= secondsPerDay;

  time_ += static_cast<int>(new_time);

  if (time_ < 0) {
    date_--;
    time_ += secondsPerDay;
  } else if (time_ >= secondsPerDay) {
    date_++;
    time_ -= secondsPerDay;
  }
}

// -----------------------------------------------------------------------------

DateTime& DateTime::operator+=(const util::Duration & s) {
  if (s.toSeconds()<std::numeric_limits<int64_t>::min() ||
      s.toSeconds()>std::numeric_limits<int64_t>::max()) ABORT("out of range");
  failIfUnset();
  int64_t  secs = s.toSeconds();
  addSeconds(secs);
  return *this;
}

// -----------------------------------------------------------------------------

DateTime& DateTime::operator-=(const util::Duration & s) {
  failIfUnset();
  util::Duration negs = s;
  negs.negate();
  this->operator+=(negs);
  return *this;
}

// -----------------------------------------------------------------------------

const DateTime DateTime::operator+(const util::Duration & s) const {
  failIfUnset();
  DateTime result = *this;  // Make a copy of myself.
  result += s;         // Use += to add s to the copy.
  return result;
}

// -----------------------------------------------------------------------------

const DateTime DateTime::operator-(const util::Duration & s) const {
  failIfUnset();
  DateTime result = *this;  // Make a copy of myself.
  result -= s;         // Use -= to subtract s from the copy
  return result;
}

// -----------------------------------------------------------------------------

const util::Duration DateTime::operator-(const DateTime& other) const {
  failIfUnset();
  other.failIfUnset();
  int secondsPerDay = 86400;
  return util::Duration(static_cast<int64_t>(
            ((date_ - other.date_)*secondsPerDay + (time_ -other.time_))) );
}

// -----------------------------------------------------------------------------

bool DateTime::operator==(const DateTime& other) const {
  failIfUnset();
  other.failIfUnset();
  return (date_ == other.date_) && (time_ == other.time_);
}

// -----------------------------------------------------------------------------

bool DateTime::operator!=(const DateTime& other) const {
  failIfUnset();
  other.failIfUnset();
  return !(*this == other);
}

// -----------------------------------------------------------------------------

bool DateTime::operator<(const DateTime& other) const {
  failIfUnset();
  other.failIfUnset();
  return (*this - other) < util::Duration(0);
}

// -----------------------------------------------------------------------------

bool DateTime::operator<=(const DateTime& other) const {
  failIfUnset();
  other.failIfUnset();
  return (*this - other) <= util::Duration(0);
}

// -----------------------------------------------------------------------------

bool DateTime::operator>(const DateTime& other) const {
  failIfUnset();
  other.failIfUnset();
  return (*this - other) > util::Duration(0);
}

// -----------------------------------------------------------------------------

bool DateTime::operator>=(const DateTime& other) const {
  failIfUnset();
  other.failIfUnset();
  return (*this - other) >= util::Duration(0);
}

// -----------------------------------------------------------------------------

void DateTime::failIfUnset() const {
  if (date_ == 0) {
    ABORT("DateTime was default-constructed and never set");
  }
}

// -----------------------------------------------------------------------------

std::size_t DateTime::timestamp() const {
  return date_ * 25 + time_;
}

// -----------------------------------------------------------------------------
size_t DateTime::serialSize() const {
  return 2;
}
// -----------------------------------------------------------------------------
void DateTime::serialize(std::vector<double> & vect) const {
  vect.push_back(static_cast<double>(date_));
  vect.push_back(static_cast<double>(time_));
}

// -----------------------------------------------------------------------------
void DateTime::deserialize(const std::vector<double> & vect, size_t & current) {
  date_ = std::lround(vect.at(current));
  time_ = std::lround(vect.at(current+1));
  current += 2;
}

// -----------------------------------------------------------------------------

/// To use a date as key in a boost::map
std::size_t hash_value(const util::DateTime& d) {
  return d.timestamp();
}

// -----------------------------------------------------------------------------

}  // namespace util

