 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       David Salinas
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
  *
  *    This program is free software: you can redistribute it and/or modify
  *    it under the terms of the GNU General Public License as published by
  *    the Free Software Foundation, either version 3 of the License, or
  *    (at your option) any later version.
  *
  *    This program is distributed in the hope that it will be useful,
  *    but WITHOUT ANY WARRANTY; without even the implied warranty of
  *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *    GNU General Public License for more details.
  *
  *    You should have received a copy of the GNU General Public License
  *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  */

#ifndef GUDHI_CLOCK_H_
#define GUDHI_CLOCK_H_


#include <boost/date_time/posix_time/posix_time.hpp>

class Clock{

public:
	Clock():end_called(false){
		startTime = boost::posix_time::microsec_clock::local_time( );
	}

	Clock(const std::string& msg_){
		end_called = false;
		begin();
		msg = msg_;
	}


	void begin() const{
		end_called = false;
		startTime = boost::posix_time::microsec_clock::local_time( );
	}

	void end() const{
		end_called = true;
		endTime = boost::posix_time::microsec_clock::local_time( );
	}

	void print() const{
		std::cout << *this << std::endl;
	}

	friend std::ostream& operator<< (std::ostream& stream,const Clock& clock){
		if(!clock.end_called)
			clock.end();

		if(!clock.end_called) stream << "end not called";
		else{
			stream << clock.msg <<":"<<clock.num_seconds() <<"s";
		}
		return stream;

	}

	double num_seconds() const{
		if(!end_called) return -1;
		return (endTime-startTime).total_milliseconds()/1000.;
	}

private:
	mutable boost::posix_time::ptime startTime, endTime;
	mutable bool end_called;
	std::string msg;

};


#endif /* GUDHI_CLOCK_H_ */
