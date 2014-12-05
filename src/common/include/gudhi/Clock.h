/*
 * Clock.h
 *
 *  Created on: Jun 17, 2014
 *      Author: dsalinas
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
