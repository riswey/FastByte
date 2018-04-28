#pragma once

//FIX to_string Bug ///////////////////////////////
namespace patch
{
	template < typename T > std::string to_string(const T& n)
	{
		std::ostringstream stm;
		stm << n;
		return stm.str();
	}
}

#define forl(i,s,f) for(int i=s;i<f;i++)

