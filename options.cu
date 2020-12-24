#include <string>
#include <fstream>

#include <iostream>

#include "options.hpp"

//__constant__ Options dopts;
//const Options& opts = options_details::options;

//Options opts;
void trim(std::string& s)
{
	const char* ws = " \t\n\r\f\v";
	s.erase(0, s.find_first_not_of(ws));
	auto last = s.find_last_not_of(ws);
	if(last != std::string::npos) 
		s.erase(s.find_last_not_of(ws) + 1);
}

void parse(std::string line, std::string& key, std::string& value)
{
	trim(line);
	auto index = line.find('=');
	if(!line.empty() && line[0] != '#' && index != std::string::npos)
	{
		key = line.substr(0, index);
		value = line.substr(index + 1);
	}
	trim(key);
	trim(value);
}

void load_options(Options& opts, const std::string& fname)
{
	std::ifstream f(fname.c_str());
	if(!f) throw std::runtime_error("Can't open " + fname);
	std::string line;
	while(std::getline(f, line))
	{
		std::string key, value;
		parse(line, key, value);
		if(key == "nx") opts.nx = std::stoi(value);
		else if(key == "ny") opts.ny = std::stoi(value);
		else if(key == "nz") opts.nz = std::stoi(value);
		else if(key == "npx") opts.npx = std::stoi(value);
	}
	opts.nyz = opts.ny*opts.nz;
	opts.nxyz = opts.nx*opts.nyz;
//	options_details::options.nx = 16;

//	cudaMemcpyToSymbol(dopts, &op, sizeof(Options));
}

void print_options(const Options& opts) 
{
	using std::cout;
	using std::endl;
	cout << "N_X: 	" << opts.nx << endl;
	cout << "N_Y: 	" << opts.ny << endl;
	cout << "N_Z: 	" << opts.nz << endl;
	cout << "N_YZ: 	" << opts.nyz << endl;
	cout << "N_XYZ: 	" << opts.nxyz << endl;
	cout << "N_PX: 	" << opts.npx << endl;
}

