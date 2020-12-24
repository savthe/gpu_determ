#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <string>

struct Options
{
	int nx;
	int ny;
	int nz;
	int nyz;
	int nxyz;
	int npx;
};

void load_options(Options&, const std::string&);
void print_options(const Options&);

/*
namespace options_details
{
	__managed__ Options options;
}

__managed__ const Options* opts = &options_details::options;
*/
/*
extern Options opts;
static __constant__ Options dopts;

void load_options(std::string fname = "gpu.conf");
void print_options();
*/

/*
class Options
{
public:
	struct Storage
	{
		int nx;

	};

	static Options& instance()
	{
		static Options opts;
		return opts;
	}

	static const Storage storage() 
	{
		return storage_;
	}

private:
	Options() 
	{
		storage_.nx = 16;
	}
	static Storage storage_;
};

extern const Options::Storage& opts;
*/

#endif // OPTIONS_HPP

