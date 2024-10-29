#include <gmp.h>
#include <iostream>
#include <cstdio>
#include <limits>
#include <csignal>
#include <boost/chrono.hpp>

static volatile int interrupted{0};

void signal_handler(int signal)
{
  interrupted = signal;
}

void mersenne_init(mpz_t& n, unsigned long int power)
{
    mpz_init(n);
    mpz_ui_pow_ui(n,2,power);
    mpz_sub_ui(n,n,1);
}

unsigned int collatz(mpz_t& n, unsigned long int c)
{
    while (!interrupted && mpz_cmp_ui(n,1)>0)
    {
        if (mpz_odd_p(n))
        {
            mpz_mul_ui(n,n,3);
            mpz_add_ui(n,n,1);
            c+=1;
        }
        mpz_tdiv_q_2exp(n,n,1);
        c+=1;
    }
    return c;
}

bool get_positive_number(char *power, unsigned long int& res)
{
  char *pos;
  long long int val = strtoll(power, &pos, 10);
  if (pos == power)
  {
    std::cerr << "Needs a positive number as parameter" << std::endl;
    return false;
  }
  if (val<=0)
  {
    std::cerr << "Needs a positive number as parameter" << std::endl;
    return false;
  }

  if (val>std::numeric_limits<unsigned int>::max())
  {
    std::cerr << "Number too big" << std::endl;
    return false;
  }
  res = val;
  return true;
}

bool save_cache(const std::string& cache, const mpz_t &n, double count, unsigned int steps)
{
    FILE* fp = std::fopen(cache.c_str(), "w");
    if (!fp) return false;
    if (std::fwrite(static_cast<void*>(&count),sizeof(count),1,fp) != 1) goto error_handler;
    if (std::fwrite(static_cast<void*>(&steps),sizeof(steps),1,fp) != 1) goto error_handler;
    if (!mpz_out_raw(fp, n)) goto error_handler;
    std::fclose(fp);
    return true;
error_handler:
    std::cerr << "Error saving cache" << std::endl;
    std::fclose(fp);
    std::remove(cache.c_str());
    return false;
}

bool load_cache(const std::string& cache, mpz_t &n, double &count, unsigned int &steps)
{
    double tmp_count;
    unsigned int tmp_steps;
    FILE* fp = std::fopen(cache.c_str(), "r");
    if (!fp) return false;
    if (std::fread(static_cast<void*>(&tmp_count),sizeof(tmp_count),1,fp) != 1) return false;
    if (std::fread(static_cast<void*>(&tmp_steps),sizeof(tmp_steps),1,fp) != 1) return false;
    mpz_init(n);
    if (!mpz_inp_raw(n, fp)) return false;
    count = tmp_count;
    steps = tmp_steps;
    return true;
}

bool load_file(const char* filename, mpz_t &n)
{
    FILE* fp = std::fopen(filename, "r");
    if (!fp) return false;
    mpz_init(n);
    auto ret = mpz_inp_str(n, fp, 0);
    std::fclose(fp);
    return ret;
}

int main(int argc, char ** argv)
{
    // Install a signal handler
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);
    std::signal(SIGHUP, signal_handler);
    std::signal(SIGALRM, signal_handler);
    std::string cache;
    std::string sw;
    std::string type;
    unsigned int steps{};
    if (argc == 3)
    {
        cache = argv[2];
        sw = argv[1];
        cache += sw;
        cache += ".cache";
    }
    else
    {
        std::cout << "Usage: collatz_gmp {-f filename} | { -m power } | {-n value }\n";
        std::cout << "       power_of_2 will create a mersenne number with the value\n";
        std::cout << "         2**power - 1\n";
        std::cout << "       filename should contain an arbitrarily large integer\n";
        std::cout << "       value is an arbitrarily large integer"  << std::endl;
        exit(1);
    }
    boost::chrono::system_clock::time_point start;
    double already_run{0};
    mpz_t n;
    if (!load_cache(cache, n, already_run, steps))
    {
        if (sw == "-m")
        {
            unsigned long int power = 0;
            if (!get_positive_number(argv[2],power))
            {
                exit(1);
            }
            type = "mersenne";
            mersenne_init(n,power);
        }
        else if (sw == "-f")
        {
            if (!load_file(argv[2], n))
            {
                std::cerr << "Could not load file '" << argv[2] << std::endl;
                exit(1);
            }
            type = "file";
        }
        else if (sw == "-n")
        {
            if (mpz_init_set_str(n,argv[2],0) == -1)
            {
                std::cerr << "Invalid number" << std::endl;
                exit(1);
            }
            if (mpz_cmp_ui(n,1)<0)
            {
                std::cerr << "Must be a positive integer" << std::endl;
                exit(1);
            }
            type = "number";
        }
        else
        {
            std::cerr << "Unknown option '" << argv[1] << "'" << std::endl;
            exit(1);
        }
    }
    boost::chrono::duration<double,boost::ratio<1>> dur(already_run);
calculate:
    alarm(60);
    start = boost::chrono::system_clock::now();
    steps = collatz(n, steps);
    dur += boost::chrono::system_clock::now() - start;
    if (interrupted)
    {
        if (interrupted != SIGALRM)
        {
            std::cerr << "\ninterrupted, saving cache file: " << cache << std::endl;
        }
        save_cache(cache,n, dur.count(), steps);
        if (interrupted == SIGHUP || interrupted == SIGALRM)
        {
            interrupted = 0;
            goto calculate;
        }
    }
    else
    {
        auto sec = dur.count();
        int min = sec/60;
        sec = sec-(min*60);
        int hour = min/60;
        min = min - (hour*60);
        int day = hour/24;
        hour = hour - (day*24);
        std::cout << type << "," << argv[2] << "," << steps <<"," << "\"" << day << "d " << hour << ":" << min << ":" << sec << "\"," << dur.count() << std::endl;
        std::remove(cache.c_str());
    }
}
