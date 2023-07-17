#ifndef DLOG_PROGRESS_H
#define DLOG_PROGRESS_H

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace dlog
{

struct progress
{
                        progress(size_t total):
                            current_(0), total_(total)      { show_progress(); }

    progress&           operator++()                        { current_++;       if (current_ * 100 / total_ > (current_ - 1) * 100 / total_) show_progress(); check_done(); return *this; }
    progress&           operator=(size_t cur)               { current_ = cur;   show_progress(); check_done(); return *this; }
    progress&           operator()(const std::string& s)    { message_ = s;     show_progress(); check_done(); return *this; }
    template<class T>
    progress&           operator()(const T& x)              { std::ostringstream oss; oss << x; return (*this)(oss.str()); }

    inline void         show_progress() const;
    void                check_done() const                  { if (current_ >= total_) std::cout << "\n" << std::flush; }

    private:
        size_t          current_, total_;
        std::string     message_;
};

}

void
dlog::progress::
show_progress() const
{
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * current_ / total_;
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << std::setw(3) << current_ * 100 / total_ << "%";
    if (!message_.empty())
        std::cout << "  (" << message_ << ")";
    std::cout << "\r";
    std::cout.flush();
}

#endif
