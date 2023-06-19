/**
 * Author: Dmitriy Morozov <dmitriy@mrzv.org>
 * The interface is heavily influenced by GetOptPP (https://code.google.com/p/getoptpp/).
 * The parsing logic is from ProgramOptions.hxx (https://github.com/Fytch/ProgramOptions.hxx).
 *
 * History:
 *  - 2015-06-01: added Traits<...>::type_string() for long, unsigned long
 *  - ...
 *  - 2018-04-27: replace parsing logic with the one from ProgramOptions.hxx to
 *                make the parser compliant with [GNU Program Argument Syntax
 *                Conventions](https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html)
 *  - 2018-05-11: add dashed_non_option(), to accept arguments that are negative numbers
 */

#ifndef OPTS_OPTS_H
#define OPTS_OPTS_H

#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <memory>
#include <cctype>
#include <functional>

namespace opts {

// Converters
template<class T>
struct Converter
{
                    Converter()                     {}
    static
    bool            convert(const std::string& val, T& res)
    {
        std::istringstream iss(val);
        iss >> res;
        return !iss.fail() && iss.eof();
    }
};

// Type
template<class T>
struct Traits
{
    static std::string  type_string()               { return "UNKNOWN TYPE"; }
};

template<>
struct Traits<int>
{
    static std::string  type_string()               { return "INT"; }
};

template<>
struct Traits<short int>
{
    static std::string  type_string()               { return "SHORT INT"; }
};

template<>
struct Traits<long>
{
    static std::string  type_string()               { return "LONG"; }
};

template<>
struct Traits<unsigned>
{
    static std::string  type_string()               { return "UNSIGNED INT"; }
};

template<>
struct Traits<short unsigned>
{
    static std::string  type_string()               { return "SHORT UNSIGNED INT"; }
};

template<>
struct Traits<unsigned long>
{
    static std::string  type_string()               { return "UNSIGNED LONG"; }
};

template<>
struct Traits<float>
{
    static std::string  type_string()               { return "FLOAT"; }
};

template<>
struct Traits<double>
{
    static std::string  type_string()               { return "DOUBLE"; }
};

template<>
struct Traits<std::string>
{
    static std::string  type_string()               { return "STRING"; }
};


struct BasicOption
{
    using IsShort = std::function<bool(char)>;

                    BasicOption(char        s_,
                                std::string l_,
                                std::string default_,
                                std::string type_,
                                std::string help_):
                        s(s_), l(l_), d(default_), t(type_), help(help_)                    {}
    virtual         ~BasicOption()                                                          {}

    int             long_size() const                           { return l.size() + 1 + t.size(); }

    void            output(std::ostream& out, int max_long) const
    {
        out << "   ";
        if (s)
            out << '-' << s << ", ";
        else
            out << "    ";

        out << "--" << l << ' ';

        if (!t.empty())
            out << t;

        for (int i = long_size(); i < max_long; ++i)
            out <<  ' ';

        out << "   " << help;

        if (!d.empty())
        {
            out << " [default: " << d << "]";
        }
        out << '\n';
    }

    virtual bool    flag() const                                { return false; }
    virtual bool    parse(int argc, char** argv, int& i, int j, IsShort is_short);
    virtual bool    set(std::string arg) =0;

    char            s;
    std::string     l;
    std::string     d;
    std::string     t;
    std::string     help;
};

// Option
template<class T>
struct OptionContainer: public BasicOption
{
                    OptionContainer(char               s_,
                                    const std::string& l_,
                                    T&                 var_,
                                    const std::string& help_,
                                    const std::string& type_ = Traits<T>::type_string()):
                        BasicOption(s_, l_, default_value(var_), type_, help_),
                        var(&var_)                  {}

    static
    std::string     default_value(const T& def)
    {
        std::ostringstream oss;
        oss << def;
        return oss.str();
    }

    bool            set(std::string s) override             { return Converter<T>::convert(s, *var); }

    T*  var;
};

template<>
struct OptionContainer<bool>: public BasicOption
{
                    OptionContainer(char               s_,
                                    const std::string& l_,
                                    bool&              var_,
                                    const std::string& help_):
                        BasicOption(s_, l_, "", "", help_),
                        var(&var_)                          { *var = false; }

    bool            parse(int, char**, int&, int, IsShort) override                 { *var = true; return true; }
    bool            set(std::string) override                                       { return true; }
    bool            flag() const override                                           { return true; }

    bool*  var;
};

template<class T>
struct OptionContainer< std::vector<T> >: public BasicOption
{
                    OptionContainer(char               s_,
                                    const std::string& l_,
                                    std::vector<T>&    var_,
                                    const std::string& help_,
                                    const std::string& type_ = "SEQUENCE"):
                        BasicOption(s_, l_, default_value(var_), type_, help_),
                        var(&var_), first(true)             { }

    static
    std::string     default_value(const std::vector<T>& def)
    {
        std::ostringstream oss;
        oss << "(";
        if (def.size())
            oss << def[0];
        for (size_t i = 1; i < def.size(); ++i)
            oss << ", " << def[i];
        oss << ")";
        return oss.str();
    }

    bool            set(std::string s) override
    {
        if (first)
        {
            var->clear();
            first = false;
        }

        T x;
        bool result = Converter<T>::convert(s,x);
        var->emplace_back(std::move(x));
        return result;
    }

    std::vector<T>* var;
    mutable bool    first;
};


template<class T>
std::unique_ptr<BasicOption>
Option(char s, const std::string& l, T& var, const std::string& help)       { return std::unique_ptr<BasicOption>{new OptionContainer<T>(s, l, var, help)}; }

template<class T>
std::unique_ptr<BasicOption>
Option(char s, const std::string& l, T& var,
       const std::string& type, const std::string& help)                    { return std::unique_ptr<BasicOption>{new OptionContainer<T>(s, l, var, help, type)}; }

template<class T>
std::unique_ptr<BasicOption>
Option(const std::string& l, T& var, const std::string& help)               { return std::unique_ptr<BasicOption>{new OptionContainer<T>(0, l, var, help)}; }

template<class T>
std::unique_ptr<BasicOption>
Option(const std::string& l, T& var,
       const std::string& type, const std::string& help)                    { return std::unique_ptr<BasicOption>{new OptionContainer<T>(0, l, var, help, type)}; }

// PosOption
template<class T>
struct PosOptionContainer
{
                PosOptionContainer(T& var_):
                    var(&var_)                                              {}

    bool        parse(std::list<std::string>& args) const
    {
        if (args.empty())
            return false;

        bool result = Converter<T>::convert(args.front(), *var);
        if (!result)
            std::cerr << "error: failed to parse " << args.front() << '\n';
        args.pop_front();
        return result;
    }

    T*          var;
};

template<class T>
PosOptionContainer<T>
PosOption(T& var)                                                           { return PosOptionContainer<T>(var); }


// Options
struct Options
{
            Options():
                failed(false)                       {}

    inline
    Options&    operator>>(std::unique_ptr<BasicOption> opt);
    template<class T>
    Options&    operator>>(const PosOptionContainer<T>& poc);

                operator bool()                     { return !failed; }


    friend
    std::ostream&
    operator<<(std::ostream& out, const Options& ops)
    {
        int max_long = 0;
        for (auto& cur : ops.options)
        {
            int cur_long = cur->long_size();
            if (cur_long > max_long)
                max_long = cur_long;
        }

        out << "Options:\n";
        for (auto& cur : ops.options)
            cur->output(out, max_long);

        return out;
    }

    bool            parse(int argc, char** argv);

    void            unrecognized_option(std::string arg) const
    {
        std::cerr << "error: unrecognized option " << arg << '\n';
    }

    static bool     dashed_non_option(char* arg, BasicOption::IsShort is_short)
    {
        return arg[ 0 ] == '-'
                && (std::isdigit(arg[ 1 ]) || arg[ 1 ] == '.')
                && !is_short(arg[ 1 ]);
    }

    private:
        std::list<std::string>                      args;
        std::list<std::unique_ptr<BasicOption>>     options;
        bool                                        failed;
};

bool
BasicOption::parse(int argc, char** argv, int& i, int j, IsShort is_short)
{
    char* argument;
    char* cur_arg = argv[i];
    // -v...
    if (argv[i][j] == '\0')
    {
        // -v data
        if (i + 1 < argc && (argv[i+1][0] != '-' || Options::dashed_non_option(argv[i+1], is_short)))
        {
            ++i;
            argument = argv[i];
        } else
        {
            std::cerr << "error: cannot find the argument; ignoring " << argv[i] << '\n';
            return false;
        }
    } else if (argv[i][j] == '=')
    {
        // -v=data
        argument = &argv[i][j+1];
    } else if( j == 2 ) { // only for short options
        // -vdata
        argument = &argv[i][j];
    } else
    {
        std::cerr << "error: unexpected character \'" << argv[i][j] << "\' ignoring " << argv[i] << '\n';
        return false;
    }
    bool result = set(argument);
    if (!result)
        std::cerr << "error: failed to parse " << argument << " in " << cur_arg << '\n';
    return result;
}

bool
Options::parse(int argc, char** argv)
{
    std::map<char, BasicOption*>                    short_opts;
    std::map<std::string, BasicOption*>             long_opts;

    for (auto& opt : options)
    {
        if (opt->s)
            short_opts[opt->s] = opt.get();

        long_opts[opt->l] = opt.get();
    }

    auto is_short = [&short_opts](char c) -> bool   { return short_opts.find(c) != short_opts.end(); };

    for (int i = 1; i < argc; ++i)
    {
        if( argv[ i ][ 0 ] == '\0' )
            continue;
        if( argv[ i ][ 0 ] != '-' || dashed_non_option(argv[i], is_short))
            args.push_back(argv[i]);
        else
        {
            // -...
            if( argv[ i ][ 1 ] == '\0' )
            {
                // -
                args.push_back(argv[i]);
            } else if( argv[ i ][ 1 ] == '-' )
            {
                if( argv[ i ][ 2 ] == '\0' )
                {
                    // --
                    while( ++i < argc )
                        args.push_back(argv[i]);
                } else {
                    // --...
                    char* first = &argv[ i ][ 2 ];
                    char* last = first;
                    for(; *last != '=' && *last != '\0'; ++last);
                    if (first == last)
                    {
                        failed = true;
                        unrecognized_option(argv[i]);
                    } else
                    {
                        auto opt_it = long_opts.find(std::string{first,last});
                        if (opt_it == long_opts.end())
                        {
                            failed = true;
                            unrecognized_option(argv[i]);
                        } else
                        {
                            failed |= !opt_it->second->parse(argc, argv, i, last - argv[i], is_short);
                        }
                    }
                }
            } else
            {
                // -f...
                auto opt_it = short_opts.find(argv[i][1]);
                if (opt_it == short_opts.end())
                {
                    failed = true;
                    unrecognized_option(argv[i]);
                } else if (opt_it->second->flag())
                {
                    opt_it->second->parse(argc, argv, i, 0, is_short);      // arguments are meaningless; just sets the flag

                    // -fgh
                    char c;
                    for(int j = 1; (c = argv[i][j]) != '\0'; ++j)
                    {
                        if (!std::isprint(c) || c == '-')
                        {
                            failed = true;
                            std::cerr << "error: invalid character\'" << c << " ignoring " << &argv[i][j] << '\n';
                            break;
                        }
                        opt_it = short_opts.find(c);
                        if (opt_it == short_opts.end())
                        {
                            failed = true;
                            unrecognized_option("-" + std::string(1, c));
                            continue;
                        }
                        if (!opt_it->second->flag())
                        {
                            failed = true;
                            std::cerr << "error: non-void options not allowed in option packs; ignoring " << c << '\n';
                            continue;
                        }
                        opt_it->second->parse(argc, argv, i, 0, is_short);     // arguments are meaningless; just sets the flag
                    }
                } else
                {
                    failed |= !opt_it->second->parse(argc, argv, i, 2, is_short);
                }
            }
        }
    }

    return !failed;
}

Options&
Options::operator>>(std::unique_ptr<BasicOption> opt)
{
    options.emplace_back(std::move(opt));
    return *this;
}

template<class T>
Options&
Options::operator>>(const PosOptionContainer<T>& poc)
{
    if (!failed)
        failed = !poc.parse(args);
    return *this;
}

}

#endif
