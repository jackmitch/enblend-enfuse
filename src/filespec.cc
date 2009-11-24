/*
 * Copyright (C) 2009 Christoph L. Spiel
 *
 * This file is part of Enblend.
 *
 * Enblend is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Enblend is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Enblend; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef _WIN32
#include <io.h>
#else
#include <glob.h>
#include <unistd.h>
#endif

#include <cerrno>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>

#include <boost/assign/list_of.hpp>

#include "vigra/imageinfo.hxx"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "global.h"
#include "error_message.h"
#include "filenameparse.h"
#include "filespec.h"


// How many lines at the beginning of a response file we check to
// guess that it actually is a response file.
#define RESPONSE_FILE_LINES_TO_CHECK 20U


// Separator of key and value in a syntactic comment
#define KEY_VALUE_SEPARATOR_CHAR ':'


typedef std::pair<std::string, std::string> key_value_pair;


extern const std::string command;
extern int Verbose;


namespace enblend
{

#ifdef _WIN32
/** Add all files which match filename to filelist.  filename can
 *  contain the wildcards '*' and '?' in the leaf position, but not in
 *  the path. */
void
glob_filename_win32(FileNameList& filelist, const std::string& filename)
{
    // There has got to be an easier way...
    char drive[_MAX_DRIVE];
    char dir[_MAX_DIR];
    char fname[_MAX_FNAME];
    char ext[_MAX_EXT];
    char newFile[_MAX_PATH];

    _splitpath(filename.c_str(), drive, dir, NULL, NULL);

    struct _finddata_t finddata;
    intptr_t findhandle = _findfirst(filename.c_str(), &finddata);
    if (findhandle != -1)
    {
        do
        {
            _splitpath(finddata.name, NULL, NULL, fname, ext);
            _makepath(newFile, drive, dir, fname, ext);
            filelist.push_back(std::string(newFile));
        }
        while (_findnext(findhandle, &finddata) == 0);
        _findclose(findhandle);
    }
}
#endif


// ANTICIPATED CHANGE: We already have the same function in
// "common.h".
static bool
can_open_file(const std::string& aFilename)
{
    errno = 0;
    std::ifstream file(aFilename.c_str());
    if (!file)
    {
        std::cerr << command <<
            ": failed to open \"" << aFilename << "\": " <<
            errorMessage(errno) << "\n";
        return false;
    }
    else
    {
        errno = 0;
        file.close();
        if (file.fail())
        {
            std::cerr << command <<
                ": info: problems when closing \"" << aFilename << "\": " <<
                errorMessage(errno) << "\n";
        }
        return true;
    }
}


// ANTICIPATED CHANGE: We already have the same function in
// "common.h".  There it is called toLower().  However, we cannot
// include this file here.  The solution is to factor out string
// functions into their own module.
std::string
lower_case(const std::string& aString)
{
    std::string result(aString);
    std::transform(aString.begin(), aString.end(), result.begin(), tolower);
    return result;
}


/** Convert a_string in a string that is printable.  This is it only
 *  contains printable characters. */
std::string
printable_string(const std::string& a_string)
{
    std::string result;

    for (std::string::const_iterator c = a_string.begin();
         c != a_string.end();
         ++c)
    {
        if (isprint(*c))
        {
            result.push_back(*c);
        }
        else
        {
            std::ostringstream oss;
            oss <<
                "\\x" <<
                std::hex << std::setw(2) << std::setfill('0') <<
                (static_cast<int>(*c) & 0xff);
            result.append(oss.str());
        }
    }

    return result;
}


/** Answer line will leading and trailing whitespace removed.
 *  Normalize blank lines and lines starting with a comment character
 *  to an empty string.  Remove all whitespace between the
 *  response-file character and the following filename. */
std::string
normalize_response_file_line(const std::string& line)
{
    std::string result;
    std::string::size_type begin = line.find_first_not_of(WHITESPACE_CHARS);

    if (begin != std::string::npos && line[begin] != RESPONSE_FILE_COMMENT_CHAR)
    {
        const std::string::size_type end = line.find_last_not_of(WHITESPACE_CHARS);

        if (line[begin] == RESPONSE_FILE_PREFIX_CHAR)
        {
            begin = line.find_first_not_of(WHITESPACE_CHARS, begin + 1);
            result += RESPONSE_FILE_PREFIX_CHAR;
        }
        if (begin != std::string::npos)
        {
            result += line.substr(begin, end - begin + 1);
        }
    }

    return result;
}


static const key_value_pair empty_pair(std::make_pair(std::string(), std::string()));


/** Answer the key-value pair in line or a pair of null strings, if
 * line does not define a syntactic comment.
 *
 * A syntactic comment line matches
 *         ^[ \t]* # [ \t]* ([A-Za-z-]*) [ \t]* : [ \t]* ([^ \t]*) [ \t]*$
 * where the first captured string is the key and second one is the
 * value.
 */
key_value_pair
get_syntactic_comment(const std::string& line)
{
    std::string::size_type begin = line.find_first_not_of(WHITESPACE_CHARS);

    if (begin != std::string::npos && line[begin] == RESPONSE_FILE_COMMENT_CHAR)
    {
        begin = line.find_first_not_of(WHITESPACE_CHARS, begin + 1);
        if (begin != std::string::npos)
        {
            std::string::size_type end = begin + 1;
            while (end < line.size() && (isalpha(line[end]) || line[end] == '-'))
            {
                ++end;
            }
            const std::string key = line.substr(begin, end - begin);

            begin = line.find_first_not_of(WHITESPACE_CHARS, end);
            if (begin != std::string::npos && line[begin] == KEY_VALUE_SEPARATOR_CHAR)
            {
                begin = line.find_first_not_of(WHITESPACE_CHARS, begin + 1);
                if (begin == std::string::npos)
                {
                    return std::make_pair(key, std::string());
                }
                else
                {
                    end = line.find_last_not_of(WHITESPACE_CHARS);
                    return std::make_pair(key, line.substr(begin, end - begin + 1));
                }
            }
            else
            {
                return empty_pair;
            }
        }
        else
        {
            return empty_pair;
        }
    }
    else
    {
        return empty_pair;
    }
}


void
unroll_trace(const FilePositionTrace& trace)
{
    for (FilePositionTrace::const_iterator t = trace.begin(); t != trace.end(); ++t)
    {
        std::cerr << command <<
            ": info:    included from \"" << t->first << "\" line " << t->second << '\n';
    }
}


class AbstractGlobbingAlgorithm
{
public:
    virtual FileNameList do_glob(const std::string& a_filename,
                                 const FilePositionTrace& trace) = 0;
    virtual ~AbstractGlobbingAlgorithm() {}
    virtual std::string description() const = 0;
};


class LiteralGlobbingAlgorithm: public AbstractGlobbingAlgorithm
{
public:
    FileNameList do_glob(const std::string& a_filename,
                         const FilePositionTrace& trace)
    {
        FileNameList result;
        result.push_back(a_filename);
        return result;
    }

    std::string description() const
    {
        return "Do not glob.  Treat filenames as literals.";
    }
};


#ifdef _WIN32

class WildcardGlobbingAlgorithm: public AbstractGlobbingAlgorithm
{
public:
    FileNameList do_glob(const std::string& a_filespec,
                         const FilePositionTrace& trace)
    {
        FileNameList result;
        glob_filename_win32(result, a_filespec);
        return result;
    }

    std::string description() const
    {
        return "Glob only filenames, but not paths.";
    }
};

#else

class WildcardGlobbingAlgorithm: public AbstractGlobbingAlgorithm
{
public:
    WildcardGlobbingAlgorithm() : result_vector_(new glob_t) {}

    ~WildcardGlobbingAlgorithm() {delete result_vector_;}

    void run_glob(const std::string& a_filespec,
                  const FilePositionTrace& trace,
                  int flags)
    {
        errno = 0;
        if (glob(a_filespec.c_str(),
                 flags,
                 NULL, // (*errfunc)(const char* filename, int error_code)
                 result_vector_) != 0)
        {
            std::cerr << command <<
                ": warning: globbing \"" << a_filespec << "\" failed: " <<
                errorMessage(errno) << "\n";
            unroll_trace(trace);
        }
    }

    void convert_to_list()
    {
        result_.clear();

        for (char** path = result_vector_->gl_pathv; *path != NULL; ++path)
        {
            result_.push_back(*path);
        }
    }

    FileNameList do_glob(const std::string& a_filespec,
                         const FilePositionTrace& trace)
    {
        run_glob(a_filespec, trace, GLOB_ERR);
        convert_to_list();
        globfree(result_vector_);

        return result_;
    }

    std::string description() const
    {
        return "Glob with wildcards '?', '*', '[', and ']'.  See glob(7).";
    }

protected:
    glob_t* result_vector_;
    FileNameList result_;
};


class ShellGlobbingAlgorithm : public WildcardGlobbingAlgorithm
{
public:
    FileNameList do_glob(const std::string& a_filespec,
                         const FilePositionTrace& trace)
    {
        run_glob(a_filespec, trace, GLOB_ERR | GLOB_BRACE | GLOB_TILDE);
        convert_to_list();
        globfree(result_vector_);

        return result_;
    }

    std::string description() const
    {
        return "Glob like UN*X shells do.  Like \"wildcard\" plus '{', '}', and '~'.  See glob(7).";
    }
};

#endif // _WIN32


#define MAKE_ALGORITHM(m_algorithm_pointer) \
    std::make_pair(false, static_cast<AbstractGlobbingAlgorithm*>(m_algorithm_pointer))

#define MAKE_ALIAS(m_algorithm_pointer) \
    std::make_pair(true, m_algorithm_pointer)


class Globbing
{
    // Map from the names of the algorithms or aliases to the
    // algorithm itself.  The boolean flag indicates whether we point
    // to the algorithm proper or the name is just an alias to an
    // existing algorithm.
    typedef std::map<std::string, std::pair<bool, AbstractGlobbingAlgorithm*> > algorithm_map;

public:
    Globbing() : algorithm_name_("literal"), algorithm_(NULL)
    {
        installed_algorithms_ =
            boost::assign::map_list_of
            ("literal", MAKE_ALGORITHM(new LiteralGlobbingAlgorithm))
            ("wildcard", MAKE_ALGORITHM(new WildcardGlobbingAlgorithm))
#ifndef _WIN32
            ("shell", MAKE_ALGORITHM(new ShellGlobbingAlgorithm))
#endif
            ;

        setup_alias("literal", "none");
#ifndef _WIN32
        setup_alias("shell", "sh");
#endif
    }

    ~Globbing()
    {
        for (algorithm_map::const_iterator i = installed_algorithms_.begin();
             i != installed_algorithms_.end();
             ++i)
        {
            if (!i->second.first) // not an alias
            {
                delete i->second.second;
            }
        }

        algorithm_ = NULL;
    }

    const std::string& get_algorithm() const
    {
        return algorithm_name_;
    }

    bool set_algorithm(const std::string& an_algorithm_name)
    {
        algorithm_map::const_iterator a = installed_algorithms_.find(an_algorithm_name);
        if (a != installed_algorithms_.end())
        {
            algorithm_name_ = a->first;
            algorithm_ = a->second.second;
            return true;
        }
        else
        {
            return false;
        }
    }

    FileNameList expand(const std::string& a_filename, const FilePositionTrace& trace)
    {
        if (algorithm_ == NULL)
        {
            set_algorithm(algorithm_name_);
        }
        return algorithm_->do_glob(a_filename, trace);
    }

    algorithm_list get_known_algorithms() const
    {
        algorithm_list result;
        for (algorithm_map::const_iterator i = installed_algorithms_.begin();
             i != installed_algorithms_.end();
             ++i)
        {
            std::string description = i->second.second->description();
            if (i->second.first)
            {
                description += " (alias)";
            }
            result.push_back(make_pair(i->first, description));
        }
        return result;
    }

    void setup_alias(const std::string& an_algorithm_name, const std::string& an_alias)
    {
        algorithm_map::const_iterator a = installed_algorithms_.find(an_algorithm_name);
        installed_algorithms_[an_alias] = MAKE_ALIAS(a->second.second);
    }

private:
    std::string algorithm_name_;
    AbstractGlobbingAlgorithm* algorithm_;

    algorithm_map installed_algorithms_;
};


algorithm_list
known_globbing_algorithms()
{
    Globbing glob;
    return glob.get_known_algorithms();
}


//< src::RESPONSE_FILE_MAX_NESTING_LEVEL 63
#define RESPONSE_FILE_MAX_NESTING_LEVEL 63U


void
unfold_filename_iter(TraceableFileNameList& result,
                     unsigned nesting_level, const std::string& current_directory, FilePositionTrace& trace,
                     const std::string& filename)
{
    // Checking the nesting_lavel acts as an emergency break in the
    // case our usual recusion detection fails.
    if (nesting_level > RESPONSE_FILE_MAX_NESTING_LEVEL)
    {
        std::cerr <<
            command << ": warning: excessive nesting of " << nesting_level <<
            " levels of response files;\n" <<
            command << ": warning:     possible infinite recursion in \"" <<
            printable_string(filename) << "\"\n";
        unroll_trace(trace);
        return;
    }

    if (filename.empty())
    {
        std::cerr << command << ": info: empty filename\n";
        unroll_trace(trace);
        return;
    }
    else if (filename[0] == RESPONSE_FILE_PREFIX_CHAR)
    {
        const std::string response_filename(filename, 1); // filename alone
#ifdef DEBUG_FILESPEC
        std::cout <<
            "+ unfold_filename_iter: concatPath(current_directory, response_filename) = <" <<
            concatPath(current_directory, response_filename) << ">\n";
#endif
        const std::string response_filepath = // filename in the right path
            enblend::canonicalizePath(enblend::isRelativePath(response_filename) ?
                                      concatPath(current_directory, response_filename) :
                                      response_filename,
                                      false);

        for (FilePositionTrace::const_iterator t = trace.begin(); t != trace.end(); ++t)
        {
            if (t->first == response_filepath)
            {
                std::cerr << command <<
                    ": warning: response file \"" << printable_string(response_filepath) <<
                    "\" recursively\n";
                unroll_trace(trace);
                return;
            }
        }

#ifdef DEBUG_FILESPEC
        std::cout <<
            "+ unfold_filename_iter: current_directory = " << current_directory << "\n" <<
            "+ unfold_filename_iter: response_filename = " << response_filename << "\n" <<
            "+ unfold_filename_iter: response_filepath = " << response_filepath << "\n";
#endif

        if (Verbose >= VERBOSE_RESPONSE_FILES)
        {
            std::cerr << command <<
                ": info: consulting response file \"" <<
                printable_string(response_filepath) << "\"\n";
        }

        if (!can_open_file(response_filepath))
        {
            std::cerr << command <<
                ": failed to open response file \"" <<
                printable_string(response_filepath) << "\": " <<
                errorMessage(errno) << "\n";
            unroll_trace(trace);
            exit(1);
        }

        if (!maybe_response_file(response_filepath))
        {
            std::cerr << command <<
                ": warning: malformed response file \"" <<
                printable_string(response_filepath) << "\"\n";
            unroll_trace(trace);
            return;
        }

        errno = 0;
        std::ifstream response_file(response_filepath.c_str());
        Globbing glob;
        unsigned line_number = 0U;

        while (true)
        {
            std::string buffer;

            errno = 0;
            getline(response_file, buffer);
            if (!response_file.good())
            {
                if (!buffer.empty())
                {
                    std::cerr <<
                        command << ": warning: ignoring last line of response file \"" <<
                        printable_string(response_filepath) << "\" line " << line_number + 1U << ",\n" <<
                        command << ": warning:     because it does not end with a newline\n";
                    unroll_trace(trace);
                }
                break;
            }
            ++line_number;

            const key_value_pair comment = get_syntactic_comment(buffer);
            if ((comment.first == "glob" ||
                 comment.first == "globbing" ||
                 comment.first == "filename-globbing") &&
                !comment.second.empty())
                // We silently ignore all other keys or empty
                // values with the right key.
            {
                if (!glob.set_algorithm(lower_case(comment.second)))
                {
                    std::cerr <<
                        command <<
                        ": warning: requested unknown globbing algorithm \"" <<
                        comment.second <<
                        "\" in response file \"" << printable_string(response_filepath) <<
                        "\" line " << line_number << "\n" <<
                        command <<
                        ": warning:     will stick with algorithm \"" <<
                        glob.get_algorithm() << "\"\n";
                    unroll_trace(trace);
                }
            }

            const std::string line(normalize_response_file_line(buffer));
            if (line.empty())
            {
                continue;
            }

            const std::string response_directory = extractDirname(response_filepath);
#ifdef DEBUG_FILESPEC
            std::cout << "+ unfold_filename_iter: response_directory = " << response_directory << "\n";
#endif
            TraceableFileNameList partial_result;
            trace.push_front(std::make_pair(response_filepath, line_number));
            unfold_filename_iter(partial_result,
                                 nesting_level + 1U, response_directory, trace,
                                 line);

            for (TraceableFileNameList::const_iterator p = partial_result.begin();
                 p != partial_result.end();
                 ++p)
            {
                const FileNameList expanded_partial_result = glob.expand(p->first, trace);
                for (FileNameList::const_iterator q = expanded_partial_result.begin();
                     q != expanded_partial_result.end();
                     ++q)
                {
                    if (enblend::isRelativePath(*q))
                    {
                        const std::string path =
                            enblend::canonicalizePath(concatPath(response_directory, *q), false);
#ifdef DEBUG_FILESPEC
                        std::cout <<
                            "+ unfold_filename_iter: relative path\n" <<
                            "+ unfold_filename_iter:     q = " << *q << "\n" <<
                            "+ unfold_filename_iter:     path = <" << path << ">\n";
#endif
                        result.insert(result.end(), std::make_pair(path, trace));
                    }
                    else
                    {
#ifdef DEBUG_FILESPEC
                        std::cout << "+ unfold_filename_iter: absolute path\n";
#endif
                        result.insert(result.end(), std::make_pair(*q, trace));
                    }
                }
            }

            trace.pop_front();
        }
        if (!response_file.eof())
        {
            std::cerr << command <<
                ": warning: filesystem signals problems in response file \"" <<
                printable_string(response_filepath) << "\" line " << line_number << ": " <<
                errorMessage(errno) << "\n";
            unroll_trace(trace);
        }
    }
    else
    {
        result.push_back(std::make_pair(filename, trace));
    }
}


void
unfold_filename(TraceableFileNameList& result, const std::string& filename)
{
#ifdef _WIN32
    // This stupid wannbe OS, lacking a shell that deserves the name,
    // passes unexpanded wildcards in filename.  We ameliorate the
    // situation by doing some of the globbing ourselves; see
    // glob_filename_win32().
    FileNameList initial_files;

    if (!filename.empty() && filename[0] == RESPONSE_FILE_PREFIX_CHAR)
    {
        initial_files.push_back(filename);
    }
    else
    {
        // Not a response file, so we expand wildcards.
        glob_filename_win32(initial_files, filename);
    }

    for (FileNameList::const_iterator i = initial_files.begin();
         i != initial_files.end();
         ++i)
    {
        FilePositionTrace trace;
        unfold_filename_iter(result, 0U, "", trace, *i);
    }
#else
    FilePositionTrace trace;
    unfold_filename_iter(result, 0U, "", trace, filename);
#endif
}


bool
is_known_extension_to_vigra(const std::string& filename)
{
    const std::string known_extensions = vigra::impexListExtensions();
    const std::string extension =
        lower_case(filename.substr(filename.rfind(".") + 1));

    return known_extensions.find(extension) != std::string::npos;
}


bool
maybe_response_file(const std::string& filename)
{
    std::ifstream file(filename.c_str());
    if (file)
    {
        unsigned line_number = 0U;
        unsigned score = 0U;

        while (line_number < RESPONSE_FILE_LINES_TO_CHECK)
        {
            std::string buffer;

            getline(file, buffer);
            if (!file.good())
            {
                break;
            }
            ++line_number;

            if (line_number == 1U)
            {
                const key_value_pair comment = get_syntactic_comment(buffer);
                if ((comment.first == "response-file" ||
                     comment.first == "enblend-response-file" ||
                     comment.first == "enfuse-response-file") &&
                    lower_case(comment.second) == "true")
                {
                    return true;
                }
            }

            const std::string line(normalize_response_file_line(buffer));
            if (line.empty())
            {
                ++score;
                continue;
            }

            if (line[0] == RESPONSE_FILE_PREFIX_CHAR ||
                is_known_extension_to_vigra(line))
            {
                score += 2U;
            }
        }

#ifdef DEBUG_FILESPEC
        std::cout <<
            "+ maybe_response_file: score = " << score <<
            ", lines = " << line_number << "\n";
#endif
        return score >= line_number;
    }
    else
    {
        return false;
    }
}

} // namespace enblend
