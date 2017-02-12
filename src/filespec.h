/*
 * Copyright (C) 2009-2017 Christoph L. Spiel
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
#ifndef FILESPEC_H_
#define FILESPEC_H_

#include <list>
#include <string>

#include "selector.h"


#define RESPONSE_FILE_PREFIX_CHAR '@'     //< response-file-prefix-char \char64
#define RESPONSE_FILE_COMMENT_CHAR '#'    //< response-file-comment-char \char35
#define WHITESPACE_CHARS "\n\r\t "


namespace enblend
{
    typedef std::list<std::string> FileNameList;

    typedef std::pair<std::string, unsigned> FilePosition; /** Filename, line number pairs */
    typedef std::list<FilePosition> FilePositionTrace;     /** Traceback to a file position */


    class TraceableFileName
    {
    public:
        TraceableFileName() = delete;
        TraceableFileName(const std::string& a_filename,
                          const FilePositionTrace& a_trace,
                          selector::Abstract* a_selector);
        virtual ~TraceableFileName();
        virtual TraceableFileName* clone() const;

        const std::string& filename() const;
        const FilePositionTrace& trace() const;
        void unroll_trace() const;
        selector::Abstract* selector() const;

    protected:
        selector::Abstract* selector_;
        const std::string filename_;
        const FilePositionTrace trace_;
    };


    class TraceableFileNameAndLayer : public TraceableFileName
    {
        typedef TraceableFileName super;

    public:
        TraceableFileNameAndLayer() = delete;
        TraceableFileNameAndLayer(const std::string& a_filename,
                                  const FilePositionTrace& a_trace,
                                  const std::string& a_layer_specification);
        TraceableFileNameAndLayer(const TraceableFileNameAndLayer&);
        TraceableFileNameAndLayer& operator=(const TraceableFileNameAndLayer&) = delete;
        ~TraceableFileNameAndLayer() override;
        TraceableFileNameAndLayer* clone() const override;

        std::string layer_spec() const;

    protected:
        const std::string layer_specification_;
    };


    typedef std::list<TraceableFileName*> TraceableFileNameList; /** */


    /** Recursively unfold filename which may be a literal name or a
     *  response file.  The result is a list of only literal
     *  filenames. */
    void unfold_filename(TraceableFileNameList& result, const std::string& filename);

    /** Answer whether we suspect filename is a response file. */
    bool maybe_response_file(const std::string& filename);


    /** List of algorithm name / algorithm description pairs */
    typedef std::list<std::pair<std::string, std::string> > algorithm_list;

    /** Answer a list of all globbing algorithms including aliases
     *  that are known, i.e., have been compiled in. */
    algorithm_list known_globbing_algorithms();
}

#endif // FILESPEC_H_


// Local Variables:
// mode: c++
// End:
