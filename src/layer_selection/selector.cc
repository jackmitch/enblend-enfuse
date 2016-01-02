/*
 * Copyright (C) 2010-2016 Dr. Christoph L. Spiel
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


#include <cassert>
#include <cctype>               // std::isspace(), std::tolower()
#include <iostream>             // WHILE DEBUGGING
#include <memory>               // std::unique_ptr

#include <boost/tokenizer.hpp>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "selector.h"


#define NUMERIC_OPTION_DELIMITERS ";:/" // FIXME: Already defined on "common.h".
extern const std::string command;


template <class iterator, class unary_function>
static std::string
mapconcat(iterator a_begin, iterator an_end, unary_function a_function, const std::string& a_separator)
{
    if (a_begin == an_end)
    {
        return std::string();
    }
    else
    {
        std::string result(a_function(*a_begin));

        while (++a_begin != an_end)
        {
            result.append(a_separator).append(a_function(*a_begin));
        }

        return result;
    }
}


namespace selector
{
    std::string
    AllLayers::name() const
    {
        return "all-layers";
    }


    std::string
    AllLayers::description() const
    {
        return "select all layers in any image;";
    }


    layer_ordered_list_t
    AllLayers::viable_layers(const ImageListInformation* an_image_info, const std::string& a_filename)
    {
        const ImageInfo* image_info = an_image_info->image_info_on(a_filename);
        layer_ordered_list_t result;

        for (unsigned i = 1; i <= image_info->number_of_layers(); ++i)
        {
            result.push_back(i);
        }

        return result;
    }


    std::string
    FirstLayer::name() const
    {
        return "first-layer";
    }


    std::string
    FirstLayer::description() const
    {
        return "select only first layer in each (multi-)layer image;";
    }


    std::string
    LastLayer::name() const
    {
        return "last-layer";
    }


    std::string
    LastLayer::description() const
    {
        return "select only last layer in each (multi-)layer image;";
    }


    layer_ordered_list_t
    LastLayer::viable_layers(const ImageListInformation* an_image_info, const std::string& a_filename)
    {
        const ImageInfo* image_info = an_image_info->image_info_on(a_filename);

        return layer_ordered_list_t(1, image_info->number_of_layers());
    }


    bool
    LastLayer::select(const ImageListInformation* an_image_info,
                      const std::string& a_filename, unsigned a_layer_index)
    {
        const ImageInfo* image_info = an_image_info->image_info_on(a_filename);

        return a_layer_index == image_info->number_of_layers();
    }


    std::string
    LargestLayer::name() const
    {
        return "largest-layer";
    }


    std::string
    LargestLayer::description() const
    {
        return "select largest layer in each (multi-)layer image;";
    }


    static unsigned
    find_largest_layer(const ImageListInformation* an_image_info, const std::string& a_filename)
    {
        const ImageInfo* image_info = an_image_info->image_info_on(a_filename);
        const unsigned number_of_layers = image_info->number_of_layers();
        unsigned max_index = 0;

        int max_area = -1;
        for (unsigned i = 0; i != number_of_layers; ++i)
        {
            const LayerInfo* layer_info = an_image_info->layer_info_on(a_filename, i);
            const int area = layer_info->width * layer_info->height;
            if (area > max_area)
            {
                max_area = area;
                max_index = i;
            }
        }

        return max_index;
    }


    bool
    LargestLayer::select(const ImageListInformation* an_image_info,
                         const std::string& a_filename,
                         unsigned a_layer_index)
    {
        cache_t::const_iterator index = cache_.find(a_filename);
        unsigned max_index;

        if (index == cache_.end())
        {
            max_index = 1 + find_largest_layer(an_image_info, a_filename);
            cache_.insert(cache_t::value_type(a_filename, max_index));
        }
        else
        {
            max_index = index->second;
        }

        return a_layer_index == max_index;
    }


    layer_ordered_list_t
    LargestLayer::viable_layers(const ImageListInformation* an_image_info, const std::string& a_filename)
    {
        cache_t::const_iterator index = cache_.find(a_filename);
        unsigned max_index;

        if (index == cache_.end())
        {
            max_index = 1 + find_largest_layer(an_image_info, a_filename);
            cache_.insert(cache_t::value_type(a_filename, max_index));
        }
        else
        {
            max_index = index->second;
        }

        return layer_ordered_list_t(1, max_index);
    }


    std::string
    NoLayer::name() const
    {
        return "no-layer";
    }


    std::string
    NoLayer::description() const
    {
        return "do not select any layer from any image;";
    }


    // Implementation Note: We cannot use the new, unified
    // initialization syntax with curlies here, because it employs
    // copy-semantics.  Unique-ptr requires move-semantics, though.
    algorithm_list algorithms;
    namespace initialization
    {
        struct initialize_algorithms
        {
            initialize_algorithms()
            {
                algorithms.emplace_back(new AllLayers);
                algorithms.emplace_back(new FirstLayer);
                algorithms.emplace_back(new LastLayer);
                algorithms.emplace_back(new LargestLayer);
                algorithms.emplace_back(new NoLayer);
                // NOTE: `IndexedLayer' does _not_ belong here!
            }
        } algorithms_initializer;
    } // namespace initialization


    algorithm_list::const_iterator
    find_by_id(id_t an_id)
    {
        algorithm_list::const_iterator algorithm =
            std::find_if(algorithms.begin(),
                         algorithms.end(),
                         [&an_id](const algorithm_list::value_type& x) {return x->id() == an_id;});
        assert(algorithm != algorithms.end());

        return algorithm;
    }


    algorithm_list::const_iterator
    find_by_name(const std::string& a_name)
    {
        algorithm_list::const_iterator algorithm =
            std::find_if(algorithms.begin(),
                         algorithms.end(),
                         [&a_name](const algorithm_list::value_type& x) {return x->name() == a_name;});

        return algorithm;
    }


    ////////////////////////////////////////////////////////////////////////


    static const std::string range_separator("..");      //< layer-range-separator ..
    static const std::string empty_index_symbol("_");    //< layer-range-empty-index-symbol \char95
    static const std::string reverse_keyword("reverse"); //< layer-range-reverse-keyword reverse


    struct Index
    {
        virtual ~Index() {}
        virtual Index* clone() const = 0;
        virtual index_list_t values(int a_maximum_index) const = 0;
        virtual std::string as_string() const = 0;
    };


    struct SingletonIndex : public Index
    {
        virtual int value(int a_forward_index) const = 0;
        virtual void set_first_position(bool is_first) {}
        virtual bool get_first_position() const {assert(false); return true;}

        index_list_t values(int a_maximum_index) const
        {
            return index_list_t(1, value(a_maximum_index));
        }
    };


    class ForwardIndex : public SingletonIndex
    {
    public:
        ForwardIndex() = delete;
        ForwardIndex(int a_forward_index) : index_(a_forward_index) {}
        ForwardIndex* clone() const {return new ForwardIndex(*this);}
        int value(int) const {return index_;}
        std::string as_string() const {return std::to_string(index_);}

    private:
        const int index_;
    };


    class BackwardsIndex : public SingletonIndex
    {
    public:
        BackwardsIndex() = delete;
        BackwardsIndex(int a_backwards_index) : index_(a_backwards_index) {}
        BackwardsIndex* clone() const {return new BackwardsIndex(*this);}
        int value(int a_maximum_index) const {return a_maximum_index + index_ + 1;}
        std::string as_string() const {return std::to_string(index_);}

    private:
        const int index_;
    };


    class EmptyIndex : public SingletonIndex
    {
    public:
        EmptyIndex* clone() const override {return new EmptyIndex(*this);}
        int value(int a_maximum_index) const override {return get_first_position() ? 1 : a_maximum_index;}
        std::string as_string() const override {return empty_index_symbol;}
        void set_first_position(bool is_first) override {is_first_ = is_first;}
        bool get_first_position() const override {return is_first_;}

    private:
        bool is_first_;
    };


    class RangeOfIndices : public Index
    {
    public:
        RangeOfIndices() = delete;
        RangeOfIndices(const SingletonIndex* a_begin, const SingletonIndex* an_end, int a_stride);
        RangeOfIndices(const RangeOfIndices& another_range_of_indices);
        RangeOfIndices& operator=(const RangeOfIndices&) = delete;
        ~RangeOfIndices();
        RangeOfIndices* clone() const {return new RangeOfIndices(*this);}
        std::string as_string() const;
        index_list_t values(int a_maximum_index) const;

    private:
        SingletonIndex* begin_;
        SingletonIndex* end_;
        const int stride_;
    };


    RangeOfIndices::RangeOfIndices(const SingletonIndex* a_begin, const SingletonIndex* an_end, int a_stride) :
        begin_(static_cast<SingletonIndex*>(a_begin->clone())),
        end_(static_cast<SingletonIndex*>(an_end->clone())),
        stride_(a_stride)
    {
        assert(a_stride != 0);
        begin_->set_first_position(true);
        end_->set_first_position(false);
    }


    RangeOfIndices::RangeOfIndices(const RangeOfIndices& another_range_of_indices) :
        begin_(static_cast<SingletonIndex*>(another_range_of_indices.begin_->clone())),
        end_(static_cast<SingletonIndex*>(another_range_of_indices.end_->clone())),
        stride_(another_range_of_indices.stride_)
    {}


    RangeOfIndices::~RangeOfIndices()
    {
        delete begin_;
        delete end_;
    }


    std::string
    RangeOfIndices::as_string() const
    {
        const std::string direction(stride_ >= 1 ? "" : reverse_keyword + " ");

        return direction + begin_->as_string() + range_separator + end_->as_string();
    }


    index_list_t
    RangeOfIndices::values(int a_maximum_index) const
    {
        const int begin = begin_->value(a_maximum_index);
        const int end = end_->value(a_maximum_index);
        const int actual_begin = stride_ >= 1 ? begin : end;
        const int actual_end = stride_ >= 1 ? end : begin;
#ifdef DEBUG_FILESPEC
        std::cout <<
            "+ RangeOfIndices::values(" << a_maximum_index << "): actual_begin = " << actual_begin <<
            ", actual_end = " << actual_end << ", stride = " << stride_ << "\n";
#endif
        index_list_t result;

        int i = actual_begin;
        while (stride_ >= 1 ? i <= actual_end : i >= actual_end)
        {
            result.push_back(i);
            i += stride_;
        }

        return result;
    }


    // http://www.boost.org/doc/libs/1_58_0/libs/tokenizer/
    typedef boost::char_separator<char> separator_t;
    typedef boost::tokenizer<separator_t> tokenizer_t;


    inline static
    bool
    is_whitespace_or_empty(const std::string& a_string)
    {
        for (auto c : a_string)
        {
            if (!(std::isspace(c) || c == empty_index_symbol[0]))
            {
                return false;
            }
        }
        return true;
    }


    static SingletonIndex*
    parse_singleton(const std::string& a_token)
    {
        if (is_whitespace_or_empty(a_token))
        {
            return new EmptyIndex();
        }

        try
        {
            std::unique_ptr<std::size_t> tail_position(new std::size_t);
            const int n = std::stoi(a_token, tail_position.get());

            // The condition below reads: ``If the tail-position,
            // i.e. the first character we did not parse, equals the
            // size of the token, we have really parsed the token in
            // its entirety.''  Note: This is safer than to code
            //     a_token[*tail_position] == 0
            // which relies on a terminal character promotable to 0.
            if (*tail_position == a_token.size())
            {
                if (n >= 1)
                {
                    return new ForwardIndex(n);
                }
                else if (n <= -1)
                {
                    return new BackwardsIndex(n);
                }
                else
                {
                    std::cerr << command << ": zero (\"" << a_token << "\") is illegal as a layer index\n";
                }
            }
            else
            {
                std::cerr <<
                    command << ": trailing garbage \"" << a_token.substr(*tail_position) <<
                    "\" after integer " << n << "\n";
            }
        }
        catch (std::invalid_argument)
        {
            std::cerr << command << ": invalid integer \"" << a_token << "\"\n";
        }
        catch (std::out_of_range)
        {
            std::cerr << command << ": number " << a_token << " is outside of the range of integers\n";
        }

        return nullptr;
    }


    LayerSpecification::LayerSpecification(const std::string& a_layer_specification)
    {
        static const separator_t separator(NUMERIC_OPTION_DELIMITERS, "");
        tokenizer_t tokenizer(a_layer_specification, separator);
        int token_number = 1;

        for (tokenizer_t::iterator token = tokenizer.begin(); token != tokenizer.end(); ++token)
        {
            const std::string::size_type range_separator_position = token->find(range_separator);
            std::pair<std::string::const_iterator, std::string::const_iterator> prefix_count =
                std::mismatch(token->begin(), token->end(), reverse_keyword.begin(),
                              [](std::string::value_type v1, std::string::value_type v2)
                              {return tolower(v1) == tolower(v2);});

            if (range_separator_position == std::string::npos)
            {
                if (prefix_count.first > token->begin())
                {
                    std::cerr <<
                        command << ": layer index, #" << token_number <<
                        " does not take \"reverse\" keyword in layer specification \"" <<
                        a_layer_specification << "\"\n";
                    exit(1);
                }

                SingletonIndex* possible_index(parse_singleton(*token));
                if (possible_index == nullptr)
                {
                    std::cerr <<
                        command << ": info: cannot parse layer index, #" << token_number <<
                        " in layer specification \"" << a_layer_specification << "\"\n";
                    exit(1);
                }
                else
                {
                    indices_.push_back(possible_index);
                }
            }
            else
            {
                const size_t prefix_length = prefix_count.first - token->begin();
                const bool forward_range = prefix_count.first == token->begin();
#ifdef DEBUG_FILESPEC
                std::cout <<
                    "+ LayerSpecification::LayerSpecification(const std::string&): " <<
                    " token = <" << *token << ">, reverse_keyword = <" << reverse_keyword << ">\n" <<
                    "+ LayerSpecification::LayerSpecification(const std::string&): " <<
                    "prefix_length = " << prefix_length <<
                    ", 1st = <" << token->substr(prefix_length, range_separator_position - prefix_length) <<
                    ">, 2nd = <" << token->substr(range_separator_position + range_separator.length()) << ">\n";
#endif

                std::unique_ptr<SingletonIndex>
                    first_index(parse_singleton(token->substr(prefix_length,
                                                              range_separator_position - prefix_length)));
                std::unique_ptr<SingletonIndex>
                    second_index(parse_singleton(token->substr(range_separator_position +
                                                               range_separator.length())));

                if (first_index.get() == nullptr)
                {
                    std::cerr <<
                        command << ": info: cannot parse first bound of layer range, #" << token_number <<
                        " in layer specification \"" << a_layer_specification << "\"\n";
                    exit(1);
                }
                else if (second_index.get() == nullptr)
                {
                    std::cerr <<
                        command << ": info: cannot parse second bound of layer range, #" << token_number <<
                        " in layer specification \"" << a_layer_specification << "\"\n";
                    exit(1);
                }
                else
                {
                    if (forward_range)
                    {
                        indices_.push_back(new RangeOfIndices(first_index.get(), second_index.get(), 1));
                    }
                    else
                    {
                        indices_.push_back(new RangeOfIndices(first_index.get(), second_index.get(), -1));
                    }
                }
            }

            ++token_number;
        }
    }


    LayerSpecification::LayerSpecification(const LayerSpecification& another_layer_specification) :
        indices_(abstract_index_list_t(another_layer_specification.indices_.size(), nullptr))
    {
        std::transform(another_layer_specification.indices_.begin(),
                       another_layer_specification.indices_.end(),
                       indices_.begin(),
                       [](Index* x) {return x->clone();});
    }


    LayerSpecification::~LayerSpecification()
    {
        std::for_each(indices_.begin(), indices_.end(), [](Index* x) {delete x;});
    }


    index_list_t
    LayerSpecification::values(int a_maximum_index) const
    {
        index_list_t result;

        for (auto i : indices_)
        {
            const index_list_t is(i->values(a_maximum_index));
            result.insert(result.end(), is.begin(), is.end());
        }

        return result;
    }


    std::string
    LayerSpecification::as_string() const
    {
        return
            indices_.empty() ?
            std::string("<empty>") :
            mapconcat(indices_.begin(), indices_.end(), [](const Index* x) {return x->as_string();}, ":");
    }


    ////////////////////////////////////////////////////////////////////////


    IndexedLayer::IndexedLayer(const std::string& a_layer_specification) :
        layer_spec_(LayerSpecification(a_layer_specification))
    {}


    std::string
    IndexedLayer::name() const
    {
        return std::string("indexed-layer ") + layer_spec_.as_string();
    }


    std::string
    IndexedLayer::description() const
    {
        return "select layers from image by a given tuple of indexes;";
    }


    LayerSpecification
    IndexedLayer::layer_spec() const
    {
        return layer_spec_;
    }


    bool
    IndexedLayer::select(const ImageListInformation* an_image_info,
                         const std::string& a_filename, unsigned a_layer_index)
    {
        const ImageInfo* image_info = an_image_info->image_info_on(a_filename);
        const unsigned number_of_layers = image_info->number_of_layers();
        const index_list_t selection = layer_spec_.values(number_of_layers);

        if (selection.empty())
        {
            std::cerr <<
                command << ": info: indexed layer selection " << layer_spec_.as_string() <<
                " contains no index\n";
        }

        if (a_layer_index > number_of_layers)
        {
            std::cerr <<
                command << ": warning: " <<
                "layer index " << a_layer_index << " is larger than number of layers " <<
                number_of_layers << "\n";
        }

        return std::find(selection.begin(), selection.end(), a_layer_index) != selection.end();
    }


    layer_ordered_list_t
    IndexedLayer::viable_layers(const ImageListInformation* an_image_info,
                                const std::string& a_filename)
    {
        const ImageInfo* image_info = an_image_info->image_info_on(a_filename);
        const unsigned number_of_layers = image_info->number_of_layers();
        const index_list_t selection = layer_spec_.values(number_of_layers);

        layer_ordered_list_t result;

        std::copy(selection.begin(), selection.end(), std::back_inserter(result));

        return result;
    }
} // end namespace selector
