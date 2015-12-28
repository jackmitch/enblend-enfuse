/*
 * Copyright (C) 2010-2015 Dr. Christoph L. Spiel
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
#ifndef SELECTOR_H_
#define SELECTOR_H_


#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "info.h"
#include "layer_selection.h"


namespace selector
{
    enum class id_t
    {
        AllLayersId,
        FirstLayerId,
        LastLayerId,
        LargestLayerId,
        NoLayerId,
        IndexedLayerId
    };


    typedef std::vector<unsigned> layer_ordered_list_t;


    struct Abstract
    {
        Abstract() {}
        virtual ~Abstract() {}
        virtual id_t id() const = 0;

        virtual std::string name() const = 0;
        virtual std::string description() const = 0;

        virtual bool select(const ImageListInformation* an_image_info,
                            const std::string& a_filename,
                            unsigned a_layer_index) = 0;
        virtual layer_ordered_list_t viable_layers(const ImageListInformation* an_image_info,
                                                   const std::string& a_filename) = 0;
    };


    class AllLayers : public Abstract
    {
    public:
        id_t id() const {return id_t::AllLayersId;}
        std::string name() const;
        std::string description() const;

        bool select(const ImageListInformation*, const std::string&, unsigned) {return true;}
        layer_ordered_list_t viable_layers(const ImageListInformation* an_image_info,
                                           const std::string& a_filename);
    };


    class FirstLayer : public Abstract
    {
    public:
        id_t id() const {return id_t::FirstLayerId;}
        std::string name() const;
        std::string description() const;

        bool select(const ImageListInformation*, const std::string&, unsigned a_layer_index)
        {
            return a_layer_index == 1;
        }

        layer_ordered_list_t viable_layers(const ImageListInformation*, const std::string&)
        {
            return layer_ordered_list_t(1, 1);
        }
    };


    class LastLayer : public Abstract
    {
    public:
        id_t id() const {return id_t::LastLayerId;}
        std::string name() const;
        std::string description() const;

        bool select(const ImageListInformation* an_image_info,
                    const std::string& a_filename, unsigned a_layer_index);

        layer_ordered_list_t viable_layers(const ImageListInformation* an_image_info,
                                           const std::string& a_filename);
    };


    class LargestLayer : public Abstract
    {
        typedef std::map<std::string, unsigned> cache_t;

    public:
        id_t id() const {return id_t::LargestLayerId;}
        std::string name() const;
        std::string description() const;

        bool select(const ImageListInformation* an_image_info,
                    const std::string& a_filename, unsigned a_layer_index);
        layer_ordered_list_t viable_layers(const ImageListInformation* an_image_info,
                                           const std::string& a_filename);

    private:
        cache_t cache_;
    };


    class NoLayer : public Abstract
    {
    public:
        id_t id() const {return id_t::NoLayerId;}
        std::string name() const;
        std::string description() const;

        bool select(const ImageListInformation*, const std::string&, unsigned) {return false;}

        layer_ordered_list_t viable_layers(const ImageListInformation*, const std::string&)
        {
            return layer_ordered_list_t();
        }
    };


    typedef std::list<std::unique_ptr<Abstract> > algorithm_list;

    extern algorithm_list algorithms;

    algorithm_list::const_iterator find_by_id(id_t an_id);
    algorithm_list::const_iterator find_by_name(const std::string& a_name);


    ////////////////////////////////////////////////////////////////////////


    typedef std::vector<int> index_list_t;

    struct Index;               // forward declaration

    class LayerSpecification
    {
    public:
        LayerSpecification() = delete;
        LayerSpecification(const std::string& a_layer_specification);
        LayerSpecification(const LayerSpecification&);
        LayerSpecification& operator=(const LayerSpecification&) = delete;
        ~LayerSpecification();

        index_list_t values(int a_maximum_index) const;
        std::string as_string() const;

    private:
        typedef std::vector<Index*> abstract_index_list_t;
        abstract_index_list_t indices_;
    };


    class IndexedLayer : public Abstract
    {
    public:
        IndexedLayer() = delete;
        explicit IndexedLayer(const std::string& a_layer_specification);

        id_t id() const {return id_t::IndexedLayerId;}
        std::string name() const;
        std::string description() const;
        LayerSpecification layer_spec() const;

        bool select(const ImageListInformation* an_image_info,
                    const std::string& a_filename, unsigned a_layer_index);
        layer_ordered_list_t viable_layers(const ImageListInformation* an_image_info,
                                           const std::string& a_filename);

    private:
        const LayerSpecification layer_spec_;
    };
} // end namespace selector


#endif // SELECTOR_H_


// Local Variables:
// mode: c++
// End:
