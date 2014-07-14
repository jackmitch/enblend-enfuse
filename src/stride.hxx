#ifndef STRIDE_HXX_INCLUDED
#define STRIDE_HXX_INCLUDED


#include <vigra/imageiterator.hxx>


namespace vigra_ext
{
    namespace detail
    {
        inline static vigra::Diff2D
        strided_size(int an_xstride, int a_ystride, vigra::Diff2D a_size)
        {
            if (a_size.x % an_xstride != 0)
            {
                a_size.x += an_xstride - a_size.x % an_xstride;
            }
            if (a_size.y % a_ystride != 0)
            {
                a_size.y += a_ystride - a_size.y % a_ystride;
            }

            a_size.x /= an_xstride;
            a_size.y /= a_ystride;

            return a_size;
        }


        template <class iterator>
        inline static int
        iterator_width(iterator an_image_iterator)
        {
            iterator next_line(an_image_iterator);

            next_line.y += 1;

            return next_line[0] - an_image_iterator[0];
        }
    } // namespace detail


    template <typename pixel_type, typename accessor, typename iterator>
    inline static vigra::triple<vigra::StridedImageIterator<pixel_type>,
                                vigra::StridedImageIterator<pixel_type>,
                                accessor>
    stride(int an_xstride, int a_ystride,
           const vigra::triple<vigra::BasicImageIterator<pixel_type, iterator>,
                               vigra::BasicImageIterator<pixel_type, iterator>,
                               accessor>& an_image)
    {
        typedef vigra::StridedImageIterator<pixel_type> strided_iterator;

        const vigra::Diff2D size(detail::strided_size(an_xstride, a_ystride, an_image.second - an_image.first));
        const strided_iterator base(strided_iterator(an_image.first[0],
                                                     detail::iterator_width(an_image.first),
                                                     an_xstride, a_ystride));

        return vigra::make_triple(base, base + size, an_image.third);
    }


    template <typename pixel_type, typename accessor>
    inline static vigra::triple<vigra::ConstStridedImageIterator<pixel_type>,
                                vigra::ConstStridedImageIterator<pixel_type>,
                                accessor>
    stride(int an_xstride, int a_ystride,
           const vigra::triple<vigra::ConstImageIterator<pixel_type>,
                               vigra::ConstImageIterator<pixel_type>,
                               accessor>& an_image)
    {
        typedef vigra::ConstStridedImageIterator<pixel_type> constant_strided_iterator;

        const vigra::Diff2D size(detail::strided_size(an_xstride, a_ystride, an_image.second - an_image.first));
        const constant_strided_iterator base(an_image.first[0], detail::iterator_width(an_image.first),
                                             an_xstride, a_ystride);

        return vigra::make_triple(base, base + size, an_image.third);
    }


    template <typename pixel_type, typename accessor, typename iterator>
    inline static std::pair<vigra::StridedImageIterator<pixel_type>, accessor>
    stride(int an_xstride, int a_ystride,
           const std::pair<vigra::BasicImageIterator<pixel_type, iterator>, accessor>& an_image)
    {
        typedef vigra::StridedImageIterator<pixel_type> strided_iterator;

        return std::make_pair(strided_iterator(an_image.first[0],
                                               detail::iterator_width(an_image.first),
                                               an_xstride, a_ystride),
                              an_image.second);
    }


    template <typename pixel_type, typename accessor>
    inline static std::pair<vigra::StridedImageIterator<pixel_type>, accessor>
    stride(int an_xstride, int a_ystride,
           const std::pair<vigra::ImageIterator<pixel_type>, accessor>& an_image)
    {
        typedef vigra::StridedImageIterator<pixel_type> strided_iterator;

        return std::make_pair(strided_iterator(an_image.first[0], detail::iterator_width(an_image.first),
                                               an_xstride, a_ystride),
                              an_image.second);
    }


    template <typename pixel_type, typename accessor, typename iterator>
    inline static vigra::triple<vigra::ConstStridedImageIterator<pixel_type>,
                                vigra::ConstStridedImageIterator<pixel_type>,
                                accessor>
    stride(int an_xstride, int a_ystride,
           const vigra::triple<vigra::ConstBasicImageIterator<pixel_type, iterator>,
                               vigra::ConstBasicImageIterator<pixel_type, iterator>,
                               accessor>& an_image)
    {
        typedef vigra::ConstStridedImageIterator<pixel_type> constant_strided_iterator;

        const vigra::Diff2D size(detail::strided_size(an_xstride, a_ystride, an_image.second - an_image.first));
        const constant_strided_iterator base(constant_strided_iterator(an_image.first[0],
                                                                       detail::iterator_width(an_image.first),
                                                                       an_xstride, a_ystride));

        return vigra::make_triple(base, base + size, an_image.third);
    }


    template <typename pixel_type, typename accessor, typename iterator>
    inline static std::pair<vigra::ConstStridedImageIterator<pixel_type>, accessor>
    stride(int an_xstride, int a_ystride,
           const std::pair<vigra::ConstBasicImageIterator<pixel_type, iterator>, accessor>& an_image)
    {
        typedef vigra::ConstStridedImageIterator<pixel_type> constant_strided_iterator;

        return std::make_pair(constant_strided_iterator(an_image.first[0],
                                                        detail::iterator_width(an_image.first),
                                                        an_xstride, a_ystride),
                              an_image.second);
    }


    template <typename pixel_type, typename accessor>
    inline static std::pair<vigra::ConstStridedImageIterator<pixel_type>, accessor>
    stride(int an_xstride, int a_ystride,
           const std::pair<vigra::ConstImageIterator<pixel_type>, accessor>& an_image)
    {
        typedef vigra::ConstStridedImageIterator<pixel_type> constant_strided_iterator;

        return std::make_pair(constant_strided_iterator(an_image.first[0],
                                                        detail::iterator_width(an_image.first),
                                                        an_xstride, a_ystride),
                              an_image.second);
    }
} // namespace vigra_ext


#endif // STRIDE_HXX_INCLUDED
