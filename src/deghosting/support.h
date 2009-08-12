
#ifndef SUPPORT_H_
#define SUPPORT_H_

#include "deghosting.h"

#include <vigra/functorexpression.hxx>
// used in hugin_hdrmerge
// FIXME: move it to the hugin_hdrmerge
#include <vigra/combineimages.hxx>

namespace deghosting {

using namespace vigra;
using namespace vigra::functor;

/** Logarithm functor
 */
template <class PixelType>
class LogarithmFunctor {
    public:
        LogarithmFunctor(PixelType off=0) : offset(off)  {}
        
        PixelType operator()(PixelType const& v) const {
            return std::log(v + offset);
        }
    protected:
        PixelType offset;
};

/** Logarithm functor
 * specialization for RGBValue
 */
template <class ComponentType>
class LogarithmFunctor<RGBValue<ComponentType> > {
    public:
        LogarithmFunctor(ComponentType off=0) : offset(off) {}
        
        RGBValue<ComponentType> operator()(RGBValue<ComponentType> const& v) const {
            RGBValue<ComponentType> retVal;
            retVal[0] = log(v[0] + offset);
            retVal[1] = log(v[1] + offset);
            retVal[2] = log(v[2] + offset);
            //cout << retVal[0] << "," << retVal[1] << "," << retVal[2] << endl;
            return retVal;
        }
    protected:
        ComponentType offset;
};

/** Fuctor to normalize values
 */
template <class PixelType>
class NormalizeFunctor {
    public:
        NormalizeFunctor(PixelType f) : factor(f) {}
        NormalizeFunctor(PixelType oldMaxValue, PixelType newMaxValue) : factor(newMaxValue/oldMaxValue) {}
        
        PixelType operator()(PixelType const &v) const {
            return v*factor;
        }
    protected:
        PixelType factor;
};

/** Fuctor to normalize values
 * specialization for RGBValue
 */
template <class ComponentType>
class NormalizeFunctor<RGBValue<ComponentType> > {
    public:
        NormalizeFunctor(RGBValue<ComponentType> oldMaxValue, RGBValue<ComponentType> newMaxValue) {
            // TODO
        }
        
        RGBValue<ComponentType> operator()(RGBValue<ComponentType> const &v) {
            // TODO
        }
    protected:
        RGBValue<ComponentType> foo;
};

/** Fixed size vector with scalar multiplication and element-wise substraction and addition
 */
template <class T, int SIZE>
class AlgTinyVector
{
    public:
        const T operator[](int i) const {
            return content[i];
        }
        
        T& operator[](int i) {
            return content[i];
        }
        
        const T operator*(const AlgTinyVector<T,SIZE> t) const {
            T retVal = 0;
            for (unsigned int i = 0; i < SIZE; ++i) {
                retVal += t[i] * content[i];
            }
            return retVal;
        }
        
        const AlgTinyVector operator*(const int t) const {
            AlgTinyVector<T,SIZE> retVal;
            for (unsigned int i = 0; i < SIZE; ++i) {
                retVal[i] = t * content[i];
            }
            return retVal;
        }
        
        const AlgTinyVector operator/(const int t) const {
            AlgTinyVector<T,SIZE> retVal;
            for (unsigned int i = 0; i < SIZE; ++i) {
                retVal[i] = content[i] / t;
            }
            return retVal;
        }
        
        const AlgTinyVector operator-(const AlgTinyVector<T,SIZE> t) const {
            AlgTinyVector<T,SIZE> retVal;
            for (unsigned int i = 0; i < SIZE; ++i) {
                retVal[i] = t[i] - content[i];
            }
            return retVal;
        }
        
        const AlgTinyVector operator+(const AlgTinyVector<T,SIZE> t) const {
            AlgTinyVector<T,SIZE> retVal;
            for (unsigned int i = 0; i < SIZE; ++i) {
                retVal[i] = t[i] + content[i];
            }
            return retVal;
        }        
        
        AlgTinyVector & operator=(const AlgTinyVector<T,SIZE> & t) {
            if (this == &t)
                return *this;
            for (unsigned int i = 0; i < SIZE; ++i) {
                content[i] = t[i];
            }
            return *this;
        }
        
        AlgTinyVector & operator=(const TinyVector<T,SIZE> & t) {
            for (unsigned int i = 0; i < SIZE; ++i) {
                content[i] = t[i];
            }
            return *this;
        }
        
    private:
        T content[SIZE];
};

}  // namespace deghosting

#endif /* SUPPORT_H_ */
