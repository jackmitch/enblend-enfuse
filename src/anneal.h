/*
 * Copyright (C) 2004 Andrew Mihal
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
#ifndef __ANNEAL_H__
#define __ANNEAL_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <vector>

#include "vigra/diff2d.hxx"

using std::pair;
using std::vector;

using vigra::Point2D;

namespace enblend {

template <typename DistanceCostImage, typename StitchCostImage>
class AnnealConfiguration {
public:
    typedef DistanceCostImage DistanceCostImageType;
    typedef StitchCostImage StitchCostImageType;

    AnnealConfiguration()
    : snake(), dci(NULL), sci(NULL) { }

    AnnealConfiguration(const DistanceCostImage *d,
            const StitchCostImage *s,
            const vector<pair<bool, Point2D> > &v)
    : snake(v), dci(d), sci(s) { }

    AnnealConfiguration(const AnnealConfiguration &c)
    : snake(c.snake), dci(c.dci), sci(c.sci) { }

    ~AnnealConfiguration() { }

    AnnealConfiguration& operator=(const AnnealConfiguration& c) {
        snake.clear();
        snake.insert(snake.end(), c.snake.begin(), c.snake.end());
        dci = c.dci;
        sci = c.sci;
        return *this;
    }

    double eFunc() const;
    void step(const gsl_rng *r, double step_size);
    double operator-(const AnnealConfiguration& c) const;
    void print() const;

    const vector<pair<bool, Point2D> >& getSnake() const {
        return snake;
    }

protected:
    vector<pair<bool, Point2D> > snake;
    const DistanceCostImage *dci;
    const StitchCostImage *sci;
};

template <typename DistanceCostImage, typename StitchCostImage>
double AnnealConfiguration<DistanceCostImage, StitchCostImage>::eFunc() const {
    return 0.0;
};

template <typename DistanceCostImage, typename StitchCostImage>
void AnnealConfiguration<DistanceCostImage, StitchCostImage>::step(const gsl_rng *r, double step_size) {
    for (unsigned int i = 0; i < snake.size(); i++) {
        pair<bool, Point2D> p = snake[i];
        if (!p.first) {
            // Point is not frozen.
            // FIXME do not allow point to leave dci bounds.
            int deltaX = (int)floor((gsl_rng_uniform(r) * step_size) - (step_size / 2.0));
            int deltaY = (int)floor((gsl_rng_uniform(r) * step_size) - (step_size / 2.0));
            snake[i] = pair<bool, Point2D>(false, p.second + Diff2D(deltaX, deltaY));
        }
    }
    //cout << "generated step:";
    //print();
    //cout << endl;
    return;
};

template <typename DistanceCostImage, typename StitchCostImage>
void AnnealConfiguration<DistanceCostImage, StitchCostImage>::print() const {
    for (unsigned int i = 0; i < snake.size(); i++) {
        pair<bool, Point2D> p = snake[i];
        if (!p.first) {
            cout << " (" << p.second.x << ", " << p.second.y << ")";
        }
    }
    return;
};

template <typename DistanceCostImage, typename StitchCostImage>
double AnnealConfiguration<DistanceCostImage, StitchCostImage>::operator-(const AnnealConfiguration &c) const {
    cout << "operator-" << endl;
    double totalDistance = 0.0;
    for (unsigned int i = 0; i < snake.size(); i++) {
        pair<bool, Point2D> p1 = snake[i];
        pair<bool, Point2D> p2 = c.snake[i];
        Diff2D diff = p1.second - p2.second;
        totalDistance += diff.magnitude();
    }
    return totalDistance;
};
    
template <typename AnnealConfigurationType>
double snakeEFunc(void *xp) {
    AnnealConfigurationType *ac = (AnnealConfigurationType*)xp;
    return ac->eFunc();
};

template <typename AnnealConfigurationType>
void snakeStep(const gsl_rng *r, void *xp, double step_size) {
    AnnealConfigurationType *ac = (AnnealConfigurationType*)xp;
    ac->step(r, step_size);
};

template <typename AnnealConfigurationType>
double snakeMetric(void *xp, void *yp) {
    AnnealConfigurationType *acx = (AnnealConfigurationType*)xp;
    AnnealConfigurationType *acy = (AnnealConfigurationType*)yp;
    return (*acx) - (*acy);
};

template <typename AnnealConfigurationType>
void snakePrint(void *xp) {
    AnnealConfigurationType *ac = (AnnealConfigurationType*)xp;
    ac->print();
};

template <typename AnnealConfigurationType>
void snakeCopy(void *source, void *dest) {
    AnnealConfigurationType *acx = (AnnealConfigurationType*)source;
    AnnealConfigurationType *acy = (AnnealConfigurationType*)dest;
    (*acy) = (*acx);
};

template <typename AnnealConfigurationType>
void *snakeCopyConstruct(void *xp) {
    AnnealConfigurationType *ac = (AnnealConfigurationType*)xp;
    return new AnnealConfigurationType(*ac);
};

template <typename AnnealConfigurationType>
void snakeDestroy(void *xp) {
    AnnealConfigurationType *ac = (AnnealConfigurationType*)xp;
    delete ac;
};

template <typename DistanceCostImage, typename StitchCostImage>
void annealSnake(const DistanceCostImage* const dci,
        const StitchCostImage* const sci,
        vector<pair<bool, Point2D> > *snake) {

    typedef AnnealConfiguration<DistanceCostImage, StitchCostImage> AnnealConfigurationType;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_default_seed++;

    AnnealConfigurationType initConfig(dci, sci, *snake);

    int n_tries = 200;
    int iters_fixed_T = 10;
    double step_size  = 10.0;
    double k = 1.0;
    double t_initial = 10.0; //initConfig.eFunc() / 2.0;
    double t_final = 0.1; //initConfig.eFunc() / 100.0;
    double mu_t = 1.05;
    
    gsl_siman_params_t params = {n_tries, iters_fixed_T, step_size, k, t_initial, mu_t, t_final};

    gsl_siman_solve(rng,
            &initConfig,
            snakeEFunc<AnnealConfigurationType>,
            snakeStep<AnnealConfigurationType>,
            snakeMetric<AnnealConfigurationType>,
            snakePrint<AnnealConfigurationType>,
            snakeCopy<AnnealConfigurationType>,
            snakeCopyConstruct<AnnealConfigurationType>,
            snakeDestroy<AnnealConfigurationType>,
            0, params);

    // Copy result to input snake.
    snake->clear();
    snake->insert(snake->end(), initConfig.getSnake().begin(), initConfig.getSnake().end());

    gsl_rng_free(rng);
};

} // namespace enblend

#endif /* __ANNEAL_H__ */
