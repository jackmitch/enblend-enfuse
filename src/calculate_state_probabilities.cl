// Copyright (C) 2015 Christoph L. Spiel
//
// This file is part of Enblend.
//
// Enblend is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Enblend is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Enblend; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


inline static
void
new_state_probabilities(const int k, local float *sp, constant float *e, local float *pi)
{
    const int gid = get_global_id(0);

    local float recip_exp_delta_e[32];

    for (int j = 0; j < k; j++)
    {
        const float pi_t_j = sp[j];
        float pi_j = pi_t_j + pi[j];

        recip_exp_delta_e[gid] = native_divide(1.0f, 1.0f + native_exp(e[j] - e[gid]));
        barrier(CLK_LOCAL_MEM_FENCE);

        if (gid == 0)
        {
            for (int i = j + 1; i < k; i++)
            {
                const float pi_t = sp[i] + pi_t_j;
                float pi_t_an = pi_t * recip_exp_delta_e[i];
#ifndef __FAST_RELAXED_MATH__
                pi_t_an = isnan(pi_t_an) ? (e[j] > e[i] ? 0.0f : pi_t) : pi_t_an;
#endif
                pi_j += pi_t_an;
                pi[i] += pi_t - pi_t_an;
            }

            pi[j] = pi_j;
            sp[j] = pi_j / (float) k;
        }
    }
}


#pragma OPENCL EXTENSION cl_khr_fp64: enable


kernel void
calculate_state_probabilities(const int k,
                              global double *restrict global_sp, local float *sp,
                              constant float *e,
                              global float *restrict global_pi, local float *pi)
{
    const int gid = get_global_id(0);

    if (gid >= k)
    {
        return;
    }

    sp[gid] = (float) global_sp[gid];
    pi[gid] = global_pi[gid];

    barrier(CLK_LOCAL_MEM_FENCE);

    new_state_probabilities(k, sp, e, pi);

    barrier(CLK_LOCAL_MEM_FENCE);

    global_sp[gid] = (double) sp[gid];
    global_pi[gid] = pi[gid];
}
