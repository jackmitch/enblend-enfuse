// Copyright (C) 2015, 2016 Christoph L. Spiel
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
new_state_probabilities(const int k, local float *sp, constant float *e, local float *pi,
                        const int total_scratch_size, local float *scratch)
{
    const int gid = get_global_id(0);

    // We use `scratch' for two arrays: `reciprocal_exp_delta_e' and
    // `reduce_scratch'.  Each gets half of the available space.
    const int scratch_size = total_scratch_size / 2;
    local float *reciprocal_exp_delta_e = scratch;
    local float *reduce_scratch = scratch + scratch_size;

    for (int j = 0; j < k; j++)
    {
        const float x = native_divide(1.0f, 1.0f + native_exp(e[j] - e[gid]));
#ifdef __FAST_RELAXED_MATH__
        reciprocal_exp_delta_e[gid] = x;
#else
        reciprocal_exp_delta_e[gid] = isnan(x) ? (e[j] > e[gid] ? 0.0f : 1.0f) : x;
#endif
        barrier(CLK_LOCAL_MEM_FENCE);

        const float sp_j = sp[j] + sp[gid];
        reduce_scratch[gid] = (gid <= j || gid >= k) ? 0.0f : sp_j * reciprocal_exp_delta_e[gid];
        pi[gid] += (gid <= j || gid >= k) ? 0.0f : sp_j - reduce_scratch[gid];
        barrier(CLK_LOCAL_MEM_FENCE);

        for (int s = scratch_size / 2; s > 0; s >>= 1)
        {
            if (gid < s)
            {
                reduce_scratch[gid] += reduce_scratch[gid + s];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        if (gid == 0)
        {
            const float pi_j = sp[j] + pi[j] + reduce_scratch[0];

            pi[j] = pi_j;
            sp[j] = pi_j / (float) k;
        }
    }
}


#ifdef HAVE_EXTENSION_CL_KHR_FP64

#pragma OPENCL EXTENSION cl_khr_fp64: enable

kernel void
calculate_state_probabilities(const int k,
                              global double *restrict global_sp, local float *sp,
                              constant float *e,
                              global float *restrict global_pi, local float *pi,
                              const int total_scratch_size, local float *scratch)
{
    const int gid = get_global_id(0);

    if (gid < k)
    {
        sp[gid] = (float) global_sp[gid];
        pi[gid] = global_pi[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    new_state_probabilities(k, sp, e, pi, total_scratch_size, scratch);

    barrier(CLK_LOCAL_MEM_FENCE);
    if (gid < k)
    {
        global_sp[gid] = (double) sp[gid];
        global_pi[gid] = pi[gid];
    }
}

#else

kernel void
calculate_state_probabilities(const int k,
                              global float *restrict global_sp, local float *sp,
                              constant float *e,
                              global float *restrict global_pi, local float *pi,
                              const int total_scratch_size, local float *scratch)
{
    const int gid = get_global_id(0);

    if (gid < k)
    {
        sp[gid] = global_sp[gid];
        pi[gid] = global_pi[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    new_state_probabilities(k, sp, e, pi, total_scratch_size, scratch);

    barrier(CLK_LOCAL_MEM_FENCE);
    if (gid < k)
    {
        global_sp[gid] = sp[gid];
        global_pi[gid] = pi[gid];
    }
}

#endif // HAVE_EXTENSION_CL_KHR_FP64
