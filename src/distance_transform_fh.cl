// Copyright (C) 2013-2017 Christoph L. Spiel
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


void
manhattan_1d(global const float *restrict f, const int n, global float *restrict d)
{
    d[0] = f[0];

#pragma unroll 8
    for (int q = 1; q < n; ++q)
    {
        d[q] = min(f[q], d[q - 1] + 1.0f);
    }

#pragma unroll 8
    for (int q = n - 2; q >= 0; --q)
    {
        d[q] = min(d[q], d[q + 1] + 1.0f);
    }
}


kernel void
manhattan_2d_columns(global float *restrict output,
                     const int width, const int height,
                     global float *restrict f_base, global float *restrict d_base)
{
    const int x = get_global_id(0);

    if (x >= width)
    {
        return;
    }

    const int offset = x * max(width, height);
    global float *f = f_base + offset;
    global float *d = d_base + offset;

    for (int y = 0; y < height; y++)
    {
        f[y] = output[x + y * width];
    }

    manhattan_1d(f, height, d);

    for (int y = 0; y < height; y++)
    {
        output[x + y * width] = d[y];
    }
}


kernel void
manhattan_2d_rows(global float *restrict output,
                  const int width, const int height,
                  global float *restrict f_base, global float *restrict d_base)
{
    const int y = get_global_id(0);

    if (y >= height)
    {
        return;
    }

    const int offset = y * max(width, height);
    global float *f = f_base + offset;
    global float *d = d_base + offset;

    for (int x = 0; x < width; x++)
    {
        f[x] = output[x + y * width];
    }

    manhattan_1d(f, width, d);

    for (int x = 0; x < width; x++)
    {
        output[x + y * width] = d[x];
    }
}


////////////////////////////////////////////////////////////////////////////////


float
square(const int n)
{
    const float x = (float) n;

    return x * x;
}


void
euclidean_1d(global const float *restrict f, const int n, global float *restrict d,
             global int *restrict v, global float *restrict z)
{
    v[0] = 0;
    z[0] = -FLT_MAX;
    z[1] = +FLT_MAX;

    int k = 0;
#pragma unroll 8
    for (int q = 1; q < n; q++)
    {
        const float a = f[q] + square(q);
        float s = 0.5f * (a - (f[v[k]] + square(v[k]))) / (float) (q - v[k]);
        while (s <= z[k])
        {
            k--;
            s = 0.5f * (a - (f[v[k]] + square(v[k]))) / (float) (q - v[k]);
        }
        k++;

        v[k] = q;
        z[k] = s;
        z[k + 1] = +FLT_MAX;
    }

    k = 0;
#pragma unroll 8
    for (int q = 0; q < n; q++)
    {
        while (z[k + 1] < (float) q)
        {
            k++;
        }
        d[q] = square(q - v[k]) + f[v[k]];
    }
}


kernel void
euclidean_2d_columns(global float *restrict output,
                     const int width, const int height,
                     global float *restrict f_base, global float *restrict d_base,
                     global int *restrict v_base, global float *restrict z_base)
{
    const int x = get_global_id(0);

    if (x >= width)
    {
        return;
    }

    const int offset = x * max(width, height);
    global float *f = f_base + offset;
    global float *d = d_base + offset;
    global int *v = v_base + offset;
    global float *z = z_base + offset;

    for (int y = 0; y < height; y++)
    {
        f[y] = output[x + y * width];
    }

    euclidean_1d(f, height, d, v, z);

    for (int y = 0; y < height; y++)
    {
        output[x + y * width] = d[y];
    }
}


kernel void
euclidean_2d_rows(global float *restrict output,
                  const int width, const int height,
                  global float *restrict f_base, global float *restrict d_base,
                  global int *restrict v_base, global float *restrict z_base)
{
    const int y = get_global_id(0);

    if (y >= height)
    {
        return;
    }

    const int offset = y * max(width, height);
    global float *f = f_base + offset;
    global float *d = d_base + offset;
    global int *v = v_base + offset;
    global float *z = z_base + offset;

    for (int x = 0; x < width; x++)
    {
        f[x] = output[x + y * width];
    }

    euclidean_1d(f, width, d, v, z);

    for (int x = 0; x < width; x++)
    {
        output[x + y * width] = native_sqrt(d[x]);
    }
}
