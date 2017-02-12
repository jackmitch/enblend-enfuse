/*
 * Copyright (C) 2015-2017 Christoph L. Spiel
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
#ifndef INTROSPECTION_H_INCLUDED_
#define INTROSPECTION_H_INCLUDED_


namespace introspection
{
    void printVersion(int argc, char** argv);
    void printImageFormats();
    void printSignature();
    void printGlobbingAlgos();
    void printSoftwareComponents();
} // end namespace introspection


#endif // INTROSPECTION_H_INCLUDED_


// Local Variables:
// mode: c++
// End:
