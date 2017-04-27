// Copyright (C) 2010, Steffen Knollmann
// Released under the terms of the GNU General Public License version 3.

#ifndef UTIL_CONFIG_H
#define UTIL_CONFIG_H


#ifdef ENABLE_DEBUG
#  undef NDEBUG
#  define DEBUG
#  define XMEM_TRACK_MEM
#else
#  define NDEBUG
#  undef DEBUG
#  undef XMEM_TRACK_MEM
#endif

#endif
