/*
 * This file is part of the LAIK library.
 * Copyright (c) 2017, 2018 Josef Weidendorfer <Josef.Weidendorfer@gmx.de>
 *
 * LAIK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, version 3 or later.
 *
 * LAIK is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "laik-internal.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

// this file implements a layout for 1d vectors

typedef struct _Laik_Layout_Vector Laik_Layout_Vector;
struct _Laik_Layout_Vector {
    Laik_Layout h;
    uint64_t localLength;            // Layout data specific to this custom layout << this equals *count* in Laik_Layout h
    uint64_t numberOfExternalValues; // external values stored after localLength elements
    int32_t offset;                  // offset into allocation buffer. Mapping under the hood, but application programmer needs to specify the offset for now
};


//--------------------------------------------------------------
// interface implementation of lexicographical layout
//

// helper for copy_vector/pack_vector: lexicographical traversal
static 
bool next_idx(Laik_Range *range, Laik_Index *idx)
{
    idx->i[0]++;
    if (idx->i[0] < range->to.i[0])
        return true;
    if (range->space->dims == 1)
        return false;
}

// forward decl
static int64_t offset_vector(Laik_Layout* l, int n, Laik_Index* idx);

// return compact vector layout
static
Laik_Layout_Vector* laik_is_layout_vector(Laik_Layout* l)
{
    if (l->offset == offset_vector)
        return (Laik_Layout_Vector*) l;

    return 0; // not a lexicographical layout
}

// return map number whose ranges contains index <idx>
static
int section_vector(Laik_Layout* l, Laik_Index* idx)
{
    Laik_Layout_Vector* lv = laik_is_layout_vector(l);
    assert(lv);

    // TODO: How to test the upper bound
    if(idx->i[0] >= 0)
        return 0;

    return -1; // not found
}

// section is allocation number
static
int mapno_vector(Laik_Layout* l, int n)
{
    assert(n < l->map_count);
    assert(n == 0); // This layout supports only one mapping (allocation buffer)
    return n;
}


// return offset for <idx> in map <n> of this layout
static
int64_t offset_vector(Laik_Layout* l, int n, Laik_Index* idx)
{
    Laik_Layout_Vector *lv = laik_is_layout_vector(l);
    assert(lv);

    int64_t off = idx->i[0];

    // TODO: How to test the upper bound
    assert(off >= 0);

    return off;
}


static
char* describe_vector(Laik_Layout* l)
{
    static char s[200];

    Laik_Layout_Vector *lv = laik_is_layout_vector(l);
    assert(lv);

    int o;
    o = sprintf(s, "compact vector (%dd, %d maps, %lu localLength ",
                1, l->map_count, lv->localLength);


    o += sprintf(s+o, ")");
    assert(o < 200);

    return s;
}

static
bool reuse_vector(Laik_Layout* l, int n, Laik_Layout* old, int nold)
{
    Laik_Layout_Vector *lv_new = laik_is_layout_vector(l);
    assert(lv_new);
    Laik_Layout_Vector *lv_old = laik_is_layout_vector(old);
    assert(lv_old);
    assert((n >= 0) && (n < l->map_count));

    if (laik_log_begin(1)) {
        laik_log_append("reuse_vector: check reuse for map %d in %s",
                        n, describe_vector(l));
        laik_log_flush(" using map %d in old %s", nold, describe_vector(old));
    }

   
    if (!(lv_new->localLength <= lv_old->localLength)) {
        // no, cannot reuse
        return false;
    }
    laik_log(1, "reuse_vector: old map %d can be reused (length %lu -> %lu)",
             nold,
             lv_new->localLength,
             lv_old->localLength);

    l->count = old->count; // TODO check if this is correct
    return true;
}


// copy vector
void copy_vector(Laik_Range* range,
                          Laik_Mapping* from, Laik_Mapping* to)
{
    printf("Calling copy_vector\n");
    Laik_Layout_Vector* fromLayout = laik_is_layout_vector(from->layout);
    Laik_Layout_Vector* toLayout = laik_is_layout_vector(to->layout);
    assert(fromLayout!=0);
    assert(toLayout != 0);
    unsigned int elemsize = from->data->elemsize;
    assert(elemsize == to->data->elemsize);

    if (laik_log_begin(1)) {
        laik_log_append("copy_vector of range ");
        laik_log_Range(range);
        laik_log_append(" (count %llu, elemsize %d) from mapping %p",
            laik_range_size(range), elemsize, from->start);
        laik_log_append(" (data '%s'/%d, %s) ",
            from->data->name, from->mapNo,
            describe_vector(from->layout));
        laik_log_flush("to mapping %p (data '%s'/%d, layout %s): ",
            to->start, to->data->name, to->mapNo,
             describe_vector(to->layout));
    }

    Laik_Index idx = range->from;
    uint64_t count = 0;
    do {
        int64_t fromOffset = offset_vector(from->layout, from->layoutSection, &idx);
        int64_t toOffset = offset_vector(to->layout, to->layoutSection, &idx);
        void* fromPtr = from->start + fromOffset * elemsize;
        void* toPtr = to->start + toOffset * elemsize;
#if 0
        if (laik_log_begin(1)) {
            laik_log_append(" copy idx ");
            laik_log_Index(range->space->dims, &idx);
            laik_log_flush(" from off %lu (ptr %p) to %lu (%p)",
                           fromOffset, fromPtr, toOffset, toPtr);
        }
#endif
        memcpy(toPtr, fromPtr, elemsize);
        count++;
    } while(next_idx(range, &idx));
    assert(count == laik_range_size(range));
}

// generic pack just using offset function from layout interface.
// this packs data according to lexicographical traversal
// return number of elements packed into provided buffer
unsigned int pack_vector(Laik_Mapping* m, Laik_Range* range,
                                  Laik_Index* idx, char* buf, unsigned int size)
{
    printf("Calling pack_vector\n");
    unsigned int elemsize = m->data->elemsize;
    Laik_Layout* layout = m->layout;
    int dims = m->layout->dims;

    if (laik_index_isEqual(dims, idx, &(range->to))) {
        // nothing left to pack
        return 0;
    }

    // range to pack must within local valid range of mapping
    assert(laik_range_within_range(range, &(m->requiredRange)));

    if (laik_log_begin(1)) {
        laik_log_append("        generic packing of range ");
        laik_log_Range(range);
        laik_log_append(" (count %llu, elemsize %d) from mapping %p",
            laik_range_size(range), elemsize, m->start);
        laik_log_append(" (data '%s'/%d, %s) at idx ",
            m->data->name, m->mapNo, layout->describe(layout));
        laik_log_Index(dims, idx);
        laik_log_flush(" into buf (size %d)", size);
    }

    unsigned int count = 0;
    while(size >= elemsize) {
        int64_t off = layout->offset(layout, m->layoutSection, idx);
        void* idxPtr = m->start + off * elemsize;
#if 0
        if (laik_log_begin(1)) {
            laik_log_append(" idx ");
            laik_log_Index(dims, &idx);
            laik_log_flush(": off %lu (ptr %p), left %d", off, ptr, size);
        }
#endif
        // copy element into buffer
        memcpy(buf, idxPtr, elemsize);
        size -= elemsize;
        buf += elemsize;
        count++;

        if (!next_idx(range, idx)) {
            *idx = range->to;
            break;
        }
    }

    if (laik_log_begin(1)) {
        laik_log_append("        packed '%s': end (", m->data->name);
        laik_log_Index(dims, idx);
        laik_log_flush("), %lu elems = %lu bytes, %d left",
                       count, count * elemsize, size);
    }

    return count;
}

// generic unpack just using offset function from layout interface
// this expects provided data to be packed according to lexicographical traversal
// return number of elements unpacked from provided buffer
unsigned int unpack_vector(Laik_Mapping* m, Laik_Range* range,
                                    Laik_Index* idx, char* buf, unsigned int size)
{
    printf("Calling unpack_vector\n");
    unsigned int elemsize = m->data->elemsize;
    Laik_Layout* layout = m->layout;
    int dims = m->layout->dims;

    // there should be something to unpack
    assert(size > 0);
    assert(!laik_index_isEqual(dims, idx, &(range->to)));

    // range to unpack into must be within local valid range of mapping
    assert(laik_range_within_range(range, &(m->requiredRange)));

    if (laik_log_begin(1)) {
        laik_log_append("        generic unpacking of range ");
        laik_log_Range(range);
        laik_log_append(" (count %llu, elemsize %d) into mapping %p",
            laik_range_size(range), elemsize, m->start);
        laik_log_append(" (data '%s'/%d, %s) at idx ",
            m->data->name, m->mapNo, layout->describe(layout));
        laik_log_Index(dims, idx);
        laik_log_flush(" from buf (size %d)", size);
    }

    unsigned int count = 0;
    while(size >= elemsize) {
        int64_t off = layout->offset(layout, m->layoutSection, idx);
        void* idxPtr = m->start + off * elemsize;
#if 0
        if (laik_log_begin(1)) {
            laik_log_append(" idx ");
            laik_log_Index(dims, &idx);
            laik_log_flush(": off %lu (ptr %p), left %d", off, ptr, size);
        }
#endif
        // copy element from buffer into mapping
        memcpy(idxPtr, buf, elemsize);
        size -= elemsize;
        buf += elemsize;
        count++;

        if (!next_idx(range, idx)) {
            *idx = range->to;
            break;
        }
    }

    if (laik_log_begin(1)) {
        laik_log_append("        unpacked '%s': end (", m->data->name);
        laik_log_Index(dims, idx);
        laik_log_flush("), %lu elems = %lu bytes, %d left",
                       count, count * elemsize, size);
    }

    return count;
}


// create layout for compact vector layout covering 1 range
Laik_Layout *laik_new_layout_vector(int n, Laik_Range *range, void *layout_data)
{   
    assert(n == 1); // This layout supports only 1 mapping
    int dims = range->space->dims;
    assert(dims == 1); // This layout supports only 1d spaces

    Laik_Layout_Vector* lv = malloc(sizeof(Laik_Layout_Vector));
    if (!lv) {
        laik_panic("Out of memory allocating Laik_Layout_Vector object");
        exit(1); // not actually needed, laik_panic never returns
    }
    // count calculated later
    laik_init_layout(&(lv->h), dims, n, laik_range_size(range),
                     section_vector,
                     mapno_vector,
                     offset_vector,
                     reuse_vector,
                     describe_vector,
                     pack_vector,
                     unpack_vector,
                     copy_vector);

    return (Laik_Layout*) lv;
}


uint64_t laik_get_length_vector(Laik_Layout *l)
{
    Laik_Layout_Vector* lv = laik_is_layout_vector(l);
    assert(lv != 0);

    return lv->localLength;
}
