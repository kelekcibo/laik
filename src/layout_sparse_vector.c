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

// For calculating mapping for the sparse vector layout
typedef struct _Laik_Interval_Vector Laik_Interval_Vector;
struct _Laik_Interval_Vector
{
    int64_t from;
    int64_t to; // including "from", but excluding "to"
};

typedef struct _Laik_Map_Vector Laik_Map_Vector;
struct _Laik_Map_Vector
{
    uint64_t size;                      // number of intervals stored in the map
    Laik_Interval_Vector * intervals; // All intervals are stored here

    int64_t lower_bound; // including
    int64_t upper_bound; // excluding
};

typedef struct _Laik_Layout_Vector Laik_Layout_Vector;
struct _Laik_Layout_Vector {
    Laik_Layout h;
    int id;                                   // For debug purposes
    int64_t localLength;                      // length of sparse vector without external values  
    uint64_t numberOfExternalValues;          // external values stored after localLength elements
    uint64_t currentExternalVallue;           // Mapping for external indexes to allocation indexes (offset in the allocated buffer by LAIK)
    Laik_Map_Vector * globalToLocalMap;       // Mapping needed from own global indexes to allocation indexes (offset in the allocated buffer by LAIK)
    uint64_t allocatedRangeCount;          // quick solution for implementing reuse_vector
};


//--------------------------------------------------------------
// interface implementation of sparse vector layout
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

    // Check if idx is in map // only local idx handled here and not properly... need to go trough intervals as it has gaps
    if (idx->i[0] <= lv->globalToLocalMap->upper_bound && idx->i[0] >= lv->globalToLocalMap->lower_bound)
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
    assert(n == 0);
    Laik_Layout_Vector *lv = laik_is_layout_vector(l);
    assert(lv);
    Laik_Map_Vector * m = lv->globalToLocalMap;
    assert(m);

    int64_t idx_val = idx->i[0];
    Laik_Interval_Vector * intervals = m->intervals;
    int64_t localOffset = 0;
    bool wasInLocalRange = false;
    for (uint64_t i = 0; i < m->size; i++)
    {
        if (idx_val >= intervals[i].from && idx_val < intervals[i].to)
        {
            localOffset += idx_val - intervals[i].from;
            wasInLocalRange = true;
            break; 
        } 
        else if (idx_val < intervals[i].from)
        {
            // we catched an idx which is not locally owned by this proc
            // break so we do not iterate over all intervals unnecessary
            break;
        }
        localOffset += intervals[i].to - intervals[i].from;
    }

    // we calculated the offset for locally owned global index
    if(wasInLocalRange)
    {
        assert(localOffset >= 0 && localOffset < lv->localLength);
        return localOffset;
    }

    // if a Laik_Data container does not receive external values at all, there went something wrong
    if(lv->numberOfExternalValues == 0 && !wasInLocalRange)
        laik_panic("A vector with no external values should not map a global index of an external value.");

    // if <idx> is an external value, we need to calculate the correct offset
    // we have following layout: | Local values | external values |
    // Thus, we copy external values to the end of local values

    localOffset = lv->localLength; // offset starts there, if <idx> was not a locally owned idx

    if(lv->currentExternalVallue == lv->numberOfExternalValues)
        lv->currentExternalVallue = 0;

    return localOffset + lv->currentExternalVallue++;
}


static
char* describe_vector(Laik_Layout* l)
{
    static char s[200];

    Laik_Layout_Vector *lv = laik_is_layout_vector(l);
    assert(lv);

    int o;
    o = sprintf(s, "sparse vector (%dd, %d map, %lu localLength, %lu numberOfExternalValue in external partitioning, %ld count , %lu allocatedRangeCount",
                l->dims, l->map_count, lv->localLength, lv->numberOfExternalValues, l->count, lv->allocatedRangeCount);

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
    // printf("LAIK %d\tCalling reuse_vector\n", lv_new->id);

    if (laik_log_begin(1)) {
        laik_log_append("reuse_vector: check reuse for map %d in %s",
                        n, describe_vector(l));
        laik_log_flush(" using map %d in old %s", nold, describe_vector(old));
    }

    bool new_totalSize_fits = lv_new->allocatedRangeCount <= lv_old->allocatedRangeCount; // do not reuse a vector, if memory of new layout does not fit into old
    bool vector_size_changed = lv_new->localLength != lv_old->localLength;  // do not reuse a vector, if localLength changed (due to repartitioning)
    if (!new_totalSize_fits || vector_size_changed)
    {
        // if vector size is the same, then this means we switched from a local to external partitioning.
        // this means, that we need the Map from the local partitioning as it is not calculated for external partitioning
        if(!(vector_size_changed))
            lv_new->globalToLocalMap = lv_old->globalToLocalMap; // when we repartition, we are not able to make the optimisation of switching to the external partitioning first

        laik_log(1, "reuse_vector: old map %d cannot be reused (allocatedRangeCount %lu -> %lu \t localLength %lu -> %lu)\n",
                 nold,
                 lv_new->allocatedRangeCount,
                 lv_old->allocatedRangeCount,
                 lv_new->localLength,
                 lv_old->localLength);
        return false;
    }

    laik_log(1, "reuse_vector: old map %d can be reused (length/count %lu/%lu -> %lu/%lu)",
             nold,
             lv_new->localLength,
             lv_new->h.count,
             lv_old->localLength,
             lv_old->h.count);

    // l->count = old->count; // TODO check if this is correct

    // if new layout is layout for external partitioning, we do not calculate mapping, thus get mapping from local partitioning (old layout)
    if((int64_t)l->count != lv_new->localLength)
        lv_new->globalToLocalMap = lv_old->globalToLocalMap;

    // update allocatedRangeCount as mentioned in laik_new_vector_layout()
    lv_new->allocatedRangeCount = lv_old->allocatedRangeCount;  

    return true;
}

// copy vector
void copy_vector(Laik_Range* range,
                          Laik_Mapping* from, Laik_Mapping* to)
{
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
        laik_log_append("        vector packing of range ");
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
        laik_log_append("        vector unpacking of range ");
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

/**
 * @brief LAIK will automatically calculate the mapping for a vector.
 *
 * The mapping will be explained in regard to the HPCG Benchmark.
 * We have a local partitioning and external partitioning and we are switching between them to communicate values.
 * In order to eliminate the gaps in the memory, we provide this layout. This means we need to calculate the mapping of global/external indexes to local indexes.
 * 
 * Approach:
 * Laik_Map_Vector map contains all global ranges/intervals, which are assigned to the process in the global problem. 
 * In context of HPCG, it means the global rows assigned to process <i>, as process <i> updates the global vector at index <i> after each iteration.
 * 
 * When we have all global intervals, we can calculate the exact offset in the allocated buffer @see offset_vector().
 * 
 * We only calculate the mapping when the local partitioning is active. The external partitioning gets the mapping from the layout for the local partitioning in reuse_vector().
 *
 * @param l Sparse Vector Layout
 * @param list All ranges this process owns in current active partitioning
 * @param map_size number of intervals the mapping will have
 * @param myid rank id
 */
void calculate_mapping(Laik_Layout *l, Laik_RangeList *list, uint64_t map_size, int myid)
{
    Laik_Layout_Vector *lv = laik_is_layout_vector(l);
    assert(lv);
    assert(map_size != 0);

    Laik_Map_Vector *m = (Laik_Map_Vector *) malloc(sizeof(Laik_Map_Vector));
    m->size = map_size;

    // Allocate memory for <map_size> chunks
    m->intervals = (Laik_Interval_Vector *) malloc(map_size * sizeof(Laik_Interval_Vector));

    uint64_t chunksInitialised = 0; // for testing correctness (minimum size is 1)
    int mapNo = 0;                  // Sparse vector layout supports only one mapping
    Laik_Range* current_range;
    for (unsigned int o = list->off[myid]; o < list->off[myid + 1]; o++)
    {
        // range covering all task ranges for the map 0
        // set lower bound of map containing intervals (it is in the first range, as ranges are sorted in ascending order)
        m->lower_bound = list->trange[o].range.from.i[0];

        // For calculating the number of needed map entries
        Laik_Range *startRangeOfInterval = &(list->trange[o].range); // store beginning range of a new interval
        bool initialiseInterval = false;
        unsigned int off = list->off[myid + 1];
        while ((o + 1 < off) && (list->trange[o + 1].mapNo == mapNo))
        {
            o++;
            current_range = &(list->trange[o].range);

            // check that current range is not neighbour of previous range
            if (list->trange[o - 1].range.to.i[0] != current_range->from.i[0])
                initialiseInterval = true; // if yes, initialise chunk/interval

            // check that we are not in the last iteration
            if(!((o + 1 < off)) && !initialiseInterval)
            {
                // if yes, we need to initialize interval[].to differently
                m->intervals[chunksInitialised].from = startRangeOfInterval->from.i[0];
                m->intervals[chunksInitialised].to = current_range->to.i[0];
                chunksInitialised++;
                break;
            }

            if(initialiseInterval)
            {
                initialiseInterval = false;
                m->intervals[chunksInitialised].from = startRangeOfInterval->from.i[0];
                m->intervals[chunksInitialised].to = list->trange[o - 1].range.to.i[0];
                // start new interval
                startRangeOfInterval = current_range;
                // increase counter for chunks
                chunksInitialised++;

                // check that we are not in the last iteration
                if (!(o + 1 < off))
                {
                    // if yes, we need to initialize interval[].to differently, if
                    // current range is not neighbour of previous range
                    m->intervals[chunksInitialised].from = current_range->from.i[0];
                    m->intervals[chunksInitialised].to = current_range->to.i[0];
                    // increase counter for chunks
                    chunksInitialised++;
                }
            }
        }
    }
    assert(chunksInitialised == map_size);

    m->upper_bound = current_range->to.i[0];
    lv->globalToLocalMap = m;

    return;
}

// create layout for compact vector layout covering 1 range
Laik_Layout *laik_new_layout_vector(int n, Laik_Range *range, void * layout_data)
{
    assert(n == 1); // This layout supports only 1 mapping
    assert(layout_data != 0); // Need layout data
    int dims = range->space->dims;
    assert(dims == 1); // This layout supports only 1d spaces

    Laik_Layout_Vector* lv = malloc(sizeof(Laik_Layout_Vector));
    if (!lv) {
        laik_panic("Out of memory allocating Laik_Layout_Vector object");
        exit(1); // not actually needed, laik_panic never returns
    }

    Laik_vector_data * vd = (Laik_vector_data *) layout_data;

    lv->localLength = vd->localLength;
    lv->numberOfExternalValues = vd->numberOfExternalValues;
    lv->id = vd->id;
    lv->currentExternalVallue = 0;
    lv->allocatedRangeCount = laik_range_size(range); // if we reuse, we need to update this value to the greater value

    laik_init_layout(&(lv->h), dims, n, laik_range_size(range),
                     section_vector,
                     mapno_vector,
                     offset_vector,
                     reuse_vector,
                     describe_vector,
                     pack_vector,
                     unpack_vector,
                     copy_vector);

    laik_log(1, "laik_new_layout_vector: New sparse vector layout has been created: %s)", describe_vector(&(lv->h)));
    return (Laik_Layout*) lv;
}


uint64_t laik_get_length_vector(Laik_Layout *l)
{
    Laik_Layout_Vector* lv = laik_is_layout_vector(l);
    if(!lv)
        return 0;

    return lv->localLength;
}

uint64_t laik_get_numberOfExternalValues_vector(Laik_Layout *l)
{
    Laik_Layout_Vector *lv = laik_is_layout_vector(l);
    assert(lv != 0);

    return lv->numberOfExternalValues;
}

uint64_t laik_get_id_vector(Laik_Layout *l)
{
    Laik_Layout_Vector *lv = laik_is_layout_vector(l);
    assert(lv != 0);

    return lv->id;
}

// void laik_set_length_vector(Laik_Data *d, uint64_t localLength)
// {
//     Laik_Layout * l = laik_get_map(d, 0)->layout; // This layout makes use of only one mapping
//     Laik_Layout_Vector *lv = laik_is_layout_vector(l);
//     assert(lv != 0);

//     lv->localLength = localLength;
//     lv->numberOfExternalValues = l->count - localLength;
//     return;
// }

// DEBUG
void laik_print_local_Map(Laik_Data *d, int id)
{
    Laik_Layout_Vector *lv = laik_is_layout_vector(laik_get_map(d, 0)->layout);
    assert(lv != 0);

    if (lv->id != id)
        return;

    // print
    Laik_Map_Vector * m = lv->globalToLocalMap;
    if(!m)
    {
        printf("No Map active\n");
        return;
    }
    printf("######## LAIK %d Mapping of %s\n", lv->id, describe_vector(&(lv->h)));
    printf("chunks %lu\tupper bound %lu\tlower bound %lu\n", m->size, m->upper_bound, m->lower_bound);
    for (uint64_t i = 0; i < m->size; i++)
        printf("Interval %ld\t[%ld;%ld[\n", i, m->intervals[i].from, m->intervals[i].to);
    return;
}

void laik_print_local_Map2(Laik_Layout *l, int id)
{
    Laik_Layout_Vector *lv = laik_is_layout_vector(l);
    assert(lv != 0);

    if (lv->id != id)
        return;

    // print
    Laik_Map_Vector *m = lv->globalToLocalMap;
    if (!m)
    {
        printf("No Map active\n");
        return;
    }
    printf("######## LAIK %d Mapping of %s\n", lv->id, describe_vector(l));
    printf("chunks %lu\tupper bound %lu\tlower bound %lu\n", m->size, m->upper_bound, m->lower_bound);
    for (uint64_t i = 0; i < m->size; i++)
        printf("Interval %ld\t[%ld;%ld[\n", i, m->intervals[i].from, m->intervals[i].to);
    return;
}