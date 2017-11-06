/* This file is part of the LAIK parallel container library.
 * Copyright (c) 2017 Josef Weidendorfer
 *
 * LAIK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, version 3.
 *
 * LAIK is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Distributed Markov chain example.
 */

#ifdef USE_MPI
#include "laik-backend-mpi.h"
#else
#include "laik-backend-single.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct _MGraph {
    int n;   // number of states
    int in;  // fan-in
    int* cm; // connectivity
    double* pm; // probabilities
} MGraph;

// Produce a graph with <n> nodes and some arbitrary connectivity
// with a fan-in <in>. The resulting graph will be stored in
// <cm>[i,c], which is a <n> * (<in> +1) matrix storing the incoming nodes
// of node i in row i, using columns 1 .. <in> (column 0 is set to i).
// <pm>[i,j] is initialized with the probability of the transition
// from node <cm>[i,j] to node i, with cm[i,0] the prob for staying.
void init(MGraph* mg)
{
    int n = mg->n;
    int in = mg->in;
    int* cm = mg->cm;
    double* pm = mg->pm;

    // for normalization of probabilites
    double* sum = malloc(n * sizeof(double));
    for(int i=0; i < n; i++) sum[i] = 0.0;

    // some kind of ring structure
    for(int i=0; i < n; i++) {
        int step = 1;
        cm[i * (in + 1) + 0] = i; // stay in i
        pm[i * (in + 1) + 0] = 5;
        sum[i] += 5;
        for(int j = 1; j <= in; j++) {
            int fromNode = (i + step) % n;
            double prob = (double) ((j+i) % (5 * in)) + 1;
            sum[fromNode] += prob;
            cm[i * (in + 1) + j] = fromNode;
            pm[i * (in + 1) + j] = prob;
            step = 2 * step + j;
            while(step > n) step -= n;
        }
    }
    // normalization. never should do divide-by-0
    for(int i=0; i < n; i++) {
        for(int j = 0; j <= in; j++)
            pm[i * (in + 1) + j] /= sum[ cm[i * (in + 1) + j] ];
    }
}

void print(MGraph* mg)
{
    int n = mg->n;
    int in = mg->in;
    int* cm = mg->cm;
    double* pm = mg->pm;

    for(int i = 0; i < n; i++) {
        printf("State %2d: stay %.3f ", i, pm[i * (in + 1)]);
        for(int j = 1; j <= in; j++)
            printf("<=(%.3f)=%-2d  ",
                   pm[i * (in + 1) + j], cm[i * (in + 1) + j]);
        printf("\n");
    }
}


void run_markovPartitioner(Laik_Partitioner* pr,
                           Laik_BorderArray* ba, Laik_BorderArray* otherBA)
{
    MGraph* mg = laik_partitioner_data(pr);
    int in = mg->in;
    int* cm = mg->cm;

    Laik_Slice slc;

    // go over states and add itself and incoming states to new partitioning
    int sliceCount = laik_borderarray_getcount(otherBA);
    for(int i = 0; i < sliceCount; i++) {
        Laik_TaskSlice* ts = laik_borderarray_get_tslice(otherBA, i);
        const Laik_Slice* s = laik_taskslice_getslice(ts);
        int task = laik_taskslice_gettask(ts);
        for(int st = s->from.i[0]; st < s->to.i[0]; st++) {
            int off = st * (in + 1);
            // j=0: state itself
            for(int j = 0; j <= in; j++) {
                int inNode = cm[off + j];
                slc.from.i[0] = inNode;
                slc.to.i[0] = inNode + 1;
                laik_append_slice(ba, task, &slc, 1, 0);
            }
        }
    }
}



// iteratively calculate probability distribution, return last written data
Laik_Data* run(MGraph* mg, int miter,
               Laik_Data* data1, Laik_Data* data2,
               Laik_Partitioning* pWrite, Laik_Partitioning* pRead)
{
    int n = mg->n;
    int in = mg->in;
    int* cm = mg->cm;
    double* pm = mg->pm;

    // start reading from data1, writing to data2
    Laik_Data *dRead = data1, *dWrite = data2;
    double *src, *dst;
    uint64_t srcCount, dstCount, srcFrom, dstFrom, dstTo;

    int iter = 0;
    while(1) {
        // switch dRead to pRead, dWrite to pWrite
        laik_switchto(dRead,  pRead,  LAIK_DF_CopyIn);
        laik_map_def1(dRead, (void**) &src, &srcCount);
        srcFrom = laik_local2global_1d(dRead, 0);

        laik_switchto(dWrite, pWrite, LAIK_DF_CopyOut);
        laik_map_def1(dWrite, (void**) &dst, &dstCount);
        laik_my_slice_1d(pWrite, 0, &dstFrom, &dstTo);
        assert(dstCount == dstTo - dstFrom);

        // spread values according to probability distribution
        for(int i = dstFrom; i < dstTo; i++) {
            int off = i * (in + 1);
            double v = src[i - srcFrom] * pm[off];
            for(int j = 1; j <= in; j++)
                v += src[cm[off + j] - srcFrom] * pm[off + j];
            dst[i - dstFrom] = v;
        }

        iter++;
        if (iter == miter) break;

        // swap role of data1 and data2
        if (dRead == data1) { dRead = data2; dWrite = data1; }
        else                { dRead = data1; dWrite = data2; }
    }

    return dWrite;
}

int main(int argc, char* argv[])
{
#ifdef USE_MPI
    Laik_Instance* inst = laik_init_mpi(&argc, &argv);
#else
    Laik_Instance* inst = laik_init_single();
#endif
    Laik_Group* world = laik_world(inst);

    int n = 1000000;
    int in = 10;
    int miter = 10;
    int doPrint = 0;

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) in = atoi(argv[2]);
    if (argc > 3) miter = atoi(argv[3]);
    if (argc > 4) doPrint = 1;
    if (n == 0) n = 1000000;
    if (in == 0) in = 10;

    if (laik_myid(world) == 0) {
        printf("Init Markov chain with %d states, max fan-in %d\n", n, in);
        printf("Run %d iterations each.\n", miter);
    }

    MGraph mg;
    mg.n = n;
    mg.in = in;
    mg.cm = malloc(n * (in + 1) * sizeof(int));
    mg.pm = malloc(n * (in + 1) * sizeof(double));

    init(&mg);
    if (doPrint) print(&mg);

    // two 1d arrays, using same space
    Laik_Space* space = laik_new_space_1d(inst, n);
    Laik_Data* data1 = laik_new_data(world, space, laik_Double);
    Laik_Data* data2 = laik_new_data(world, space, laik_Double);

    // partitionings:
    // - pWrite: distributes the state values to update
    // - pRead : provides access to values of incoming states
    // - pMaster: all data at master, for checksum
    // pWrite/pRead partitionings are assigned to either data1/data2,
    // exchanged after every iteration
    Laik_Partitioning *pWrite, *pRead, *pMaster;
    Laik_Partitioner* pr;
    pWrite = laik_new_partitioning(world, space,
                                   laik_new_block_partitioner1(), 0);
    pr = laik_new_partitioner("markovin", run_markovPartitioner, &mg);
    pRead = laik_new_partitioning(world, space, pr, pWrite);
    pMaster = laik_new_partitioning(world, space, laik_Master, 0);

    if (laik_myid(world) == 0)
        printf("Start with state 0 prob 1 ...\n");

    // distributed initialization of data1
    double *v;
    uint64_t count, off;
    laik_switchto(data1, pWrite, LAIK_DF_CopyOut);
    laik_map_def1(data1, (void**) &v, &count);
    for(uint64_t i = 0; i < count; i++)
        v[i] = (double) 0.0;
    // set state 0 to probability 1
    if (laik_global2local_1d(data1, 0, &off)) {
        // if global index 0 is local, it must be at local index 0
        assert(off == 0);
        v[off] = 1.0;
    }

    Laik_Data* dRes;
    dRes = run(&mg, miter, data1, data2, pWrite, pRead);

    laik_switchto(dRes, pMaster, LAIK_DF_CopyIn);
    laik_map_def1(dRes, (void**) &v, &count);
    if (laik_myid(world) == 0) {
        assert(count == n);
        double sum = 0.0;
        for(int i=0; i < n; i++)
            sum += v[i];
        printf("  result probs: p0 = %g, p1 = %g, p2 = %g, Sum: %f\n",
               v[0], v[1], v[2], sum);
    }

    if (laik_myid(world) == 0)
        printf("Start with state 1 prob 1 ...\n");
    laik_switchto(data1, pWrite, LAIK_DF_CopyOut);
    laik_map_def1(data1, (void**) &v, &count);
    for(uint64_t i = 0; i < count; i++)
        v[i] = (double) 0.0;
    if (laik_global2local_1d(data1, 1, &off))
        v[off] = 1.0;
    dRes = run(&mg, miter, data1, data2, pWrite, pRead);
    laik_switchto(dRes, pMaster, LAIK_DF_CopyIn);
    laik_map_def1(dRes, (void**) &v, &count);
    if (laik_myid(world) == 0) {
        double sum = 0.0;
        for(int i=0; i < n; i++)
            sum += v[i];
        printf("  result probs: p0 = %g, p1 = %g, p2 = %g, Sum: %f\n",
               v[0], v[1], v[2], sum);
    }

    if (laik_myid(world) == 0)
        printf("Start with all probs equal ...\n");
    laik_switchto(data1, pWrite, LAIK_DF_CopyOut);
    laik_map_def1(data1, (void**) &v, &count);
    double p = 1.0 / n;
    for(uint64_t i = 0; i < count; i++)
        v[i] = p;
    dRes = run(&mg, miter, data1, data2, pWrite, pRead);
    laik_switchto(dRes, pMaster, LAIK_DF_CopyIn);
    laik_map_def1(dRes, (void**) &v, &count);
    if (laik_myid(world) == 0) {
        double sum = 0.0;
        for(int i=0; i < n; i++)
            sum += v[i];
        printf("  result probs: p0 = %g, p1 = %g, p2 = %g, Sum: %f\n",
               v[0], v[1], v[2], sum);
    }

    laik_finalize(inst);
    return 0;
}