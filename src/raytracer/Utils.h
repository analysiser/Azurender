//
//  Utils.h
//  P3
//
//  Created by Xiao Li on 5/3/14.
//
//

//#include <stdio.h>


struct cPhoton
{
	/* data */
    float position[3];              // 12
    unsigned int index;             // 4   // index in photon list
    int splitAxis;                  // 4     // char and short cannot be compiled by ispc
    int padding;
};

struct metaCPhoton
{
    unsigned int cphotonIndex;      // index in cphoton list
    unsigned int sortedIndex[3];    // index in sorted index lists
};

struct KDNode
{
    int cphotonIndex;
    int splitAxis;
    float splitValue;
    int head;           // including
    int tail;           // excluding
    int size;
    bool isLeaf;
    
    KDNode *left;
    KDNode *right;
};
