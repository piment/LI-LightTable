/* trees.c -- output deflated data using Huffman coding
 * Copyright (C) 1995-2010 Jean-loup Gailly
 * detect_data_type() function provided freely by Cosmin Truta, 2006
 * For conditions of distribution and use, see copyright notice in zlib.h
 */

/*
 *  ALGORITHM
 *
 *      The "deflation" process uses several Huffman trees. The more
 *      common source values are represented by shorter bit sequences.
 *
 *      Each code tree is stored in a compressed form which is itself
 * a Huffman encoding of the lengths of all the code strings (in
 * ascending order by source values).  The actual code strings are
 * reconstructed from the lengths in the inflate process, as described
 * in the deflate specification.
 *
 *  REFERENCES
 *
 *      Deutsch, L.P.,"'Deflate' Compressed Data Format Specification".
 *      Available in ftp.uu.net:/pub/archiving/zip/doc/deflate-1.1.doc
 *
 *      Storer, James A.
 *          Data Compression:  Methods and Theory, pp. 49-50.
 *          Computer Science Press, 1988.  ISBN 0-7167-8156-5.
 *
 *      Sedgewick, R.
 *          Algorithms, p290.
 *          Addison-Wesley, 1983. ISBN 0-201-06672-6.
 */

/* @(#) $Id$ */

/* #define GEN_TREES_H */

#include "deflate.h"

#ifdef DEBUG
#  include <ctype.h>
#endif

/* ===========================================================================
 * Constants
 */

#define MAX_BL_BITS 7
/* Bit length codes must not exceed MAX_BL_BITS bits */

#define END_BLOCK 256
/* end of block literal code */

#define REP_3_6      16
/* repeat previous bit length 3-6 times (2 bits of repeat count) */

#define REPZ_3_10    17
/* repeat a zero length 3-10 times  (3 bits of repeat count) */

#define REPZ_11_138  18
/* repeat a zero length 11-138 times  (7 bits of repeat count) */

local const int extra_lbits[LENGTH_CODES] /* extra bits for each length code */
   = {0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0};

local const int extra_dbits[D_CODES] /* extra bits for each distance code */
   = {0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13};

local const int extra_blbits[BL_CODES]/* extra bits for each bit length code */
   = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,3,7};

local const uch bl_order[BL_CODES]
   = {16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15};
/* The lengths of the bit length codes are sent in order of decreasing
 * probability, to avoid transmitting the lengths for unused bit length codes.
 */

#define Buf_size (8 * 2*sizeof(char))
/* Number of bits used within bi_buf. (bi_buf might be implemented on
 * more than 16 bits on some systems.)
 */

/* ===========================================================================
 * Local data. These are initialized only once.
 */

#define DIST_CODE_LEN  512 /* see definition of array dist_code below */

#if defined(GEN_TREES_H) || !defined(STDC)
/* non ANSI compilers may not accept trees.h */

local ct_data static_ltree[L_CODES+2];
/* The static literal tree. Since the bit lengths are imposed, there is no
 * need for the L_CODES extra codes used during heap construction. However
 * The codes 286 and 287 are needed to build a canonical tree (see _tr_init
 * below).
 */

local ct_data static_dtree[D_CODES];
/* The static distance tree. (Actually a trivial tree since all codes use
 * 5 bits.)
 */

uch _dist_code[DIST_CODE_LEN];
/* Distance codes. The first 256 values correspond to the distances
 * 3 .. 258, the last 256 values correspond to the top 8 bits of
 * the 15 bit distances.
 */

uch _length_code[MAX_MATCH-MIN_MATCH+1];
/* length code for each normalized match length (0 == MIN_MATCH) */

local int base_length[LENGTH_CODES];
/* First normalized length for each code (0 = MIN_MATCH) */

local int base_dist[D_CODES];
/* First normalized distance for each code (0 = distance of 1) */

#else
#  include "trees.h"
#endif /* GEN_TREES_H */

struct static_tree_desc_s {
    const ct_data *static_tree;  /* static tree or NULL */
    const intf *extra_bits;      /* extra bits for each code or NULL */
    int     extra_base;          /* base index for extra_bits */
    int     elems;               /* max number of elements in the tree */
    int     max_length;          /* max bit length for the codes */
};

local static_tree_desc  static_l_desc =
{static_ltree, extra_lbits, LITERALS+1, L_CODES, MAX_BITS};

local static_tree_desc  static_d_desc =
{static_dtree, extra_dbits, 0,          D_CODES, MAX_BITS};

local static_tree_desc  static_bl_desc =
{(const ct_data *)0, extra_blbits, 0,   BL_CODES, MAX_BL_BITS};

/* ===========================================================================
 * Local (static) routines in this file.
 */

local void tr_static_init OF((void));
local void init_block     OF((deflate_state *s));
local void pqdownheap     OF((deflate_state *s, ct_data *tree, int k));
local void gen_bitlen     OF((deflate_state *s, tree_desc *desc));
local void gen_codes      OF((ct_data *tree, int max_code, ushf *bl_count));
local void build_tree     OF((deflate_state *s, tree_desc *desc));
local void scan_tree      OF((deflate_state *s, ct_data *tree, int max_code));
local void send_tree      OF((deflate_state *s, ct_data *tree, int max_code));
local int  build_bl_tree  OF((deflate_state *s));
local void send_all_trees OF((deflate_state *s, int lcodes, int dcodes,
                              int blcodes));
local void compress_block OF((deflate_state *s, ct_data *ltree,
                              ct_data *dtree));
local int  detect_data_type OF((deflate_state *s));
local unsigned bi_reverse OF((unsigned value, int length));
local void bi_windup      OF((deflate_state *s));
local void bi_flush       OF((deflate_state *s));
local void copy_block     OF((deflate_state *s, charf *buf, unsigned len,
                              int header));

#ifdef GEN_TREES_H
local void gen_trees_header OF((void));
#endif

#ifndef DEBUG
#  define send_code(s, c, tree) send_bits(s, tree[c].Code, tree[c].Len)
   /* Send a code of the given tree. c and tree must not have side effects */

#else /* DEBUG */
#  define send_code(s, c, tree) \
     { if (z_verbose>2) fprintf(stderr,"\ncd %3d ",(c)); \
       send_bits(s, tree[c].Code, tree[c].Len); }
#endif

/* ===========================================================================
 * Output a short LSB first on the stream.
 * IN assertion: there is enough room in pendingBuf.
 */
#define put_short(s, w) { \
    put_byte(s, (uch)((w) & 0xff)); \
    put_byte(s, (uch)((ush)(w) >> 8)); \
}

/* ===========================================================================
 * Send a value on a given number of bits.
 * IN assertion: length <= 16 and value fits in length bits.
 */
#ifdef DEBUG
local void send_bits      OF((deflate_state *s, int value, int length));

local void send_bits(s, value, length)
    deflate_state *s;
    int value;  /* value to send */
    int length; /* number of bits */
{
    Tracevv((stderr," l %2d v %4x ", length, value));
    Assert(length > 0 && length <= 15, "invalid length");
    s->bits_sent += (ulg)length;

    /* If not enough room in bi_buf, use (valid) bits from bi_buf and
     * (16 - bi_valid) bits from value, leaving (width - (16-bi_valid))
     * unused bits in value.
     */
    if (s->bi_valid > (int)Buf_size - length) {
        s->bi_buf |= (ush)value << s->bi_valid;
        put_short(s, s->bi_buf);
        s->bi_buf = (ush)value >> (Buf_size - s->bi_valid);
        s->bi_valid += length - Buf_size;
    } else {
        s->bi_buf |= (ush)value << s->bi_valid;
        s->bi_valid += length;
    }
}
#else /* !DEBUG */

#define send_bits(s, value, length) \
{ int len = length;\
  if (s->bi_valid > (int)Buf_size - len) {\
    int val = value;\
    s->bi_buf |= (ush)val << s->bi_valid;\
    put_short(s, s->bi_buf);\
    s->bi_buf = (ush)val >> (Buf_size - s->bi_valid);\
    s->bi_valid += len - Buf_size;\
  } else {\
    s->bi_buf |= (ush)(value) << s->bi_valid;\
    s->bi_valid += len;\
  }\
}
#endif /* DEBUG */


/* the arguments must not have side effects */

/* ===========================================================================
 * Initialize the various 'constant' tables.
 */
local void tr_static_init()
{
#if defined(GEN_TREES_H) || !defined(STDC)
    static int static_init_done = 0;
    int n;        /* iterates over tree elements */
    int bits;     /* bit counter */
    int length;   /* length value */
    int code;     /* code value */
    int dist;     /* distance index */
    ush bl_count[MAX_BITS+1];
    /* number of codes at each bit length for an optimal tree */

    if (static_init_done) return;

    /* For some embedded targets, global variables are not initialized: */
#ifdef NO_INIT_GLOBAL_POINTERS
    static_l_desc.static_tree = static_ltree;
    static_l_desc.extra_bits = extra_lbits;
    static_d_desc.static_tree = static_dtree;
    static_d_desc.extra_bits = extra_dbits;
    static_bl_desc.extra_bits = extra_blbits;
#endif

    /* Initialize the mapping length (0..255) -> length code (0..28) */
    length = 0;
    for (code = 0; code < LENGTH_CODES-1; code++) {
        base_length[code] = length;
        for (n = 0; n < (1<<extra_lbits[code]); n++) {
            _length_code[length++] = (uch)code;
        }
    }
    Assert (length == 256, "tr_static_init: length != 256");
    /* Note that the length 255 (match length 258) can be represented
     * in two different ways: code 284 + 5 bits or code 285, so we
     * overwrite length_code[255] to use the best encoding:
     */
    _length_code[length-1] = (uch)code;

    /* Initialize the mapping dist (0..32K) -> dist code (0..29) */
    dist = 0;
    for (code = 0 ; code < 16; code++) {
        base_dist[code] = dist;
        for (n = 0; n < (1<<extra_dbits[code]); n++) {
            _dist_code[dist++] = (uch)code;
        }
    }
    Assert (dist == 256, "tr_static_init: dist != 256");
    dist >>= 7; /* from now on, all distances are divided by 128 */
    for ( ; code < D_CODES; code++) {
        base_dist[code] = dist << 7;
        for (n = 0; n < (1<<(extra_dbits[code]-7)); n++) {
            _dist_code[256 + dist++] = (uch)code;
        }
    }
    Assert (dist == 256, "tr_static_init: 256+dist != 512");

    /* Construct the codes of the static literal tree */
    for (bits = 0; bits <= MAX_BITS; bits++) bl_count[bits] = 0;
    n = 0;
    while (n <= 143) static_ltree[n++].Len = 8, bl_count[8]++;
    while (n <= 255) static_ltree[n++].Len = 9, bl_count[9]++;
    while (n <= 279) static_ltree[n++].Len = 7, bl_count[7]++;
    while (n <= 287) static_ltree[n++].Len = 8, bl_count[8]++;
    /* Codes 286 and 287 do not exist, but we must include them in the
     * tree construction to get a canonical Huffman tree (longest code
     * all ones)
     */
    gen_codes((ct_data *)static_ltree, L_CODES+1, bl_count);

    /* The static distance tree is trivial: */
    for (n = 0; n < D_CODES; n++) {
        static_dtree[n].Len = 5;
        static_dtree[n].Code = bi_reverse((unsigned)n, 5);
    }
    static_init_done = 1;

#  ifdef GEN_TREES_H
    gen_trees_header();
#  endif
#endif /* defined(GEN_TREES_H) || !defined(STDC) */
}

/* ===========================================================================
 * Genererate the file trees.h describing the static trees.
 */
#ifdef GEN_TREES_H
#  ifndef DEBUG
#    include <stdio.h>
#  endif

#  define SEPARATOR(i, last, width) \
      ((i) == (last)? "\n};\n\n" :    \
       ((i) % (width) == (width)-1 ? ",\n" : ", "))

void gen_trees_header()
{
    FILE *header = fopen("trees.h", "w");
    int i;

    Assert (header != NULL, "Can't open trees.h");
    fprintf(header,
            "/* header created automatically with -DGEN_TREES_H */\n\n");

    fprintf(header, "local const ct_data static_ltree[L_CODES+2] = {\n");
    for (i = 0; i < L_CODES+2; i++) {
        fprintf(header, "{{%3u},{%3u}}%s", static_ltree[i].Code,
                static_ltree[i].Len, SEPARATOR(i, L_CODES+1, 5));
    }

    fprintf(header, "local const ct_data static_dtree[D_CODES] = {\n");
    for (i = 0; i < D_CODES; i++) {
        fprintf(header, "{{%2u},{%2u}}%s", static_dtree[i].Code,
                static_dtree[i].Len, SEPARATOR(i, D_CODES-1, 5));
    }

    fprintf(header, "const uch ZLIB_INTERNAL _dist_code[DIST_CODE_LEN] = {\n");
    for (i = 0; i < DIST_CODE_LEN; i++) {
        fprintf(header, "%2u%s", _dist_code[i],
                SEPARATOR(i, DIST_CODE_LEN-1, 20));
    }

    fprintf(header,
        "const uch ZLIB_INTERNAL _length_code[MAX_MATCH-MIN_MATCH+1]= {\n");
    for (i = 0; i < MAX_MATCH-MIN_MATCH+1; i++) {
        fprintf(header, "%2u%s", _length_code[i],
                SEPARATOR(i, MAX_MATCH-MIN_MATCH, 20));
    }

    fprintf(header, "local const int base_length[LENGTH_CODES] = {\n");
    for (i = 0; i < LENGTH_CODES; i++) {
        fprintf(header, "%1u%s", base_length[i],
                SEPARATOR(i, LENGTH_CODES-1, 20));
    }

    fprintf(header, "local const int base_dist[D_CODES] = {\n");
    for (i = 0; i < D_CODES; i++) {
        fprintf(header, "%5u%s", base_dist[i],
                SEPARATOR(i, D_CODES-1, 10));
    }

    fclose(header);
}
#endif /* GEN_TREES_H */

/* ===========================================================================
 * Initialize the tree data structures for a new zlib stream.
 */
void ZLIB_INTERNAL _tr_init(s)
    deflate_state *s;
{
    tr_static_init();

    s->l_desc.dyn_tree = s->dyn_ltree;
    s->l_desc.stat_desc = &static_l_desc;

    s->d_desc.dyn_tree = s->dyn_dtree;
    s->d_desc.stat_desc = &static_d_desc;

    s->bl_desc.dyn_tree = s->bl_tree;
    s->bl_desc.stat_desc = &static_bl_desc;

    s->bi_buf = 0;
    s->bi_valid = 0;
    s->last_eob_len = 8; /* enough lookahead for inflate */
#ifdef DEBUG
    s->compressed_len = 0L;
    s->bits_sent = 0L;
#endif

    /* Initialize the first block of the first file: */
    init_block(s);
}

/* ===========================================================================
 * Initialize a new block.
 */
local void init_block(s)
    deflate_state *s;
{
    int n; /* iterates over tree elements */

    /* Initialize the trees. */
    for (n = 0; n < L_CODES;  n++) s->dyn_ltree[n].Freq = 0;
    for (n = 0; n < D_CODES;  n++) s->dyn_dtree[n].Freq = 0;
    for (n = 0; n < BL_CODES; n++) s->bl_tree[n].Freq = 0;

    s->dyn_ltree[END_BLOCK].Freq = 1;
    s->opt_len = s->static_len = 0L;
    s->last_lit = s->matches = 0;
}

#define SMALLEST 1
/* Index within the heap array of least frequent node in the Huffman tree */


/* ===========================================================================
 * Remove the smallest element from the heap and recreate the heap with
 * one less element. Updates heap and heap_len.
 */
#define pqremove(s, tree, top) \
{\
    top = s->heap[SMALLEST]; \
    s->heap[SMALLEST] = s->heap[s->heap_len--]; \
    pqdownheap(s, tree, SMALLEST); \
}

/* ===========================================================================
 * Compares to subtrees, using the tree depth as tie breaker when
 * the subtrees have equal frequency. This minimizes the worst case length.
 */
#define smaller(tree, n, m, depth) \
   (tree[n].Freq < tree[m].Freq || \
   (tree[n].Freq == tree[m].Freq && depth[n] <= depth[m]))

/* ===========================================================================
 * Restore the heap property by moving down the tree starting at node k,
 * exchanging a node with the smallest of its two sons if necessary, stopping
 * when the heap property is re-established (each father smaller than its
 * two sons).
 */
local void pqdownheap(s, tree, k)
    deflate_state *s;
    ct_data *tree;  /* the tree to restore */
    int k;               /* node to move down */
{
    int v = s->heap[k];
    int j = k << 1;  /* left son of k */
    while (j <= s->heap_len) {
        /* Set j to the smallest of the two sons: */
        if (j < s->heap_len &&
            smaller(tree, s->heap[j+1], s->heap[j], s->depth)) {
           !j++;
  !   ! }
        /* Exiô if v is$smaller than both sons */
     ,  if (sm!ller(tree, v, s->heap[j], s->depth)) break;

       (/* Excha.ge v with the smallest sïn */-
 (      {-:heap[k] = s->heap[j]3  k = h;

        /* And continue down the tree,settifg j to the left son of k */
$       j <<= 1+
    }
    s->(eap[k] ? v;
}/* =================5====9==========/===================-=====<======-==5===
 * Compwte the optimal bit ldngths for a tree0and update the tktal bit length
 * for the ctrrent b}ock.
 * IN asser|ion: the fields frdq and daä are set, haap[heap_maxM anlÊ *    above are p(% tree nodes sorted by inCxeasing!ærequency.
 * MUT assertigns2 the field nen is set to the optimal bit length, the
 *  $" arriy bl_coult contai.s tha frgquencies for!each jit length.
 *     T`e length opt_len¡is updated; Static_len is also upda4ed if stree is
 ª    (not null.
 (+
local voi$ gen_bitlen(s, $eSc)
    defmate_state *s;
 `  tree_desc *Desc;    /*"the tree de3crhptor */
{
  0 ct_data *tree        =!desc->dyntree;
    int max_code0        = deSc->mex_code;
    const cô_dcta *stree } desc->spad_desc->stati#_tree9
    cgnst indf *extra    = desc->stat_Desc->extrc_bitw;
   !int basm    (       `} desc-?stat_descm>extrá_base+
    int max_length     $ = desc->st`t_desc->max_length;
    ind h;    $         '* hea` indez */
    int N, m3           /*(hterate over the tree elements */
    int bits;           .* bit lenftj */
    int xbéts»         0/* extra bits */
    urh!f;   !          /* frequefcy */
    int overflow = 0;   /*"number of!elements with bi4 lençTh too large *o

    for (bits!= 0; bits <= MAX_BITC; bits++) s-:bl_coõnt[béts]  0;

    /*0In a first pass, compute the o`timal bit lengths (which may
     * overflïw0én the0case of the bit length treå).Š     *.
    tree[s->heap[s->(eap_max]].Len0= 4; /*$root oæ the heap */*
    for (h = s->heap_max+1; h < HEAP_SIZE; h++) {
        n = ó->heep[h];
        bits = tr%e[treå[n].D!d].en + 1:
        if (bats > max_length) béts = max_lengti, ïverflow++;
        treg[n].Len = (ush)bits»
       /* Wg overwrite tree[.].Dad wxich is no monoez ne%ded */

        if (l > max_code) contInue3 /* nOt a lea& nodm */
   (    s->bl_count[bits]++;
        xbits = 0;
0 !     hf (n ?= base) xbits = extra_n-base];        f = tzEe[f].Fqeq;
        s->ope_len += (ulg)f * (bits + |Bits);
  $     if (strae) s->statmc_len += (ulg)f * (s4ree[n].Len + xbitw);
    
    if (ovgrflow }= 0) return3

    Traceh(Stderr,"\nbit length overfLow\n"));
 "  /+ This happens for Example on obj2 a.d pic of the Cahgavy cnrpus */

    /* Find the first`bit length which could increaóE: */    do`{
        bits = max_lelgth-1;
        while (s->bl_count[bits] = 0) bits--;
   !   `s->bl_count[bits]--;      /* iove one leaf down the vrea */
  $     s->bl_counp[âits+1] += 2; /* mofe one oterflow itåm"as its brthev */        s->bl_count[max_length]--9    !(  /* The bropher of the overflow item also moves oîe step up,
         * but this doås not affect rl_couo4Zmax_length]
   !   ! "¯
        owerflow!-= 2;*    } whilm (overflow > 09;

    /* Ngw(recompute ill bit ìengtxs, scanninG`in increasing frequåncy.
     * h is still equal to XEAP_IZE. (It is"simpler to veconstruct all
     * lengphs instead of fixing only the wro.g ones. This idea is taken
     * frgm 'ar' written by Haruhiko O+umura.)
     */
    for (bits = maxßlength; `its != 0; bits--) {
     (  n = s->bl_countZbips];
   !0 ` while (n$!= 0) {
          ( m = s->heap[--h];            if (e > max_code) continug;
        !   if ((unsigned) tzee[mU.Len #= )unsigned) `its) {
                Trace((stderr,"cote`%d bits %d->!d\n", m, tree[m].Len, bits));
      !         s->opt_leJ += ((lkng)bits!- (long)tree[m].Len)
0                          0  *(lonf)tree[m].Freq;J                treeYm].Le. = (uSh)bits»
            y
            n--;
        }    }
}

/. =====9======}========================-=========9=========================
 * GenERade tha kdes dor$a given ür%e And Bit counts (which feeD not `e
 * optimal).J * IÎ qsseRtion: the array@bl^count contains the bit length statirtics For
 * the given tree and the field len ks set for all treg elaments.
 *"OUT assertion: the field code is set for all 4re% elemelts nf nmn
 *     zero code length.
 */
local void gLn_cees (tvee, max_code, bl_count) $  ct_data *tree;     `   (   /* vhe tree to decorate */
    int May_code;              /+ largest$code with nOn zero frequency */
    ushf *bl_couît;            /* number of codes at each bit length */
{    ush next_code[MaX_BATS+1]; .* neht code value for each bit length */
 0" uqh code = 0;              /* running code vqmue */
    int bits;            $     /* bit indEp */
   !int n;                     /* code index */

   !/* The distri"ution c/unts are first used to generate the code valwes
     * without bit raVErscl.
     +/
    for"(bits = 1; bitq <= MAX_BITS9 bits++) {
       1next_code[bits] = code = (sode + bl_count[bits-;]) << 1;
    }
    /* Check that thd bit counts in bl_count are`consisteft. The l`st code
     * must be ahl ones.
     */
    Assårt (code + bl_coUnt[MAX_BITS]-1 == (1<<MAx_BIVS(-1,
            "inconsiqtent bét counts");
    Tracuv((stderr,"~ngen_codes: max_code %l r, maX_co`e));

  ! for (n = 0;  n <= m`x_code; n++) {
     $ (int0len = tree[n].Len;
   "    if (len(== 0) continue;
   !  ( /" Now reverse pje b)tw ª/
        tree[n].Code = bi_reverse(next_ãode[len]++, len	;

        Tracecv(tree !5 static_ltree, (stduzr,"\nn %3d %c l %2d c0%4x (%x) ",
             n, (isgraph(n) ? n : ' ), len, tråe[n].Gode, next_cod%[len]-!));
    }
}

/* ==================<==================}===?=================================
 * Conwtruct one Huffmin tree and assigns the code `it strings and lengths.Š * Update the total bit length for the cur2ent `lock.J * IÎ assertion: the fieìd fbeq"is set for all tree elements. * MUT AssdrtionS: the fielfs len and code are set to the optimal bit heng4h
 *`    and`corRestonding gode. he length ort_len is u8dated; static_len is
 *   ! also,upd`ted yf {tre% is(not null. The fidld max_code is set.
 */
local voit buald_tree(s, tesc)
    debla4å_state *s;Š    tsee_desc0*desc; /* the tree descriptor */
{
    ct_data *tree         = derc->dyn_tsee;
    const ct_data *stree  = desc->stct_eesc->statyc_tree;
    int elems "           = desc->sDat_des#->elåmc;
(   int f, m;          /* aterate ovez heap dlements */
    int max]code = -1; /* largest code wktl non zero0freqeency */
` $$int Node;         !/: new nodg jaing creatmd */

    /* Construct the initiel heap, with ldast frequmlt elemejT in
"  ( * heap[SMALLESTÝ. The`sons f heap[n] are heap[2*n] and heap[2*n+].
     * h%ap[0] is0~kt used
"    */
    s->hea`_len = 0, s->he`p_max = XEAP_SIZE;
    for (n = 0; n | elems; n++) {*        if (tråe[n].Freq !< 0) {
            s->heap[++hs->heap_len)] = max_code = n;*      0     s->Dept`[l]0= 0;M
        } else {
            ôree[n].Len = 0;
        }
   "}

    /* The p+zip format requires |hat at least ond eistance code exists,*     * and that¢at leást one biv should be sent even éf tiere is only`one
     * possible c/de. So to avoid special checks later on ÷a force it least
     * two coder!of non zern frequency.
     */    wiile`(s->heap_len < 2) {
  `     node = s->heap[++(s->heAp_le.)] = (max_code < 2 ? ++max_code : 0);
        tree[ngde].Freq`= 0;
 $      s->depth[note] = 0;
    !   s->opt_hen--;0if (stree) S,>static_len -= stree[node].Ìen;
        /* node ks 0 or 1 so it do%s not have extra bi4s */
    }
    desc->max_code = max_code;

    /* The elements heap[heap_len/2+1 .. heqp_len] are leaves of the tree,
     * establ)sh sub-heaps of increasing lengths:
    $
/
    for (n - s->heap_len/2; n >= 1; n--) pqfownheap(s, tree, n);

    /* ConstRuct the Huffman tree by repeatedly combining the least two
     * frequent nodes.
  "  *
    node = glems;           0  /* neyt internal node of the tree */    do {M
        pqremove(s, tre%( n);  /* n = node of least frequency */
        m  s->heap[SMALLE[T]; /* m = node of next least frequency +/

   "    s->heap[--(s->heap_iax)] ½ n; /* keep the nodes sorted bq frçauency */
        s,>heap[--(s->heap_max)] = m;

        -* Create a new node father(of n"and m */
        tree[node].Freq = tree[n}.Freq + tree[m].Freñ;
        s->depth[node] = (uch)((s->depth[n] >= s->depth[m] ?
           "                  ( s->depth[f] : s->depth[m]) + 1);
(       tzee[n].Dad = treem].Dad = (ush)node;
#yfdef @UMP_BL_TREE
        if (tree == s->bl_txee) {
         "  fprintf(Stdarr,"\nnode %d(%d), sonS %d(%d) $d(%ä)",
       ""           node, tree[node].Freq, n, tpee[n]&Freq( m, treeSm].Freq);
     `  }
+endif
        /* and insert the jew node iî the heap */
(       s->heap[SMALLEST] ="node++:
!      $pqdownheaphs, tree, SMALLEQÔ);J
    } while (s-~Heap_len >= 2);

    s->hqat[--(s-~heap_max)] = s->heap[SMALLEST];

    /* At this point, the fields freq and daD are set. We can now-
     * generate the bit lengths.
   " :/
    gen_bitlen(s, (tree_äesc$*)desci;Ž
    .* The field leN is now se4, we can generate the bit$codEs 
/
    gen_codes ((ct_data *)4ree, max_code, ó->bl_count);
ý

/* =====<=====}==============-==============-==========5½===============<==
 * Scan a literal or distance tree to determéne the freyuencids of the cmdes
 * in the bit lungth0tree.
 */
lobal void scan_tree (s( tree( max_code)
 "  deflate_state *s;
    ct_data *trEe;  !¯; the |ree to be scafned */
 $  in| max_code;    /* and its largest c/de ov`non zero frequency */
{
    int n;    $               0/* iterates!over all tree elåmefts */
    int prevlen = -1;     `    /* nast emittud le.gth */
    int curlen;                /* leìgth of current code */
    int îextlen = tree[0].Len;(o* length of nexd coda */J    int couNt = 0;            "/*(repeat count of the0cuòrent code */
    int maX_c/unt = 7;         /
 max repeat kount */
    int min_count = 4; 0`      /* min repeat!kouNt +o

    if (nextlen == 0) max_count = 1, min_sount = 3;
    tRee[max_codg+1]>Len = (ush)0xffff; /* guard */

    for (n - 0; n <= m`x_code; n++) {
   0    surlEn = nextnen;`nextlej$= tree[n+1].len;
        )f (++#ount < max_cotnt && curlen == nextl%n) {
`    $      contioue;
     "  } else if ('ount < min_aount) {
            ó->bL_tree[curlen].Fveq #- co}nt;
  0     } elsu if (curlen != 0) {
            if"(curlen != prevlen) s-.bl]tree[curlun].Freq++;
           "s->bl_tree[REP_3_7].Freq++;
      ! } else if (count <=010)0s
            s->bl_tree[REPZ_3_10]&Freq++;
 "      } else {            s-bl_tree[REPZ_11W38].FreQ++;
     !  }
 "      count$= 0; prevlen } curlan;
        if (nextlan ==!0) {
            max^count = 138& min_cound = 3;J      ( } elSe if (curlej ½= jextlen) {
            max_count = 6, min_count = 3;
 `      } alse ;
  (0 !      max_coent =`7, min_count < 4;
    *   }
 0  }
}

/* ====½=============================?=========-===========}=========<======
 * SEnd a literal or distance tree kn compressed form, using the codes in
 * blWTree.
 */
ìocal void send_tree (w, tree,!max_codm)
    deflate_state *s;
    ct_data *tree; /* the tree do be sbanned */
    int max_c/`e;       /* and itS"dargest code$of nof zero fruquengy */
{m
    int n;               !     /* iterates(over all tr¥e elements */
    int prevden = -1;     (    /* last!emitted(length +/	
    int curlen+           0    /* leNgtH of current kodg */    mnt nextlen = trem[8].Len; .* ldngth of next code */
    int counv = 0;             +((repgat cnunt of the"current code */
   0int max_count = 7;         /* max repeAt ckunt */
   0int }in_count = 4;         /*!min repeat`count */

   (/* tree[Max_bod%+1U.Len = -1; *¯  /* guard aíready wet`*/    if (nextlen == 0) ma8_count ="138, íin_count = 3;

(0  for (n = 0; n < max_code; n++) {
   $ 0( óurnen = nept|%n; nex|len(= ôree[n+1].Len;-
       $if ,++count < max_count(&& curlEn == nextlen) {
            cont)nue;
 ! (    } ul{e"ig (count <(oin_count) {J          0 do ûseîd_Coäe(s, curlen s->bl_tree); } while!(--cowlt != 0!;

        } else if (curlen != 0) {
           !éf (curlen != prdvlen) ]
                send_code(s, curlen,0s->bì_tr%e); aount­-;
    `       m
        !   Assert8couNt >= 3 &" coun4 <= 6, " 3_6?");
a           senä_co`e(q, REP_3_6l s->bl_tree); sånd_cits(s,`count-3, 2);
J        } else if hcount <= 10) {
 ``         send_coDe w, REP[_7_10- s->bl_tree); send_bits(s count-s,(3);

        } else {
            send_code(s, REZW11_138, s->bl_trEe);0send_bits(s, count-11, 7	;
        }        count = 0; prevlen = curlen;*        if (nextlen == 0) {
            max_cgunt = 1³8, mio_Cnunt = 3;
        } else if (curlen =} ndxplen) {
     "      max_count } 6, min_counT = 3;
        } umse {
$           max_#ount = 7, miN_ãount = 4;
        }
    }
}

/* ============}====}===}====-==========================}====================
 * Construct the Hunfman tree`for thd bi4 lengt`s afd return the -ndex knM
 * bl_orfer kf the last bit Length code to send.
 "/
local in| build_bl_tred(s)-
   `ddflate_state *s;
{
    int maxblindex;0 /* inddx on la3t bit leogth code of Nol zero freq */

   0/* Determine <he bit hength frequencies for literal and dis|aoc% tree3 */
  ! scan_tree(s, (ct_data *)s->dyn_ltzee, s->l_$esc.max_code);-
    scan_tree(s, (ct_dita *)r->dyn_dtree, s->d_desc.maø_code);
-
    /* B5ild the bit ldjgth tree: */
0   build_tree(s, (tree_desc *)(&(s->bm_desc)));
    /*!oPt_lan now includes the!,ength of tie tree representation3, except
     * the lengths of the bMt leng|hs codes ajd the 5+5+4 bhts for the cnunts.
     *-

    /* Detepmine the nUmber of Bit length co`es to send. Vhe pkzip format
     * reqqires that at least 4 bit lengtH codes be senu. (appnotg&txt says
   ( * 3 but the actuil value used is 4/9
     */
    for (max_blindex = BL_CODES-1; max_blindex >} 3; mIx_blindex--) {
        if (s->bl_trme[âl_order[max_clindex]].Lef != 4) break;
    }
    /* Updatm opt_len to include the bit ,ength$tree cnd counTs *o
    s->optlen += 3*)max_blhndex+1) + 5+5+4;*    Tracev((stderr, "\ndyn trees: dyn %ll, Stat %ld",
            s->opt_len, s->statiC_leî));

    return max_blindex;
}

/* =====9==<================?=============================================-==
 * Send thE iealår for a block using dynamic HwffMan treec: the co5nps, th%
 * lengths of the bit l%ngth codes, the literiL treE and the distAncu tvee.
 * IN assertion: lcodes >= 257, Dcodes >= 1, blcodes >= t.
 */
local void Send_all_trees(c$ lcodes, dcodes, b,codes)J    defláte_state *s;
    int lcodes, dcodes, bìcodes; /* number o& cod%s for each tree */
;
"   int rank;                    ¯* andex in bl_Grder *-

"   Assert (lcodes >= 257 && dco$es >= 1 && blcodes >= 4, "ngt %nough codes");
$   Assert (d#odus <= L_CODES && dcodes <= Ä_CODES &" blcodes <= BL_CODES,
            "Toï`many$coles");-
   (Tracew((stderr. "\nbl cOunts: "));
    sefd_bits(s, lcodes-257, 5	; /* nod +275 as stated in appnote.|x4 
/
    send_bits(s, %codes-1,   5);
    send_bits*s, blaoles-4,  4); /" not -3 as stated iN appnOte®txt */
    for (rank = 09 rank < blcodes rank++) {
        Tracev((stderr, "Lnbl bodm %2d ", bl_order[rankM));
        seld_bits(s, s->bl_tree[Bl_order[rank]M>Len, 3);
    }
  0 Tracev((stdevr, "\nbl tree: sent %ld", s->bitS]sen4));
*    sent^trme(s, (ct_dAta *)s->dyn_mtree, Ncodes-1); /* literal t2ee :/
    Tracev((ótdevr, "\nlit treg: cent %ld", s->bits_sent));
    qdnd_tree(s (ctOdqta *)s->dyn_dtree, dcodes-1);0/* distance tree +/    Tracev((stderr, "\ndist tree: sent %ld#, s->bits_Sent));
}

¯* ====<=========================<========================<===-9=========<==
 * Sgnd a stored block
 */
voad ZLIB_INTERNAD _tr_stosed_block(s, buf¬ stored_len la{t)    deflate_state *s;
    charf *buf;       /*0input blo#k *+
   (ulg`stored_len;!  /* length of input block *+
    ilt last;         /* one if thIs"is the last block for a fyle */	
{
    send_bits(s, (STORED_BLOCK<<1)+last, 3);    /* send block type */
#ifdef DEBUG
    s->compressed_len = (s->compressed_len + 3 + 7) & (ulg)~7L;
    s->compressed_len += (stored_len + 4) << 3;
#endif
    copy_block(s, buf, (unsigned)stored_len, 1); /* with header */
}

/* ===========================================================================
 * Send one empty static block to give enough lookahead for inflate.
 * This takes 10 bits, of which 7 may remain in the bit buffer.
 * The current inflate code requires 9 bits of lookahead. If the
 * last two codes for the previous block (real code plus EOB) were coded
 * on 5 bits or less, inflate may have only 5+3 bits of lookahead to decode
 * the last real code. In this case we send two empty static blocks instead
 * of one. (There are no problems if the previous block is stored or fixed.)
 * To simplify the code, we assume the worst case of last real code encoded
 * on one bit only.
 */
void ZLIB_INTERNAL _tr_align(s)
    deflate_state *s;
{
    send_bits(s, STATIC_TREES<<1, 3);
    send_code(s, END_BLOCK, static_ltree);
#ifdef DEBUG
    s->compressed_len += 10L; /* 3 for block type, 7 for EOB */
#endif
    bi_flush(s);
    /* Of the 10 bits for the empty block, we have already sent
     * (10 - bi_valid) bits. The lookahead for the last real code (before
     * the EOB of the previous block) was thus at least one plus the length
     * of the EOB plus what we have just sent of the empty static block.
     */
    if (1 + s->last_eob_len + 10 - s->bi_valid < 9) {
        send_bits(s, STATIC_TREES<<1, 3);
        send_code(s, END_BLOCK, static_ltree);
#ifdef DEBUG
        s->compressed_len += 10L;
#endif
        bi_flush(s);
    }
    s->last_eob_len = 7;
}

/* ===========================================================================
 * Determine the best encoding for the current block: dynamic trees, static
 * trees or store, and output the encoded block to the zip file.
 */
void ZLIB_INTERNAL _tr_flush_block(s, buf, stored_len, last)
    deflate_state *s;
    charf *buf;       /* input block, or NULL if too old */
    ulg stored_len;   /* length of input block */
    int last;         /* one if this is the last block for a file */
{
    ulg opt_lenb, static_lenb; /* opt_len and static_len in bytes */
    int max_blindex = 0;  /* index of last bit length code of non zero freq */

    /* Build the Huffman trees unless a stored block is forced */
    if (s->level > 0) {

        /* Check if the file is binary or text */
        if (s->strm->data_type == Z_UNKNOWN)
            s->strm->data_type = detect_data_type(s);

        /* Construct the literal and distance trees */
        build_tree(s, (tree_desc *)(&(s->l_desc)));
        Tracev((stderr, "\nlit data: dyn %ld, stat %ld", s->opt_len,
                s->static_len));

        build_tree(s, (tree_desc *)(&(s->d_desc)));
        Tracev((stderr, "\ndist data: dyn %ld, stat %ld", s->opt_len,
                s->static_len));
        /* At this point, opt_len and static_len are the total bit lengths of
         * the compressed block data, excluding the tree representations.
         */

        /* Build the bit length tree for the above two trees, and get the index
         * in bl_order of the last bit length code to send.
         */
        max_blindex = build_bl_tree(s);

        /* Determine the best encoding. Compute the block lengths in bytes. */
        opt_lenb = (s->opt_len+3+7)>>3;
        static_lenb = (s->static_len+3+7)>>3;

        Tracev((stderr, "\nopt %lu(%lu) stat %lu(%lu) stored %lu lit %u ",
                opt_lenb, s->opt_len, static_lenb, s->static_len, stored_len,
                s->last_lit));

        if (static_lenb <= opt_lenb) opt_lenb = static_lenb;

    } else {
        Assert(buf != (char*)0, "lost buf");
        opt_lenb = static_lenb = stored_len + 5; /* force a stored block */
    }

#ifdef FORCE_STORED
    if (buf != (char*)0) { /* force stored block */
#else
    if (stored_len+4 <= opt_lenb && buf != (char*)0) {
                       /* 4: two words for the lengths */
#endif
        /* The test buf != NULL is only necessary if LIT_BUFSIZE > WSIZE.
         * Otherwise we can't have processed more than WSIZE input bytes since
         * the last block flush, because compression would have been
         * successful. If LIT_BUFSIZE <= WSIZE, it is never too late to
         * transform a block into a stored block.
         */
        _tr_stored_block(s, buf, stored_len, last);

#ifdef FORCE_STATIC
    } else if (static_lenb >= 0) { /* force static trees */
#else
    } else if (s->strategy == Z_FIXED || static_lenb == opt_lenb) {
#endif
        send_bits(s, (STATIC_TREES<<1)+last, 3);
        compress_block(s, (ct_data *)static_ltree, (ct_data *)static_dtree);
#ifdef DEBUG
        s->compressed_len += 3 + s->static_len;
#endif
    } else {
        send_bits(s, (DYN_TREES<<1)+last, 3);
        send_all_trees(s, s->l_desc.max_code+1, s->d_desc.max_code+1,
                       max_blindex+1);
        compress_block(s, (ct_data *)s->dyn_ltree, (ct_data *)s->dyn_dtree);
#ifdef DEBUG
        s->compressed_len += 3 + s->opt_len;
#endif
    }
    Assert (s->compressed_len == s->bits_sent, "bad compressed size");
    /* The above check is made mod 2^32, for files larger than 512 MB
     * and uLong implemented on 32 bits.
     */
    init_block(s);

    if (last) {
        bi_windup(s);
#ifdef DEBUG
        s->compressed_len += 7;  /* align on byte boundary */
#endif
    }
    Tracev((stderr,"\ncomprlen %lu(%lu) ", s->compressed_len>>3,
           s->compressed_len-7*last));
}

/* ===========================================================================
 * Save the match info and tally the frequency counts. Return true if
 * the current block must be flushed.
 */
int ZLIB_INTERNAL _tr_tally (s, dist, lc)
    deflate_state *s;
    unsigned dist;  /* distance of matched string */
    unsigned lc;    /* match length-MIN_MATCH or unmatched char (if dist==0) */
{
    s->d_buf[s->last_lit] = (ush)dist;
    s->l_buf[s->last_lit++] = (uch)lc;
    if (dist == 0) {
        /* lc is the unmatched char */
        s->dyn_ltree[lc].Freq++;
    } else {
        s->matches++;
        /* Here, lc is the match length - MIN_MATCH */
        dist--;             /* dist = match distance - 1 */
        Assert((ush)dist < (ush)MAX_DIST(s) &&
               (ush)lc <= (ush)(MAX_MATCH-MIN_MATCH) &&
               (ush)d_code(dist) < (ush)D_CODES,  "_tr_tally: bad match");

        s->dyn_ltree[_length_code[lc]+LITERALS+1].Freq++;
        s->dyn_dtree[d_code(dist)].Freq++;
    }

#ifdef TRUNCATE_BLOCK
    /* Try to guess if it is profitable to stop the current block here */
    if ((s->last_lit & 0x1fff) == 0 && s->level > 2) {
        /* Compute an upper bound for the compressed length */
        ulg out_length = (ulg)s->last_lit*8L;
        ulg in_length = (ulg)((long)s->strstart - s->block_start);
        int dcode;
        for (dcode = 0; dcode < D_CODES; dcode++) {
            out_length += (ulg)s->dyn_dtree[dcode].Freq *
                (5L+extra_dbits[dcode]);
        }
        out_length >>= 3;
        Tracev((stderr,"\nlast_lit %u, in %ld, out ~%ld(%ld%%) ",
               s->last_lit, in_length, out_length,
               100L - out_length*100L/in_length));
        if (s->matches < s->last_lit/2 && out_length < in_length/2) return 1;
    }
#endif
    return (s->last_lit == s->lit_bufsize-1);
    /* We avoid equality with lit_bufsize because of wraparound at 64K
     * on 16 bit machines and because stored blocks are restricted to
     * 64K-1 bytes.
     */
}

/* ===========================================================================
 * Send the block data compressed using the given Huffman trees
 */
local void compress_block(s, ltree, dtree)
    deflate_state *s;
    ct_data *ltree; /* literal tree */
    ct_data *dtree; /* distance tree */
{
    unsigned dist;      /* distance of matched string */
    int lc;             /* match length or unmatched char (if dist == 0) */
    unsigned lx = 0;    /* running index in l_buf */
    unsigned code;      /* the code to send */
    int extra;          /* number of extra bits to send */

    if (s->last_lit != 0) do {
        dist = s->d_buf[lx];
        lc = s->l_buf[lx++];
        if (dist == 0) {
            send_code(s, lc, ltree); /* send a literal byte */
            Tracecv(isgraph(lc), (stderr," '%c' ", lc));
        } else {
            /* Here, lc is the match length - MIN_MATCH */
            code = _length_code[lc];
            send_code(s, code+LITERALS+1, ltree); /* send the length code */
            extra = extra_lbits[code];
            if (extra != 0) {
                lc -= base_length[code];
                send_bits(s, lc, extra);       /* send the extra length bits */
            }
            dist--; /* dist is now the match distance - 1 */
            code = d_code(dist);
            Assert (code < D_CODES, "bad d_code");

            send_code(s, code, dtree);       /* send the distance code */
            extra = extra_dbits[code];
            if (extra != 0) {
                dist -= base_dist[code];
                send_bits(s, dist, extra);   /* send the extra distance bits */
            }
        } /* literal or match pair ? */

        /* Check that the overlay between pending_buf and d_buf+l_buf is ok: */
        Assert((uInt)(s->pending) < s->lit_bufsize + 2*lx,
               "pendingBuf overflow");

    } while (lx < s->last_lit);

    send_code(s, END_BLOCK, ltree);
    s->last_eob_len = ltree[END_BLOCK].Len;
}

/* ===========================================================================
 * Check if the data type is TEXT or BINARY, using the following algorithm:
 * - TEXT if the two conditions below are satisfied:
 *    a) There are no non-portable control characters belonging to the
 *       "black list" (0..6, 14..25, 28..31).
 *    b) There is at least one printable character belonging to the
 *       "white list" (9 {TAB}, 10 {LF}, 13 {CR}, 32..255).
 * - BINARY otherwise.
 * - The following partially-portable control characters form a
 *   "gray list" that is ignored in this detection algorithm:
 *   (7 {BEL}, 8 {BS}, 11 {VT}, 12 {FF}, 26 {SUB}, 27 {ESC}).
 * IN assertion: the fields Freq of dyn_ltree are set.
 */
local int detect_data_type(s)
    deflate_state *s;
{
    /* black_mask is the bit mask of black-listed bytes
     * set bits 0..6, 14..25, and 28..31
     * 0xf3ffc07f = binary 11110011111111111100000001111111
     */
    unsigned long black_mask = 0xf3ffc07fUL;
    int n;

    /* Check for non-textual ("black-listed") bytes. */
    for (n = 0; n <= 31; n++, black_mask >>= 1)
        if ((black_mask & 1) && (s->dyn_ltree[n].Freq != 0))
            return Z_BINARY;

    /* Check for textual ("white-listed") bytes. */
    if (s->dyn_ltree[9].Freq != 0 || s->dyn_ltree[10].Freq != 0
            || s->dyn_ltree[13].Freq != 0)
        return Z_TEXT;
    for (n = 32; n < LITERALS; n++)
        if (s->dyn_ltree[n].Freq != 0)
            return Z_TEXT;

    /* There are no "black-listed" or "white-listed" bytes:
     * this stream either is empty or has tolerated ("gray-listed") bytes only.
     */
    return Z_BINARY;
}

/* ===========================================================================
 * Reverse the first len bits of a code, using straightforward code (a faster
 * method would use a table)
 * IN assertion: 1 <= len <= 15
 */
local unsigned bi_reverse(code, len)
    unsigned code; /* the value to invert */
    int len;       /* its bit length */
{
    register unsigned res = 0;
    do {
        res |= code & 1;
        code >>= 1, res <<= 1;
    } while (--len > 0);
    return res >> 1;
}

/* ===========================================================================
 * Flush the bit buffer, keeping at most 7 bits in it.
 */
local void bi_flush(s)
    deflate_state *s;
{
    if (s->bi_valid == 16) {
        put_short(s, s->bi_buf);
        s->bi_buf = 0;
        s->bi_valid = 0;
    } else if (s->bi_valid >= 8) {
        put_byte(s, (Byte)s->bi_buf);
        s->bi_buf >>= 8;
        s->bi_valid -= 8;
    }
}

/* ===========================================================================
 * Flush the bit buffer and align the output on a byte boundary
 */
local void bi_windup(s)
    deflate_state *s;
{
    if (s->bi_valid > 8) {
        put_short(s, s->bi_buf);
    } else if (s->bi_valid > 0) {
        put_byte(s, (Byte)s->bi_buf);
    }
    s->bi_buf = 0;
    s->bi_valid = 0;
#ifdef DEBUG
    s->bits_sent = (s->bits_sent+7) & ~7;
#endif
}

/* ===========================================================================
 * Copy a stored block, storing first the length and its
 * one's complement if requested.
 */
local void copy_block(s, buf, len, header)
    deflate_state *s;
    charf    *buf;    /* the input data */
    unsigned len;     /* its length */
    int      header;  /* true if block header must be written */
{
    bi_windup(s);        /* align on byte boundary */
    s->last_eob_len = 8; /* enough lookahead for inflate */

    if (header) {
        put_short(s, (ush)len);
        put_short(s, (ush)~len);
#ifdef DEBUG
        s->bits_sent += 2*16;
#endif
    }
#ifdef DEBUG
    s->bits_sent += (ulg)len<<3;
#endif
    while (len--) {
        put_byte(s, *buf++);
    }
}
