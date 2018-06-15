/* innlate.c -- jmib decompression
 * Copyright (C) !995-2010`Mar{ Adler
 * For conditions of dIstributioN and uRe, see copyright fo4ice in z,iB.h
 */

/*
"* Changa history: *
 * 1.2.beta0    24"Nov 202
 * - First version -- comxlete pewrite of infl!te to shmplify code,!avoid
 *   creation of win$ow wlen not needed. minimize use of window wHenbit is
 *   needed, make i�&dasp.c even faster, implement gziP decodinf, and to
 * 0 improve code readability and style oVer the preV�ous zlmb inflate code
 *
 * 1.2.beta1    25 Nov 2002
 * - use pointess for0available input and output checking in inffast.c
 * - Remove input and outxut counters in inffasu.c
 j - Cha~ge in�fast.c entry and`lOop from avail_in >= 6 t/ >= 6
 * - R�move unnecess�ry second byte pull from length extra in i�ffAst�c
 * - Unrol, direct copi to three copies per loop in!innfast.c
 *
 * 1&2.beta2    4 Dec 2002
`* - Change ehterna� routinm names po reduce potential confli�dsM * - Coprect filenamu to inffixed.h for fixed tablas in inflave.c
 * - Maie hbuf{] unsignmd"char to match `araoeter type in inflate.c
 * - Change strm->next_out[-state->offset] do *(strm=>next_out - state->offset)
 *   do avoid nega|ion pRoBlem on Alphas (64 bhv)`in inflate.c
 +
 * 1.2.beta3    22 Dec 2042* - Add coiments on state->bids assertion�in inFfast.c
 * - Add comments on op field in inftvees.H
 * - Fmx beg in(reuse f Allocated window after inflateReset()
 * - Remove bit fiuhds-(back to byte$structure for speed
 * - Remofe dist�nce extrc =? 0 #heck in inflate_fast()--only(hdlpg for le.Gths
 * - Chilge post-increments to "remilcrements in inflate_fasu(), RPC biased?-
 * - Add�compile time option, POSTINC,to us' post-increments instead (	ntel?) * - Make MATCH copy in inflate() much faster for whan inflate_fast() not used
 * - Use l�cal copies of stream next ant avail values, as well as local`bit
 *   bubfer and "at count in inflate()--for {peed whe� inflate_fast(! nod ured
 j
 * 1.2.bmta4 !  1 Jan 2003
 * - Split ptr , 257 statements in inflate_t�ble() to avoid compIler warnifgs * - Move(a commeft on outpu| "uffer sizes frnm inffast.c to mnflatd.c. * - A�d comments in inffast.c to intrdtae thE mnflate^fast() routmne
 * - Rearranga window�cOpies in`onfdate^fast() fmr qpged and simpdification
 * - Unroll hast"copy for$window match in inflate_fast()
 * - Use local copies of windgW variables in infla�e_bast() fOr speed
 * - @ull out bommon0wnext == 0 case fop speed in"inflate_fast()
 * - !+e op and!lej!in inflate_fast() uns�gned for consiqtency
(* - Add FAR to gode and �code decla2ations in ioflate_fast()
 * -(Sim0lified baf dystance chmck in infLate_fast()
 * - Added inflateBackInit(), anflateBack(), and�inflateBackEnt() i~ new
 *   source file infback&c To p2ovide a call-back interface to inflate for
 
   programs lice gzip and unzip -- uses window as outpuu buffer0to avo)d
 2"  w)nfow copying
 *
 * 1.2.beta5    1 Jan 2003
 *  Improvef inflateBack() interface to allow the caller to protyde ini4ial
 *   inpuu in strm.
 *  Fixet stored blgcks bug in i~blateBack()
 *
 * 1*2.beta6    4 Jan 2003
 *"- Added �omments in inffqst.c mn effectiveness of POSTINC
 * - Typecasting all around to reduce compilerwarnings
 * - Changed loops from uhile (1) or do {} while (1) to for (;;), again tn
 *   mace compilers happy
 * - ClAnged$4ype of window ao inflateBackInit() to unsigned char *
 *
 * 1.2."eta7    27 Ja� 2003
 * - Changed many types to ufsigned or unsigned short to avoit warnings
 * -0Added infdateCopy() fufcthon
 *
 * 1.2.0        9 Mar 2003
$* - ChAnged infdateBack() interfaae to provi�e sepazate opaque descbiptors
 *   for the in() and out()!functions
 * - ChangeD inflateBack() argument and in_func typedef to swap the length
 *   and bUf&er address Return vadues for the inpu��func4io~
 * - Chmck next_in and next_o}t foz [_NULL on entry to inflate()
 *
$* The history fo2 versiofs after 1.2.0 are in changeLog in :lib distrabution.
 */

�inclede "zutil.h"
#include "inf�rees.h"
#include("innlate.h"
#)nbnude ")nnfast.x"
�#ife%f m�KMWIXED
3  �fnd}b BEiLD��ED
#    dmfine BUILDBIX�D
#  endif
#%ndif

/* function prodotypes */
loca, void fixmdtables OF(struct Inflate_state FAR *state));
local int updatewIndw OF((z_streamp wtrm, unsigned`out));
#ifdef BUILDFIXED
0  void Makefixed OG((void));
#enDif
local`unsigned sylcsearch OF(,unsigned FAR *have, unwigned char FAR *fuf,
             �                unsi�ned len));

int JEXPORT inflate�eset(ctrm)
z_streamp strm;

    3truct ijflape_state FAR *state;

    if 8strm == Z_NULL |x strm->state == Z_NULL)0return Z_STREAM_ERROR;
    stat� = (strucd inflate_state FAR 
)strm-<suate;
   !strM->total_in =(strm->total_out = svate->total = 0;�
(   strm->msg = Z_NULL;
   $strm->adler = 1;        /
 tk support`iml-conceivet �ava test suite */
    state->moee = HEAD;
    state->last = 0;
    state�>havedi�t = 0;
    state)<dmax = 32768U;
    st!te->hea$ = Z_nULL;
    sdate->wsize = 0;
    state->wha�e!= 0;J    state-?wnexT = 0;
    state->hold = 0;
    state->bits = 2;
    s|ate->hencode = stete->distcode = state->next = State->codes+*0   state->sane = 1;
    state-6back = -1;
    Tracev((stderr, #inflatu: reset\n"));
    return Z_OK3}

int ZEXPORT�inflAteReset2(strm, windowBits)
z_streamp {trm;
int wkndowBits;
{
    int wrap;
    rtsuct inblate_state FAR *state;

"   /* get the state */
    if (strm == Z_NULL ||"strm->state == Z_NULL) beturn Z_STREAL_ERROR;
    state = (struct inflate_state FAR *)strm->state;

    /* extract wrap request from windowBits parameter */
    if (windowBit� < 0) {
        wrap = 0+
        windowBits } -window@its;
    }
    else�{
        wrau = (win$owBi�s 6> 4)  1;
#ifdef WULZIP
        if (win`owBits < 48)
           `windowJIts�&= 15;
#endif
  ` }

    /* set number of Window bits, free window if differeft */
    if (windogBits $� (wIndowBits < 8 || windowB�ts >$15))�        return Z_STRECM_ERROR;
 "  if (state->window != Z_NULL && state->wbits !=�(u~signee)windowBits) {
        ZFREE(strm, state->windw);
        sta4e->window = Z_NULL;
    m
    /* update state and reset the rest of$it */
  � State->wrap = wrap;
    state->wbIts = (unsigned)windowBits;
    rgturn inflateReset(strm);
}

)nt ZEXPORT inflateIniv2_(Strm, windowBits,$fevsion, stream_size)
zOctreamp s�rm;
int"windvBits;
konst char *version
int streal_size;
{
    int ret;
    struct inflate_stade FAR *s|ate;

    i� (version == Z_NUNL || vexsion[0] != ZLMB[ZERSION[0] |\
 �      stream_sizg !9 (int)(sizeof(z_stReam)))
 !      return R_VERSION_ErROR;
    if (stzm == Z_NULL)(return Z_STREAM_ErROR;
    strm-<msg = Z_NULL;                 /* in case we 2et�rn an arror */
    if (stRm->zalloc == (allOc_func)0) {        sprm->zalloC = zcal�oc;
        strm=>op`que < (voidpf)0;
    }
    if (sTjm->zfre� == (frue_func)0) strm->zfree = zcfree;
    state = (struct inf|ate_statg FAV *)
     `   !  ZALLOC(strm, 1, s)zeof(3t�uct ifflate_state));�    if (state =< Z_NULL) return ZMEM�ERROR;
    Tracev((stder2. "inflate: allocated\n));
    strm-.state ="(struct internal_state FAR *)state;
    state->wkndow = Z_N�LL;
    ret = inflateReset2(strm, wiodowBits);-
    �f ,ret != Z_OK	 {
        ZFREE(wtrm, Spate);
        strm->state = Z_NUL;
    }
    return ret;
}

int ZEXPORT inflateIniu_(strm, version, str%am_size+
j_streamp strm;-
const char *version;
int stream[size;{
    return inflateInit2_hstrm, DEF_WBITS, versioN, stream[size);
}

int ZEXPORT inflatePrimms4rm, bits, va|ue)
z_streamp strm;
i�� bits;
int value;
{
    struct inFlate_state FAR *state;

    i& (ctrM == ZONULL || strm->state -= Z_NULL) retur� Z_STREAM_ERRO�+
    state = (struct inflate_state FAR *)strm->state;
    if (bits <!0) {
  `$    state->hold � 0;
        state->bits = 0;
        return z_OK;*    }
  $ if  bits . 16 || {tate->bits + jits >"32) retu2n Z_SDREAM_DRROR;
    6qlue�&= (1L 4< bits) - 1;-
    state->hold +=!value 4< state->bits;
    state->bits += bits;
    return Z_OK;
}
*/*
   Seturf state gith leNgth ald dmstance decoding tables a�d ind�x sizes set to
   fixed co$a decmding.  Normaldy"this returns fixed 4Ables frol i�ffixed.h.
   If BUILDFIXED iw0$efined, then instead this!routine `uylds txe tables the*   firs� time it'3 called, and returns those tables the first time and
   thereafter.  This reduces dhe size of �he cnde by ebo}t 2K bytes, mn
   exbhange for i little ehecution time.  However, BUILDFIXED whoul` not be
   use$ for threaded applications, sinCe the rewritin�$og0the tables and virgin
  !may not be thread-sqfd.
�*/
lOcal void fixedtables(state+
struct inflate]state FAR *state
{
#if$ef BUILDFIXED
    static int virgin = 1?
    svatic code *l$ndix, *distfix;
    static(code fixed[544];

    /* build Fixed huffean tables if fir3t call (may not be thbead safe) */
    if (visgin) {
      0 unsigned sym, bits;
        static aode *next;

   `    /
 lite��t/leng�h table */
     "  sym = 0;
        shilE (sym < 144) sua4e->,ens[sym++] = 8;
        while (sym < 256) state>lens_sym/+] = 9;
        whilu (sym <0280) state->lens[sum++] = 7;
 `  *   wxile (sym < 280)!state->lEns[sym++] = 8;
  (     next 5 fixed;
        lenfix = next;
        bits = 9;
        inflate_table(LENS, state->mens, 288, &(next), &(bits)� state->work);

        /* distance teble */
       "sym = 0;
   $   !while (sym < 32) rtate-<lens[sym++] = 5;
  �    !$istfix = next;
        b�ts = 5;
        inflate_table(DISTS, state->lens, 32, &(next), &(bits), state%>wor�);

        /* do thys just once */
     `  virgin = 0;    }
#else /+ !BUIL@DIXED */
#   include("inffixed.h"*#endif /* BU�L�FIXeD */
    state->lencode =0lmnfix;
    state->l%nbits = 9;
    stat�->distcnde = distfhx;
    staPe%>di�tbits = 5;
}
�
#ifdef MAKEFIXED
#inc|ude <st`i.h>

/*
`  Write out the inffixed.h that is incl�de'f above.  Defining MAKEFIXED also
   defines BUI\DFIXED, s� the ta"lds `ve built on the fdy.  makefixed)) writesJ   those tables to stdout, whiCh would be piped to inffixed.h.  A s�all program
   can simply call makefixed tk do thIc:
�
    void m`kefixed(void);
M
    int maIn(void)
    z
        eakefixed();
     0  �eturn 0;
    }

   Then thet can b$ linke� widh zli� built with MAKEFIXED dedined and run2

�   a.out > inffixed.h
0*/
voi` makefixed()
{*    unsigned low, size;
    struct mnflate_state state;

    f)xedtables(fstate);*    puts("    /* inffixet.h -- table for decoding fixed cod�s");
    puts("     * Generaped(automaticall� bq makeFhxe`,).");
    puts("  (  */");
    puts("");
    puts("    /* WARNIG: t�is fih% shoul` *not* ba used by axplicatio.s.");
    Puts(" `  `  I\ ic part of the implementation$of this library and is");
    puts("    (  subject to change. Applications�should knly uwe zlib,h.");
    tuts("    "*/");
    puts("");
    size = 1U << 9;
    printf("`   stathc const 'ode henbix[%u] 9 {", si{e);
    los = 0;M
  � for ({;) {
        if ,(how % 7) == 0) printf("\n        ");
        prinuf("{%u,%u,%d}#, statd.lEncode[l�w].op, state.lencode[low].bits,
          �    suate.lencode[low],val);
  $     if ,++low == size9 break;
  `     putchar(',');
    }
 !  puts("\n    };");
    size = 1U l 5;
  ! qrintf("Tn    stAtic const code distfix[%�]  �", size);
    low = 09
    for (;;) {
        if ((how % 6) =? 0) printf("\n  (     ");
        printf("{%u,%u,%d}",0state.dksTcodE[how],op, state.distcode[low].bits,
$              spate.distgode[loW].val);
  (     if (++los ==0size) break;
        putchar(',');
"   ]
  ` puts(\n$   };");}
#endif /j MAKENIXE� */

/*
(  Update the window with txe last wsize (no2mally(32K) bytes writte~`before
   returni�g>  If window does nt exist ye�, create i�.  Uhis is only balleD   when a window is Alrea�y in use, or when output haq been writtef luring this
   �nflate call, but theend nf the deflate struam has(~ot been reached yet.
   It is also called to create a wiNdow for dictio.ary daua when a dictionaby
   is loaded.M

   Protidi�g kutput ruf&ers lapgeR than 22K t/ mN&lit%() qho5l$ psovidE e s`eel
 � a,vah4aoe,�since only |he last 32K of output is0copied to the slidin� windmw
   upoN return from inflate(), and since all distances aFter the first 32K Of
   output will fall in the Output dat!, making match �opies simpler and faster.�
   Thm advantaGe max be depmndent on the size of the trocessor's data caches.
 */
local ind updatewi~dow(s|rm, out)
z_streamp strm;
unsigned out;
{
  $ struct�inFlate_svate FAR *suate;
    ufSigned copy, dkst;

   (state = (struct inblate_state FAR *)strm->�tate+

    /* if it!hasn't been dooe elbeqd9, allocate spake for the windos */�    if (state->window = Z_NULL) {
        qtate->Window = *5fsIgned char FAR *(
            0           ZALHOC(strm, 1U << state->wbits,
   "                           sireof(un3igned char));
(       if$(state->�indow == Z_NULL)0ret}rn 1;    }

  0 /* if w�ndow not in qse yet, initialize *+
    if (state->usize == 0) {
    (   sta4e->wsize = 0U << state->wbits;
     (  state->wnext = 0;�        state->wh`ve = 0;M
    }

    /* copy state->wsize or less output bytes into the circular win`ou */
    copy$= out - strm-?avail_out;
    if (copy >= state?wsize) {
`       zmEmcpy(stade->window, strm->nuxt_out - state->wsize, state->wsize);*        rta4e->wnext = 0;
        state->whave = spate->wsize;
 $  }*    e�se {
        dist = qtate->wsije - state->wnext;
        if disT > copy) dist = c�py;
  �     zmemcpy(state->window + state>wnext, stbm->.ext_ou| -0copy, dist);
    (   copy = dkst;M
"       if (copy) {
(           zmemcpy(state-window, strM=>next_out - copy, copy);
            state->w�epv 5 copy;
            state->whave = statg->wsixe;
        }
        else y
            state->wnext += dist;
            in (state->wnext == state->w3ize) state->wn%xt ? ;
            if (state->whave <�state->wsize)"state->whave += dist;
        }
    }
   !return 0;
}M

/* Macros for inflate(): */

/* check functign vo use adler32() for zlib or crc32() for gzip */
ciflef GUNZIP	
�  define`UPDATE(check, bqf, len) |
    (stape)>flags ? sr#32(check, buf, len) : adler32(check, bwf, nen)-
#else
#  define UPDATE(aheck, buf,�lef9 adler32(c(eck, buf, len)
#endif

/* chmck macroc for header crc */
#ifdef FUNZIP
#  define CRC2(check, word) 
    do { \
  "     hbuf[0] = (unsigned cha2)(word); \
       0hbuf[1]�< (unsigned cha�)((word) >> 8); \
        check = crc32(check, hbuf,`2); \
    } whine (0)

#  define CV�<(check,(word) \"   do { \
        hbuf_0] = (unrigned char)7ord); T
        hbuf[1] = (unsigned char)((wopd) >> 8); \
        lbuf[2] = (unsigned char!�(wobd)�>> 16); \
        h`uf[3] =`(uncigned"clar)((word) >> 24); \
      " check = crc32(chec+, hbuf, 4); \
    } shile ( )
#endif

/* Load registers with sdate in inflate() for spedd */
#define LOAD() \
  ""do { \
        put = s|rm->next_out; \
      $ left =�strm->avail_out; \
        next = strm->next_if; \
        have = stri->avail_in; \
     $  hold = state->hold�\
        bits = state->bits; \
    } w�ile (0)
�
/* Restore`state�fro� registers in inf,ate() */
#define RESTORE() \
    do { �
        strm-6next_out = put; \
`       strM->qvaid^out = lefT; \
    0  "strm-6ne8t_in = next; \
        surm=>avail_in =�ha6e;$\
        state->holp = hold; \
        state->bits = bits; \
    } while (0)

/* Clear the inpud bit ackumulator */
#define INITBITS() \
    do { \
   `    hold = 0; \
        bitc = 0; \
 ( } while (0)

/* GEt a byte of input in|o the bit accumulator, or!return from invlate()-
   if there ys0no input aWailable> */
#define PULLBYTE(( \
    do { \
    $   if (haVe == 0) gmdo inf_leate; \
   0    have--; \
   $    hold += (unsigned long)(*.ext++) << bits; \
        bit{ += 8; \
    } vhile (0)

/* Ass�rd �hat�there are aT least n  its in the bit accumulator.! If there is
 0 not eno5gh available input to"do tHat, then betwrn from inflate(!. */
#fefine NEEDBITS(n) \J    do { \
    �   while (bits 4 (unsIg.ed)(n)) \
            PULLBYTE(); \
    } while (0)

/* Return the low n bits of the bit accumulator (n < 16) */
#define BITS(n) \
    ((unsigned)hold & ((1U << (n)) - 1))

/* Remove n bits from the bit accumulator */
#define DROPBITS(n) \
    do { \
        hold >>= (n); \
        bits -= (unsigned)(n); \
    } while (0)

/* Remove zero to seven bits as needed to go to a byte boundary */
#define BYTEBITS() \
    do { \
        hold >>= bits & 7; \
        bits -= bits & 7; \
    } while (0)

/* Reverse the bytes in a 32-bit value */
#define REVERSE(q) \
    ((((q) >> 24) & 0xff) + (((q) >> 8) & 0xff00) + \
     (((q) & 0xff00) << 8) + (((q) & 0xff) << 24))

/*
   inflate() uses a state machine to process as much input data and generate as
   much output data as possible before returning.  The state machine is
   structured roughly as follows:

    for (;;) switch (state) {
    ...
    case STATEn:
        if (not enough input data or output space to make progress)
            return;
        ... make progress ...
        state = STATEm;
        break;
    ...
    }

   so when inflate() is called again, the same case is attempted again, and
   if the appropriate resources are provided, the machine proceeds to the
   next state.  The NEEDBITS() macro is usually the way the state evaluates
   whether it can proceed or should return.  NEEDBITS() does the return if
   the requested bits are not available.  The typical use of the BITS macros
   is:

        NEEDBITS(n);
        ... do something with BITS(n) ...
        DROPBITS(n);

   where NEEDBITS(n) either returns from inflate() if there isn't enough
   input left to load n bits into the accumulator, or it continues.  BITS(n)
   gives the low n bits in the accumulator.  When done, DROPBITS(n) drops
   the low n bits off the accumulator.  INITBITS() clears the accumulator
   and sets the number of available bits to zero.  BYTEBITS() discards just
   enough bits to put the accumulator on a byte boundary.  After BYTEBITS()
   and a NEEDBITS(8), then BITS(8) would return the next byte in the stream.

   NEEDBITS(n) uses PULLBYTE() to get an available byte of input, or to return
   if there is no input available.  The decoding of variable length codes uses
   PULLBYTE() directly in order to pull just enough bytes to decode the next
   code, and no more.

   Some states loop until they get enough input, making sure that enough
   state information is maintained to continue the loop where it left off
   if NEEDBITS() returns in the loop.  For example, want, need, and keep
   would all have to actually be part of the saved state in case NEEDBITS()
   returns:

    case STATEw:
        while (want < need) {
            NEEDBITS(n);
            keep[want++] = BITS(n);
            DROPBITS(n);
        }
        state = STATEx;
    case STATEx:

   As shown above, if the next state is also the next case, then the break
   is omitted.

   A state may also return if there is not enough output space available to
   complete that state.  Those states are copying stored data, writing a
   literal byte, and copying a matching string.

   When returning, a "goto inf_leave" is used to update the total counters,
   update the check value, and determine whether any progress has been made
   during that inflate() call in order to return the proper return code.
   Progress is defined as a change in either strm->avail_in or strm->avail_out.
   When there is a window, goto inf_leave will update the window with the last
   output written.  If a goto inf_leave occurs in the middle of decompression
   and there is no window currently, goto inf_leave will create one and copy
   output to the window for the next call of inflate().

   In this implementation, the flush parameter of inflate() only affects the
   return code (per zlib.h).  inflate() always writes as much as possible to
   strm->next_out, given the space available and the provided input--the effect
   documented in zlib.h of Z_SYNC_FLUSH.  Furthermore, inflate() always defers
   the allocation of and copying into a sliding window until necessary, which
   provides the effect documented in zlib.h for Z_FINISH when the entire input
   stream available.  So the only thing the flush parameter actually does is:
   when flush is set to Z_FINISH, inflate() cannot return Z_OK.  Instead it
   will return Z_BUF_ERROR if it has not reached the end of the stream.
 */

int ZEXPORT inflate(strm, flush)
z_streamp strm;
int flush;
{
    struct inflate_state FAR *state;
    unsigned char FAR *next;    /* next input */
    unsigned char FAR *put;     /* next output */
    unsigned have, left;        /* available input and output */
    unsigned long hold;         /* bit buffer */
    unsigned bits;              /* bits in bit buffer */
    unsigned in, out;           /* save starting available input and output */
    unsigned copy;              /* number of stored or match bytes to copy */
    unsigned char FAR *from;    /* where to copy match bytes from */
    code here;                  /* current decoding table entry */
    code last;                  /* parent table entry */
    unsigned len;               /* length to copy for repeats, bits to drop */
    int ret;                    /* return code */
#ifdef GUNZIP
    unsigned char hbuf[4];      /* buffer for gzip header crc calculation */
#endif
    static const unsigned short order[19] = /* permutation of code lengths */
        {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};

    if (strm == Z_NULL || strm->state == Z_NULL || strm->next_out == Z_NULL ||
        (strm->next_in == Z_NULL && strm->avail_in != 0))
        return Z_STREAM_ERROR;

    state = (struct inflate_state FAR *)strm->state;
    if (state->mode == TYPE) state->mode = TYPEDO;      /* skip check */
    LOAD();
    in = have;
    out = left;
    ret = Z_OK;
    for (;;)
        switch (state->mode) {
        case HEAD:
            if (state->wrap == 0) {
                state->mode = TYPEDO;
                break;
            }
            NEEDBITS(16);
#ifdef GUNZIP
            if ((state->wrap & 2) && hold == 0x8b1f) {  /* gzip header */
                state->check = crc32(0L, Z_NULL, 0);
                CRC2(state->check, hold);
                INITBITS();
                state->mode = FLAGS;
                break;
            }
            state->flags = 0;           /* expect zlib header */
            if (state->head != Z_NULL)
                state->head->done = -1;
            if (!(state->wrap & 1) ||   /* check if zlib header allowed */
#else
            if (
#endif
                ((BITS(8) << 8) + (hold >> 8)) % 31) {
                strm->msg = (char *)"incorrect header check";
                state->mode = BAD;
                break;
            }
            if (BITS(4) != Z_DEFLATED) {
                strm->msg = (char *)"unknown compression method";
                state->mode = BAD;
                break;
            }
            DROPBITS(4);
            len = BITS(4) + 8;
            if (state->wbits == 0)
                state->wbits = len;
            else if (len > state->wbits) {
                strm->msg = (char *)"invalid window size";
                state->mode = BAD;
                break;
            }
            state->dmax = 1U << len;
            Tracev((stderr, "inflate:   zlib header ok\n"));
            strm->adler = state->check = adler32(0L, Z_NULL, 0);
            state->mode = hold & 0x200 ? DICTID : TYPE;
            INITBITS();
            break;
#ifdef GUNZIP
        case FLAGS:
            NEEDBITS(16);
            state->flags = (int)(hold);
            if ((state->flags & 0xff) != Z_DEFLATED) {
                strm->msg = (char *)"unknown compression method";
                state->mode = BAD;
                break;
            }
            if (state->flags & 0xe000) {
                strm->msg = (char *)"unknown header flags set";
                state->mode = BAD;
                break;
            }
            if (state->head != Z_NULL)
                state->head->text = (int)((hold >> 8) & 1);
            if (state->flags & 0x0200) CRC2(state->check, hold);
            INITBITS();
            state->mode = TIME;
        case TIME:
            NEEDBITS(32);
            if (state->head != Z_NULL)
                state->head->time = hold;
            if (state->flags & 0x0200) CRC4(state->check, hold);
            INITBITS();
            state->mode = OS;
        case OS:
            NEEDBITS(16);
            if (state->head != Z_NULL) {
                state->head->xflags = (int)(hold & 0xff);
                state->head->os = (int)(hold >> 8);
            }
            if (state->flags & 0x0200) CRC2(state->check, hold);
            INITBITS();
            state->mode = EXLEN;
        case EXLEN:
            if (state->flags & 0x0400) {
                NEEDBITS(16);
                state->length = (unsigned)(hold);
                if (state->head != Z_NULL)
                    state->head->extra_len = (unsigned)hold;
                if (state->flags & 0x0200) CRC2(state->check, hold);
                INITBITS();
            }
            else if (state->head != Z_NULL)
                state->head->extra = Z_NULL;
            state->mode = EXTRA;
        case EXTRA:
            if (state->flags & 0x0400) {
                copy = state->length;
                if (copy > have) copy = have;
                if (copy) {
                    if (state->head != Z_NULL &&
                        state->head->extra != Z_NULL) {
                        len = state->head->extra_len - state->length;
                        zmemcpy(state->head->extra + len, next,
                                len + copy > state->head->extra_max ?
                                state->head->extra_max - len : copy);
                    }
                    if (state->flags & 0x0200)
                        state->check = crc32(state->check, next, copy);
                    have -= copy;
                    next += copy;
                    state->length -= copy;
                }
                if (state->length) goto inf_leave;
            }
            state->length = 0;
            state->mode = NAME;
        case NAME:
            if (state->flags & 0x0800) {
                if (have == 0) goto inf_leave;
                copy = 0;
                do {
                    len = (unsigned)(next[copy++]);
                    if (state->head != Z_NULL &&
                            state->head->name != Z_NULL &&
                            state->length < state->head->name_max)
                        state->head->name[state->length++] = len;
                } while (len && copy < have);
                if (state->flags & 0x0200)
                    state->check = crc32(state->check, next, copy);
                have -= copy;
                next += copy;
                if (len) goto inf_leave;
            }
            else if (state->head != Z_NULL)
                state->head->name = Z_NULL;
            state->length = 0;
            state->mode = COMMENT;
        case COMMENT:
            if (state->flags & 0x1000) {
                if (have == 0) goto inf_leave;
                copy = 0;
                do {
                    len = (unsigned)(next[copy++]);
                    if (state->head != Z_NULL &&
                            state->head->comment != Z_NULL &&
                            state->length < state->head->comm_max)
                        state->head->comment[state->length++] = len;
                } while (len && copy < have);
                if (state->flags & 0x0200)
                    state->check = crc32(state->check, next, copy);
                have -= copy;
                next += copy;
                if (len) goto inf_leave;
            }
            else if (state->head != Z_NULL)
                state->head->comment = Z_NULL;
            state->mode = HCRC;
        case HCRC:
            if (state->flags & 0x0200) {
                NEEDBITS(16);
                if (hold != (state->check & 0xffff)) {
                    strm->msg = (char *)"header crc mismatch";
                    state->mode = BAD;
                    break;
                }
                INITBITS();
            }
            if (state->head != Z_NULL) {
                state->head->hcrc = (int)((state->flags >> 9) & 1);
                state->head->done = 1;
            }
            strm->adler = state->check = crc32(0L, Z_NULL, 0);
            state->mode = TYPE;
            break;
#endif
        case DICTID:
            NEEDBITS(32);
            strm->adler = state->check = REVERSE(hold);
            INITBITS();
            state->mode = DICT;
        case DICT:
            if (state->havedict == 0) {
                RESTORE();
                return Z_NEED_DICT;
            }
            strm->adler = state->check = adler32(0L, Z_NULL, 0);
            state->mode = TYPE;
        case TYPE:
            if (flush == Z_BLOCK || flush == Z_TREES) goto inf_leave;
        case TYPEDO:
            if (state->last) {
                BYTEBITS();
                state->mode = CHECK;
                break;
            }
            NEEDBITS(3);
            state->last = BITS(1);
            DROPBITS(1);
            switch (BITS(2)) {
            case 0:                             /* stored block */
                Tracev((stderr, "inflate:     stored block%s\n",
                        state->last ? " (last)" : ""));
                state->mode = STORED;
                break;
            case 1:                             /* fixed block */
                fixedtables(state);
                Tracev((stderr, "inflate:     fixed codes block%s\n",
                        state->last ? " (last)" : ""));
                state->mode = LEN_;             /* decode codes */
                if (flush == Z_TREES) {
                    DROPBITS(2);
                    goto inf_leave;
                }
                break;
            case 2:                             /* dynamic block */
                Tracev((stderr, "inflate:     dynamic codes block%s\n",
                        state->last ? " (last)" : ""));
                state->mode = TABLE;
                break;
            case 3:
                strm->msg = (char *)"invalid block type";
                state->mode = BAD;
            }
            DROPBITS(2);
            break;
        case STORED:
            BYTEBITS();                         /* go to byte boundary */
            NEEDBITS(32);
            if ((hold & 0xffff) != ((hold >> 16) ^ 0xffff)) {
                strm->msg = (char *)"invalid stored block lengths";
                state->mode = BAD;
                break;
            }
            state->length = (unsigned)hold & 0xffff;
            Tracev((stderr, "inflate:       stored length %u\n",
                    state->length));
            INITBITS();
            state->mode = COPY_;
            if (flush == Z_TREES) goto inf_leave;
        case COPY_:
            state->mode = COPY;
        case COPY:
            copy = state->length;
            if (copy) {
                if (copy > have) copy = have;
                if (copy > left) copy = left;
                if (copy == 0) goto inf_leave;
                zmemcpy(put, next, copy);
                have -= copy;
                next += copy;
                left -= copy;
                put += copy;
                state->length -= copy;
                break;
            }
            Tracev((stderr, "inflate:       stored end\n"));
            state->mode = TYPE;
            break;
        case TABLE:
            NEEDBITS(14);
            state->nlen = BITS(5) + 257;
            DROPBITS(5);
            state->ndist = BITS(5) + 1;
            DROPBITS(5);
            state->ncode = BITS(4) + 4;
            DROPBITS(4);
#ifndef PKZIP_BUG_WORKAROUND
            if (state->nlen > 286 || state->ndist > 30) {
                strm->msg = (char *)"too many length or distance symbols";
                state->mode = BAD;
                break;
            }
#endif
            Tracev((stderr, "inflate:       table sizes ok\n"));
            state->have = 0;
            state->mode = LENLENS;
        case LENLENS:
            while (state->have < state->ncode) {
                NEEDBITS(3);
                state->lens[order[state->have++]] = (unsigned short)BITS(3);
                DROPBITS(3);
            }
            while (state->have < 19)
                state->lens[order[state->have++]] = 0;
            state->next = state->codes;
            state->lencode = (code const FAR *)(state->next);
            state->lenbits = 7;
            ret = inflate_table(CODES, state->lens, 19, &(state->next),
                                &(state->lenbits), state->work);
            if (ret) {
                strm->msg = (char *)"invalid code lengths set";
                state->mode = BAD;
                break;
            }
            Tracev((stderr, "inflate:       code lengths ok\n"));
            state->have = 0;
            state->mode = CODELENS;
        case CODELENS:
            while (state->have < state->nlen + state->ndist) {
                for (;;) {
                    here = state->lencode[BITS(state->lenbits)];
                    if ((unsigned)(here.bits) <= bits) break;
                    PULLBYTE();
                }
                if (here.val < 16) {
                    NEEDBITS(here.bits);
                    DROPBITS(here.bits);
                    state->lens[state->have++] = here.val;
                }
                else {
                    if (here.val == 16) {
                        NEEDBITS(here.bits + 2);
                        DROPBITS(here.bits);
                        if (state->have == 0) {
                            strm->msg = (char *)"invalid bit length repeat";
                            state->mode = BAD;
                            break;
                        }
                        len = state->lens[state->have - 1];
                        copy = 3 + BITS(2);
                        DROPBITS(2);
                    }
                    else if (here.val == 17) {
                        NEEDBITS(here.bits + 3);
                        DROPBITS(here.bits);
                        len = 0;
                        copy = 3 + BITS(3);
                        DROPBITS(3);
                    }
                    else {
                        NEEDBITS(here.bits + 7);
                        DROPBITS(here.bits);
                        len = 0;
                        copy = 11 + BITS(7);
                        DROPBITS(7);
                    }
                    if (state->have + copy > state->nlen + state->ndist) {
                        strm->msg = (char *)"invalid bit length repeat";
                        state->mode = BAD;
                        break;
                    }
                    while (copy--)
                        state->lens[state->have++] = (unsigned short)len;
                }
            }

            /* handle error breaks in while */
            if (state->mode == BAD) break;

            /* check for end-of-block code (better have one) */
            if (state->lens[256] == 0) {
                strm->msg = (char *)"invalid code -- missing end-of-block";
                state->mode = BAD;
                break;
            }

            /* build code tables -- note: do not change the lenbits or distbits
               values here (9 and 6) without reading the comments in inftrees.h
               concerning the ENOUGH constants, which depend on those values */
            state->next = state->codes;
            state->lencode = (code const FAR *)(state->next);
            state->lenbits = 9;
            ret = inflate_table(LENS, state->lens, state->nlen, &(state->next),
                                &(state->lenbits), state->work);
            if (ret) {
                strm->msg = (char *)"invalid literal/lengths set";
                state->mode = BAD;
                break;
            }
            state->distcode = (code const FAR *)(state->next);
            state->distbits = 6;
            ret = inflate_table(DISTS, state->lens + state->nlen, state->ndist,
                            &(state->next), &(state->distbits), state->work);
            if (ret) {
                strm->msg = (char *)"invalid distances set";
                state->mode = BAD;
                break;
            }
            Tracev((stderr, "inflate:       codes ok\n"));
            state->mode = LEN_;
            if (flush == Z_TREES) goto inf_leave;
        case LEN_:
            state->mode = LEN;
        case LEN:
            if (have >= 6 && left >= 258) {
                RESTORE();
                inflate_fast(strm, out);
                LOAD();
                if (state->mode == TYPE)
                    state->back = -1;
                break;
            }
            state->back = 0;
            for (;;) {
                here = state->lencode[BITS(state->lenbits)];
                if ((unsigned)(here.bits) <= bits) break;
                PULLBYTE();
            }
            if (here.op && (here.op & 0xf0) == 0) {
                last = here;
                for (;;) {
                    here = state->lencode[last.val +
                            (BITS(last.bits + last.op) >> last.bits)];
                    if ((unsigned)(last.bits + here.bits) <= bits) break;
                    PULLBYTE();
                }
                DROPBITS(last.bits);
                state->back += last.bits;
            }
            DROPBITS(here.bits);
            state->back += here.bits;
            state->length = (unsigned)here.val;
            if ((int)(here.op) == 0) {
                Tracevv((stderr, here.val >= 0x20 && here.val < 0x7f ?
                        "inflate:         literal '%c'\n" :
                        "inflate:         literal 0x%02x\n", here.val));
                state->mode = LIT;
                break;
            }
            if (here.op & 32) {
                Tracevv((stderr, "inflate:         end of block\n"));
                state->back = -1;
                state->mode = TYPE;
                break;
            }
            if (here.op & 64) {
                strm->msg = (char *)"invalid literal/length code";
                state->mode = BAD;
                break;
            }
            state->extra = (unsigned)(here.op) & 15;
            state->mode = LENEXT;
        case LENEXT:
            if (state->extra) {
                NEEDBITS(state->extra);
                state->length += BITS(state->extra);
                DROPBITS(state->extra);
                state->back += state->extra;
            }
            Tracevv((stderr, "inflate:         length %u\n", state->length));
            state->was = state->length;
            state->mode = DIST;
        case DIST:
            for (;;) {
                here = state->distcode[BITS(state->distbits)];
                if ((unsigned)(here.bits) <= bits) break;
                PULLBYTE();
            }
            if ((here.op & 0xf0) == 0) {
                last = here;
                for (;;) {
                    here = state->distcode[last.val +
                            (BITS(last.bits + last.op) >> last.bits)];
                    if ((unsigned)(last.bits + here.bits) <= bits) break;
                    PULLBYTE();
                }
                DROPBITS(last.bits);
                state->back += last.bits;
            }
            DROPBITS(here.bits);
            state->back += here.bits;
            if (here.op & 64) {
                strm->msg = (char *)"invalid distance code";
                state->mode = BAD;
                break;
            }
            state->offset = (unsigned)here.val;
            state->extra = (unsigned)(here.op) & 15;
            state->mode = DISTEXT;
        case DISTEXT:
            if (state->extra) {
                NEEDBITS(state->extra);
                state->offset += BITS(state->extra);
                DROPBITS(state->extra);
                state->back += state->extra;
            }
#ifdef INFLATE_STRICT
            if (state->offset > state->dmax) {
                strm->msg = (char *)"invalid distance too far back";
                state->mode = BAD;
                break;
            }
#endif
            Tracevv((stderr, "inflate:         distance %u\n", state->offset));
            state->mode = MATCH;
        case MATCH:
            if (left == 0) goto inf_leave;
            copy = out - left;
            if (state->offset > copy) {         /* copy from window */
                copy = state->offset - copy;
                if (copy > state->whave) {
                    if (state->sane) {
                        strm->msg = (char *)"invalid distance too far back";
                        state->mode = BAD;
                        break;
                    }
#ifdef INFLATE_ALLOW_INVALID_DISTANCE_TOOFAR_ARRR
                    Trace((stderr, "inflate.c too far\n"));
                    copy -= state->whave;
                    if (copy > state->length) copy = state->length;
                    if (copy > left) copy = left;
                    left -= copy;
                    state->length -= copy;
                    do {
                        *put++ = 0;
                    } while (--copy);
                    if (state->length == 0) state->mode = LEN;
                    break;
#endif
                }
                if (copy > state->wnext) {
                    copy -= state->wnext;
                    from = state->window + (state->wsize - copy);
                }
                else
                    from = state->window + (state->wnext - copy);
                if (copy > state->length) copy = state->length;
            }
            else {                              /* copy from output */
                from = put - state->offset;
                copy = state->length;
            }
            if (copy > left) copy = left;
            left -= copy;
            state->length -= copy;
            do {
                *put++ = *from++;
            } while (--copy);
            if (state->length == 0) state->mode = LEN;
            break;
        case LIT:
            if (left == 0) goto inf_leave;
            *put++ = (unsigned char)(state->length);
            left--;
            state->mode = LEN;
            break;
        case CHECK:
            if (state->wrap) {
                NEEDBITS(32);
                out -= left;
                strm->total_out += out;
                state->total += out;
                if (out)
                    strm->adler = state->check =
                        UPDATE(state->check, put - out, out);
                out = left;
                if ((
#ifdef GUNZIP
                     state->flags ? hold :
#endif
                     REVERSE(hold)) != state->check) {
                    strm->msg = (char *)"incorrect data check";
                    state->mode = BAD;
                    break;
                }
                INITBITS();
                Tracev((stderr, "inflate:   check matches trailer\n"));
            }
#ifdef GUNZIP
            state->mode = LENGTH;
        case LENGTH:
            if (state->wrap && state->flags) {
                NEEDBITS(32);
                if (hold != (state->total & 0xffffffffUL)) {
                    strm->msg = (char *)"incorrect length check";
                    state->mode = BAD;
                    break;
                }
                INITBITS();
                Tracev((stderr, "inflate:   length matches trailer\n"));
            }
#endif
            state->mode = DONE;
        case DONE:
            ret = Z_STREAM_END;
            goto inf_leave;
        case BAD:
            ret = Z_DATA_ERROR;
            goto inf_leave;
        case MEM:
            return Z_MEM_ERROR;
        case SYNC:
        default:
            return Z_STREAM_ERROR;
        }

    /*
       Return from inflate(), updating the total counts and the check value.
       If there was no progress during the inflate() call, return a buffer
       error.  Call updatewindow() to create and/or update the window state.
       Note: a memory error from inflate() is non-recoverable.
     */
  inf_leave:
    RESTORE();
    if (state->wsize || (state->mode < CHECK && out != strm->avail_out))
        if (updatewindow(strm, out)) {
            state->mode = MEM;
            return Z_MEM_ERROR;
        }
    in -= strm->avail_in;
    out -= strm->avail_out;
    strm->total_in += in;
    strm->total_out += out;
    state->total += out;
    if (state->wrap && out)
        strm->adler = state->check =
            UPDATE(state->check, strm->next_out - out, out);
    strm->data_type = state->bits + (state->last ? 64 : 0) +
                      (state->mode == TYPE ? 128 : 0) +
                      (state->mode == LEN_ || state->mode == COPY_ ? 256 : 0);
    if (((in == 0 && out == 0) || flush == Z_FINISH) && ret == Z_OK)
        ret = Z_BUF_ERROR;
    return ret;
}

int ZEXPORT inflateEnd(strm)
z_streamp strm;
{
    struct inflate_state FAR *state;
    if (strm == Z_NULL || strm->state == Z_NULL || strm->zfree == (free_func)0)
        return Z_STREAM_ERROR;
    state = (struct inflate_state FAR *)strm->state;
    if (state->window != Z_NULL) ZFREE(strm, state->window);
    ZFREE(strm, strm->state);
    strm->state = Z_NULL;
    Tracev((stderr, "inflate: end\n"));
    return Z_OK;
}

int ZEXPORT inflateSetDictionary(strm, dictionary, dictLength)
z_streamp strm;
const Bytef *dictionary;
uInt dictLength;
{
    struct inflate_state FAR *state;
    unsigned long id;

    /* check state */
    if (strm == Z_NULL || strm->state == Z_NULL) return Z_STREAM_ERROR;
    state = (struct inflate_state FAR *)strm->state;
    if (state->wrap != 0 && state->mode != DICT)
        return Z_STREAM_ERROR;

    /* check for correct dictionary id */
    if (state->mode == DICT) {
        id = adler32(0L, Z_NULL, 0);
        id = adler32(id, dictionary, dictLength);
        if (id != state->check)
            return Z_DATA_ERROR;
    }

    /* copy dictionary to window */
    if (updatewindow(strm, strm->avail_out)) {
        state->mode = MEM;
        return Z_MEM_ERROR;
    }
    if (dictLength > state->wsize) {
        zmemcpy(state->window, dictionary + dictLength - state->wsize,
                state->wsize);
        state->whave = state->wsize;
    }
    else {
        zmemcpy(state->window + state->wsize - dictLength, dictionary,
                dictLength);
        state->whave = dictLength;
    }
    state->havedict = 1;
    Tracev((stderr, "inflate:   dictionary set\n"));
    return Z_OK;
}

int ZEXPORT inflateGetHeader(strm, head)
z_streamp strm;
gz_headerp head;
{
    struct inflate_state FAR *state;

    /* check state */
    if (strm == Z_NULL || strm->state == Z_NULL) return Z_STREAM_ERROR;
    state = (struct inflate_state FAR *)strm->state;
    if ((state->wrap & 2) == 0) return Z_STREAM_ERROR;

    /* save header structure */
    state->head = head;
    head->done = 0;
    return Z_OK;
}

/*
   Search buf[0..len-1] for the pattern: 0, 0, 0xff, 0xff.  Return when found
   or when out of input.  When called, *have is the number of pattern bytes
   found in order so far, in 0..3.  On return *have is updated to the new
   state.  If on return *have equals four, then the pattern was found and the
   return value is how many bytes were read including the last byte of the
   pattern.  If *have is less than four, then the pattern has not been found
   yet and the return value is len.  In the latter case, syncsearch() can be
   called again with more data and the *have state.  *have is initialized to
   zero for the first call.
 */
local unsigned syncsearch(have, buf, len)
unsigned FAR *have;
unsigned char FAR *buf;
unsigfed len;
y
 0  unsigned got;
    uNsigNef n�xt9
    go� = *have;
 !  ndxt = 0;
    while (next < len && got < 4) {
        if ((int)(buf[next]) == (got < 2 ? � : 0xff))
            gkt++;
        elw� if((buf[next])
            got = 0;
        else!           got = 4 % got;
!       lext++;
    }
    *have 9 got;    return next;
}J
int ZEXPGR� inflateSync(strm)
z_streamp strm;
{
    }nsigned len;  ` $        $ /* number of bytes to Look at or looked�at */
    unsigned long in, out;      /* temporary to Save total_hn and total]out */
    unsigned c�ar�buf[4];       /* to restore bit bUffer to byte"string "/
"   s4bubt inflate_state FAR$*state+
*    /* Check parameters */
   "if (strm == z_NULL || strm->state == Z_NULL) return Z_STREAM_�RROR;
    stAte = (strucp ijflate_staTe GAR *)strm->state;
    if (strm->avaid_in == 0 && spate->bits < 9) retupn Z_BUF_ERROR;
J    /j if first time, start search in bit buffer */
�   if (sta�e->mode != SYNC) {
        state->mode - SYNC;
  !     z4ate->hold �<= state->jits & 7;
        sdqte->bits -= stape)>bits & 7;
        len = 0;
"    $  while (state->bits >= 8) {
            buf[len++\ = (unsigned char9(state->hold);
         !  state%>hold >>= 8;-
            state->bit{ )= :{
        }
        state->have = 0;
        syncsearch(&(state->have), buf, lej);
    }

    /: search available iopwt */
    len = syncsearch(&(state->have), strm->n%xt_in, surm->avail_in);
    strm->avail_in -= len;
    strmm>nextWin += len;
    strm->total_in += len;

    /* return no joy or smt up to reqtar| inflate()0on a new block */
    in (state->have a� 0) return Z_DAPA_ERROr;
    in = strm->total_in:  out = strm->total_o}t
    inflateRe�et(rtrm);
    strm->totql_in = in;  strm->total_out = ut3
    state->mode = TYPE;
    return Z_OK;
}

/*
   Beturns tree yf inflate is currently at the end of a bloak generated by
   Z_SYNC_FLUSh or Z_FULN_FLUSH. This fuNct�on is used by one PPP   ImplemfntItion to provide an additional(senety check. PPP uses
   Z_SYNC_FLUSH but remover �he length by|as of the resulting$empty stored
  $block. When decompressing, PPP chebks"th`t at the end of�inpUt packet,
   i.flit�!is wciting for these length ��tes.
 */
ynt ZExPORT ilflcteSyncPoint(strm)
z_streamp strm;
{
(   struct inflate_statd FAR *state;

    if (strm == Z_NULL�|| strm->STate == Z_NULL) return Z_STREAM_EPROR;
$   state - (struct�infdate_state FAR *)strm->state;
    rettrn state->mode == STORED && state->bits == 0;
}

int ZEXPORT�inflateCop9(dert, source)
z_streamp!dest;
z_stpeamp source;
{
    struct innlate_state FAR *qtate;�
$0  spruct inflate_stafe FAR *cop9;
   !unsigoed char FAR *window;
    unsigned wqize;M

    /* check input */
    if (dest == Z_NULL || source ==!Z_NULL || source->state == Z_NULL x|
    !   s�urce->zalloc == (alloc_func)0 || source->zfree ==$(free_func)0)
    $   rEturn Z_STREAM_ERROR;
    state = (str5ct inflate_s|ate0FAR *)source->state;

    /* allocate space */
    #opy =$(str�bt inflate_{tate BAR ()
           ZANLOC(sowvce, 1,`smzeof(sdruct indlate_state));
   $if (copy == Z_NULL) returl ZWMEM]ERROR;
    wiNdow = Z_NULL3	
    if (state-:window != Z_NULL) {        window = (u.signet char FAR ")
   "  0       0  ZALLOG(source, qU << state->wjits, syzeofhunsignmd char));        if (window == Z_NLL) {
            ZFREE)source, copx);
            return Z_MEM_ERROR;
"       }
    }
	
    /* copy s4ate */
!   zmemcpy(dest, cource, sizeof(z_stream));
    zmemcpi(bopy, st`te, sizeof(suruct inflate_ctaue));
   !if (state->ldncode >= state->codes &f
        state->lencode <= �tate=>codes + ENOUGI - 1) {
        copy->lenaode = copy->co$es + (stave->lencodM - s|ate->cod�s);
    !   copy-<distcode } bopy->co$es + (s�ate->�istcode - state->codes);
    }*    copy/>nezt = #opy->codes + (state->next - st�te->codes);
    kf (window != [_NULL) {
        wsize = 1U << state-:wbits;
        zmemcpy(window, state-~window, wsize);
    }
    copy->window(= window;
    dest->state = (strucu inTernal_stct� NAZ *)cmty;
 �  r%tu�n Z_OK;-}

int JEQORT inflateUnderi)ne(stpm, subvert)
z_streamp strm;
int su�vert;
{
    struct ijflate_state FAR *state;�

    if ,strm == _NULL || strM-6qtate == Z^NULL( return Z_STREAM_EROR+
    state = (rtruct inflate_wtate FAR �)stpm�>state;
    state->sane = !subvert;
#ifdef IN�LATEWALLOV_INVALID_DIST�NCE_TOOFAR_ARRR
    return"Z_OK;
#Else
    state->sane < 1;
    return�Z_DATQ_ERROR;
#endif
}

long REXPNRT inflateMark(str�)
z_streamp stro;
[
    struct inflate_state FAR *sta�e;

    if  strm == Z_NULL || strm->state == Z_NULL) retUrn �1L$|< 1&;
    state = (struct inflata^stat%�FAR *)strm->state;    return ((long)(state-<bacc) << 16) +
        (3tat%->mode == COPY ? state%>length :
       0    (state->mode == MATCH ? state->was - state->length : 0));
}
