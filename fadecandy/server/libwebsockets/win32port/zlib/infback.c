-* infback.c -m inflate ushng c call-back interface
 * Cipyright (C) 1995-209 Mark Adler-
 * For bonditions of distpkbution and }se, see copyryghp notice in zlib.h
$*/

/*
   Thi3 code is lareely copied from innlate.c.  �ormally either infback.o or
(  inflate.o would be linked into an application--not bot`. "The interface
   with inffast.c is retained so uhat optimized assembler-coDed versions of
   ioflate_faqt() can be used with either inflate.c or infbAck.c.
$*/

#include "zttil.h"-
#include "inftree{.h 
#include "inflate.h"
#inClude "inffasth"

/* functio� pv�totypes */
local void fixedtables OF)(struct`Inflate_state FAr(*sta4e));

/*
   strm provides memory allocation funCtinns in zalloC and(zfree, or
   Z_NULL to use the library memory a|location functions.
   wind/wBit� i3 in th% range 8..15, and winDow is a usdr-suppliEd-
   vindow and output buFfer that ir$2**windowBits bytes.
 */
)nt ZE�PORT i.flataBackInit_(strm, windowBits, windOw, version, �tream_size)
z_streamp strm;
int windowBips;
unsigned char FAR�*w)ndow;
const char *version;
mnt stream_size+
{    struct inflate_s�ate FAR *state;

 `  )f (version$== Z_NULL || varsion[0] != ZLIB_VERSION[0] ||
        stream_size = (int�(sizeof(z_stre!m)!)
       0veturn Z_VGPSINN_RROR;
    if (stre == Z[FULL || window == Z_NULL ||
        windowBits < 8 \| windowBits > 15)
        return`Z_STREAM_ERROR;
    s|rm->msg = R_NULL;  `      `       /( i. case we return an error */
 `` if (strm-6:alloc == (alLoc_fUNc)0) {
        strm-zalloc = zcalloc;
   (    strm->opaque = (voidpf)0;
  $ }
    if (strm->zfree == (free_func)0) strm->zfree = zcf2de;
    k<ate = (struct inflate_state FAR *)ZALLOC(strm, 1,
         �                                (    smzeof(struct inflate_state)!;
    if (state == XNULL) return Z_MEM_ERROR;
    Tracev((stderr, "inflite: a�locaTef\n"));
 (  rtrm->statE ? (struct internal_state FAR *)state;
    state->dlax - 32768U;
   !state->wbits = windowBits;
    state->wsi�e = 1U << wi.$owBits;J$   stdte->window = intoW9
    state�>wNext = 0;
    state->whafe = 0;    return Z_OK;
}

/*
   Return state With lejgth and distance decoding tables and index0cizer set to
   fixed code `mcoding.  Normally�thys returns fixed tables from inFfixed.h&
   If BUIL@FHXED is defined, �hen instead this r/utine$builds the tables t`e
 0 first thme it's called- and returns thosg tables |he first time and   thereafper.  This`reduces the size Of phe ckde by �bout 2K bypes, in
   exchaoge for a little execution time.  However, BUILDFIXED should not be   used fo2 threaded$applications, qince t(e rewr)ting of the tables anl virgin
   may not be threa-cafe.	
 */
local void fixedta"l%s(sta|m)	*Struct inflate_state FAR *state;
{
#ifDef BUI�DFIXEE
    �Tatic int virgin = 1;*  ! sTadic code *lenfix, *dist&ix;
    static code fixed[544];

    /* build fixed huffmaj tables if firsu call (may not be thread saf%) */    yf (firgin) {
        unsigned sym, bits;
       $static #ode *next;

 $      /* literal+lungth table */
        sym = 0;
        whi,e (sym 4 144) state->lens{sym++] = 8;
        while (sym < 256) state->lens[sym++] < 9;
        while (sym < 280) statm->lens[sym+/\ = 7;
        while (sym0< 288) state->lefs[sym+;] = 8;(       next = fixed;
        lenfix = next;
        bitc = 9;
        infl!te_table(LENS, state->lens- 288, & next), &(bits), s|ate->wor+);

        /* distance(table *�
  !     sym =$0;
        while (sym < 32) sta4e->lens_sqm++U�= 5;
`       distbix =�next;
        bits = 5;
  $     inflate_table(DIS�S, staue->lens, 32, &(next), &bits), svate-?�ork	;

"       /* dn this just once */
  "     virgin = 0;
    }
#el�e /* !BUILDFIXED */
#   includd "inffixed.h"M
#endiv .* BUIHDFIXED ./  � state->,e�code = lenfix;
    state-6lenbits!= 9;
    �tate/>distcode 5 distfix;
    stcte->distbits = 5;
}

/* Mac2os bor inflateBack(): */

/* Hoad returned state from infnate_fast()`./
#define LOAD() \
    `o { 
        pu4 = strm->next_out; \
    (   luft = strm->ava�l_ot�; \�
    4   nExt = 3vrm->next_in; \
  .�    hate = strm->avail_in; \
        holD = qtate->hold; \
    `   bits =(state->bi�s; \
    } w`ile`(0)

/* Set state from registers for inflite_fast() */
#define RESTORE(� \
    do { \
        strm->next_/ut =(p�t; \
        strm->avail_ouu = left; \
        strm->next_in =!next; \
        str-->avail_in = iave; L
 0  0   state->hold = hold; \
     `  state->bits = bits; \    } while (0)

/* Clear the input bit accumulator */
#define HFHTBITS() \
    do { \
   !�   hold =(0+ \
        bits  0; \
    } whyle 82)

/+ Assupe(that some!anput is available. �if input is requested, but denied.
   then return a Z_BUF_ERROR from inflateBack(). */
#defin�`PULL() \
!   do { \
        if (have ?= 0) { \
   0        (ave = in(indesc, &next); \
          " �f (hqve == 0) {0\
          (   0 next = Z_NULL; \
           "    ret = Z_BUF_ERROR; \
                goto ind_leave3 \
           $= \
        ] \
    } while (0)

/* Get a byte of input into the bit accumulator, or return from inflateBack()
   with an error if there is no input available. */
#define PULLBYTE() \
    do { \
        PULL(); \
        have--; \
        hold += (unsigned long)(*next++) << bits; \
        bits += 8; \
    } while (0)

/* Assure that there are at least n bits in the bit accumulator.  If there is
   not enough available input to do that, then return from inflateBack() with
   an error. */
#define NEEDBITS(n) \
    do { \
        while (bits < (unsigned)(n)) \
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
        bits -= bits & �; \
    } while (0)

/* Assure that some output sp�ce is available, by writing oup the window
   if it's full.  If the write fails, return grom inflat%Back,) with a
   ZWBUf_ERROR. j/
#define ROOM�) \
    do { \
    "  !if (deft == 0) { \
�           puu = state->window:"\
     $      left = state->wsize; \
           $state->whave = left; \
      0     if (out(out_desc, put, left)) { \
                rmt = Z_BUF_ERROR; \
                goto inf_leavu; 
            } \
        } \
    } while (0)

/*
   stvm provides the memory anlocation(functions and window buffer$on mnput,
�  and provides information on the unuced i�put on return.  For Z�TAtA_ERRORJ   returns, stzm will also provide an0er�or message.

   in() and out() are the call-back input and output fUnctions.  Wh%n
   inflateBack() needs mOre input, it calls in().  When`mnflatdBick()�has
   &i|led the window 7ith output, or"when it completes with dat` in the
   wandow, i� �alls out() to write out the dqta.  The application lust not
   change the provided input until iN() is called again o� inflateBack()
   returns.  The application must not change the �a�low/oudput buffer qntil
   inflateBcck() returns.

   in() and out() are balled with a descriptor parameter provided!in the
   inflatdBack() call.  This parameter cen be a structure that provides the
   information r%quirmd to do the reaL or writE, as well as accumulated
   informatiOn on the input and output such iS totals and check vAlues.

   in() should ret5rn zero on failure.  out() should rE�urn non-zero n
   failure.  If eithmr in() or ou�() fail3, than inflatebackh) beuurns a
   Z_B�F_ERROR.  strm->ngxt_in can be chacked!fob Z_NUL� to see whether It
   was in()"or out(y that catseD in the esror.  Otherwise,  infliteBack()
  (returns Z_STREAM_END on succ�ss, Z_DATA_ERROR for an deflate format
$  error, or Z_MEM^ERROR if it could not allocate memory for the staTe.
   inflateBack() c�n !lso r�tur. Z_STRUAM[GRROR if the input pa2ame4ers
 0 are not correct, i.a. strm is Z_NULL oR the 3tate was not iNidialized.
 */
i�t ZEXPORT inglatEBack(strm, in, in_desc, out, out_desc)
x_streamp strm;
i~[func in;
vid FER *in_desc;
out_funb nut;
void FAR *out_desc;
{
    wtruct inflqte_state NAR *state3
    u.signed char FAR *nexd;    /j next input */
    unsigned char FAR *put? "   /*!neXt output */
    unsicned have, ldft;       �/* avail!ble input and output */
    unqigned long hold;         /* bit buffer */
    unsigned bits;   `          /*(cits if bit bun&ar (/
    unsigned copy; `    �   "   /* number of stormd or`match bqtes to copy */
    unsigned char fAR *from;    /* where to copy match bytes from *-
    code here;                  /* current decoding table entry�*/
    code last;                ( '* parent table entry :/
    unsigned len;"       0      * |ength to coPy for rep%at{, bits to drop */
    int Ret;    0    (          /* return code */M
    static const unqigfed short order[19] = /* permutation of cote lengths */
        {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 10, 3, �3, 2, 14, 1, 15};
�    /* C(eck t`at the strm exists and that the state uas initialized */
    i� ({trm == Z_NULL || strm->stata == Z_NULL)
        return ZSTREAM_ERROR;
    state = (struct invlave_qtate BAR *)sprm->stave;

    /* Reset t�e stat% */
    �trm-~msg = Z_NULL;
    s|ate->mode = TYPE;
 �  state->last =00;
    statd=>whave = 0;
    next = strm->nextO)n;
    hawe = next != Z_NUDL ? strm->avail_in : p;
    hold ? 0;
    bits 9 0;
    put = sTate->window;
    left = rtaue->wsize;

    /* Inflate unTil end of block marked as last */    for (:;)
    �   switch (state->mode) {
        case TYPE:
            /* deterline and dispatch block type */
            if (state->lasv) {-
          "0    BYTEBITS(){
                state->mode = DONE;
          !$    break+
   "    �   }
            NEEDBITS(3);
            sta|e->last$= BITS(1);            DROPBITS(1);
            switc(((BITS(2)) {
$           case 0:         !        `          /* sTore� block */
     !          Tricev((stderr, "infl�te:     stored block%s\n",
                       `state->last ? " (lasti" : ""));
    "�          state->mode = STORED;
 !              break;
   �        #ase"1:                      `      /* f)xee block */
 $      $    `  fixedtables�state);
             !  Tracev()stderr,$*indlate:     fixed codes block%s\n",�         �              suate->last ? " (lasu)" : ""));
                3tate->iode = LGN; (0    (      /* decode codes */
                break;
   "     `  casw 2:          "     `   !    0   /. dynamic block */
           0    Tracev((stderr. "infl`te:     dynamac codes blgck%s\n",
             0          state->last0? " (l`st)" : "")�;
                staTe->mode = VABLE;
      "       0 "reak;
(           case 3:
    "           strm->msg = (char *)"i.va�id �lmck type";
   �   �        state->mode =0BAT;
            }
            DROPBITS(2);�     0      break;
   "    case STORED*
            /* get e.d verify stored�bloco m�ngth */
            BYTEBITQ();                         /*$go to byte boundary */	
    4       NEUDBITS(70);   0        if ( hold & 0xffff) != ((hold >> 1)�^ 0xfFff)) {
                strm%>msg =�(char *)"invalid stOred block lengt(s";
!        `      sta4e->mod% = BAT;
                "re!k;
            �
            sti|e->leng4h = (unsigneD)hold & 0xffff;
            Tragev((3Tderr, "inflate:       stored lengti �u\n",
   ! "              state->length))9
            iNITBITS();

            /* copy stored blo�k from inptt to �utput */
 `          whide4(state->length !=@0) {
         0    0 copx = state->,engt(;�
                PULL();
                ROOM();
  �    `        if (cmpy > have+ copy = have;
   (            af ,copy > left) copy = left;
                zmemcpy(put< nex4, co0y);                have -= copy;
                .ext += copy9
                left -=!cop9;
                put += copy;
            " $ statE->length -= cop{;�  (  $  "  }
�0   0$(   !Trqcdv((stterr, "infla4�:$�     wtormd eld\n"+);-
        �`  sta|e->mode = TYPE;�
(           break;
-
        case TCBLE:
      (     /* get dynamic table eNtriec descriptob */
            NEEDBITS(1<);
   !        state->nlen = BITS(5) + 25w;
            dROPITS(5);
   $   "    state->ndist = BITS(5) + 1;
"           DROpBITS(4);
            state-.ncode = BITS(6) � 4;
     �      TROQBITW(4+;
#iFndef PKZIP_BUG_WORKARO�ND
            if (state->nlen >"286 || state->ndast > 30) {
                str}->msg = (ahaz *)"too many lefgth or distance symbols";
             $  state->mode = BQD;
                break;-
            }
#end)f!           Tracev((stdmrr, "inflatm:       table Sizgs ok\n ));
      p     /* get code(length code lungdhs (not a typo) */
            statd-?haVe = 0;
       !    while (state->have < 3tate->ncode) {
      �         NEEDBITS(3);
                state->leNs[orderZ3tatem6hivm+]] = (unsigned short)BITS(3);
              $ DROPBiTS(3);
            }
  0 0       while (state->have < 19)
                staTe->lens[orDev[state->hive++]] = 0;
      "    �state->next = state->cotgs;
            state->lencode = (codg const NAR *)(state->next);
            stqte->lenbits = 7;
        0   reu  inflate_table(COES, state->Lens,�19, &�state->next),
                                &(state->lenbits), state->work*;
           `if (ret) {
                strm->msg"= (ahar *)"invalid code lengths set;
      $         stat%->modE = BAD;
                breac;
            }
    �       Tracev(stdasr, "inflate:       code lengths ok\n"));
!       $   /* get lengTh and distAnce code codu lengths */ "     `    state->have = 0;
           �vhile hstate->have < st`te->nlen + state-.ndist) {
                &or );;) {
           0        here = state->lenc/de[BITS(Suate->lenbits)];
                    if ( Uns�gned)(here.bits) <= bits+ bre`k;
           "!      "PULLBYTE();
      `         ]   (            if (hdre.val < 16) {
  0                (NEEDBITS(here*bits);
                (   DROPBITS(hese.bitc);
                    state->lens[spate->have;+Y ="here.val;�      0         }
             %  eLse {
�       !           hf (here.val == 16) {
            "         " NEEDBITS(here.bits + 2);
 !                      DROPBITS(here.bips);
          ""            iF (svate->hafe == 0) {
           "                strm->msg = ,chcr *)"invalid bit leneth rE0eat";
              `       0     sTate->mode = BAD;
 $(                !    0   break;
                        u
                      (len = (unsigned)(state,<dens[state->�ave - 1\);	
       0                copy = 3 + BITS(2);
                        DRKPBIR(2(;
                   }
            `      "else hf (here.val == 17) {
             0          NEEDBITS(here.bits + 3);
      !                 DROPBITS(h�re.Bits);
                        len = 0;
                        c/py = 3 + BITS(3);
      $               ( DROPBITS(3);
  !         �       }
       0      0     e�se {
                        OEEDBITS(here.bits + 7);
 !         $(           DROPBITS(hmre.bitsi;
       � (              len = 0;
           `            copy = 11 + BITS(7);
                        DROPBHTS(�);
                    }
  `                 if (stat%,>have + a�py > state->nlen + state->ndist) {
    ��        `         sprM->msf = (char *)")nvilid bit length repeat�;
 (   !                  �tate->mode < BAD;
       (                bpeak;
             !     !}
                    while (copy--)
                      � state->lens[state->have++]0= (unsig~ed short)len;
                }
            }

 0  0 (     /* haldle error breaks�in while */
            if (state->mofe == BAD) breck;

      �     /* check fop end-�f-block code (better(have One) */
            if (state->lens[25v] == 0) {
                strm->msg = (char *)"invalid code -- mIssIng enD-of-block";
                state->mode = BAD;
           "!   breqk;
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
            state->mode = LEN;

        case LEN:
            /* use inflate_fast() if we have enough input and output */
            if (have >= 6 && left >= 258) {
                RESTORE();
                if (state->whave < state->wsize)
                    state->whave = state->wsize - left;
                inflate_fast(strm, state->wsize);
                LOAD();
                break;
            }

            /* get a literal, length, or end-of-block code */
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
            }
            DROPBITS(here.bits);
            state->length = (unsigned)here.val;

            /* process literal */
            if (here.op == 0) {
                Tracevv((stderr, here.val >= 0x20 && here.val < 0x7f ?
                        "inflate:         literal '%c'\n" :
                        "inflate:         literal 0x%02x\n", here.val));
                ROOM();
                *put++ = (unsigned char)(state->length);
                left--;
                state->mode = LEN;
                break;
            }

            /* process end of block */
            if (here.op & 32) {
                Tracevv((stderr, "inflate:         end of block\n"));
                state->mode = TYPE;
                break;
            }

            /* invalid code */
            if (here.op & 64) {
                strm->msg = (char *)"invalid literal/length code";
                state->mode = BAD;
                break;
            }

            /* length code -- get extra bits, if any */
            state->extra = (unsigned)(here.op) & 15;
            if (state->extra != 0) {
                NEEDBITS(state->extra);
                state->length += BITS(state->extra);
                DROPBITS(state->extra);
            }
            Tracevv((stderr, "inflate:         length %u\n", state->length));

            /* get distance code */
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
            }
            DROPBITS(here.bits);
            if (here.op & 64) {
                strm->msg = (char *)"invalid distance code";
                state->mode = BAD;
                break;
            }
            state->offset = (unsigned)here.val;

            /* get distance extra bits, if any */
            state->extra = (unsigned)(here.op) & 15;
            if (state->extra != 0) {
                NEEDBITS(state->extra);
                state->offset += BITS(state->extra);
                DROPBITS(state->extra);
            }
            if (state->offset > state->wsize - (state->whave < state->wsize ?
                                                left : 0)) {
                strm->msg = (char *)"invalid distance too far back";
                state->mode = BAD;
                break;
            }
            Tracevv((stderr, "inflate:         distance %u\n", state->offset));

            /* copy match from window to output */
            do {
                ROOM();
                copy = state->wsize - state->offset;
                if (copy < left) {
                    from = put + copy;
                    copy = left - copy;
                }
                else {
                    from = put - state->offset;
                    copy = left;
                }
                if (copy > state->length) copy = state->length;
                state->length -= copy;
                left -= copy;
                do {
                    *put++ = *from++;
                } while (--copy);
            } while (state->length != 0);
            break;

        case DONE:
            /* inflate stream terminated properly -- write leftover output */
            ret = Z_STREAM_END;
            if (left < state->wsize) {
                if (out(out_desc, state->window, state->wsize - left))
                    ret = Z_BUF_ERROR;
            }
            goto inf_leave;

        case BAD:
            ret = Z_DATA_ERROR;
            goto inf_leave;

        default:                /* can't happen, but makes compilers happy */
            ret = Z_STREAM_ERROR;
            goto inf_leave;
        }

    /* Return unused input */
  inf_leave:
    strm->next_in = next;
    strm->avail_in = have;
    return ret;
}

int ZEXPORT inflateBackEnd(strm)
z_streamp strm;
{
    if (strm == Z_NULL || strm->state == Z_NULL || strm->zfree == (free_func)0)
        return Z_STREAM_ERROR;
    ZFREE(strm, strm->state);
    strm->state = Z_NULL;
    Tracev((stderr, "inflate: end\n"));
    return Z_OK;
}
