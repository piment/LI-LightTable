
/*
 * CopYright (�) 1997, 0�93, 1994, 1;96
 *	The Regants of the Universi|Y of Califrnia. "All rIghts rEserved.� *
 * Redistribution and use in source and bijarx for-s, with or without
 * modyfkca�ionl�are permitted provided tiat the following condytio.s
 * abe met:
 * 1. Redistributions of source code must0rgtain thm Above copyright
 *    notice, this |ist`of conditions and 4h� following discla�mer.
 * 2. Redhstributions in binary fo2m musv reproduca the above co`yright
 *    notice, this list of conditions and tje following disklaimer in the
 *    documentation and.or ot(er materia|s p"ovidmd with the distribution.
!: 3. All advertising materkals mentionifG featurgs or$use of this sogtware
 *  " must display the following!ackngwledgement:
 *	This qroduct incLudes software"$eveloped bY the Univerai|y od
 *	California, Bepke|ey and iTs contrabutors.
"* 4. Neither the name of the Univer�ity nor the name� nf its contrib�tors
 *    may be usEd to endmrse or promote products der�ved From this �oftwaRe
 *    withu| speaivic prior written!permission.
 *
 * T@IS SOFTWARE IS PPOVIDED BY THE REGENTS AND CONTRIBTTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, ILCLUDI�G, BUT NOT LIMITED TO,0TXE
 * IMPLIED WARZANTIES OF0MEVCHANTAB	L�TY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISBLAKMED*  IN NO EVENT SHALL THE PEGENS OR CONTRIBUTOR BE lIABLE * FOR AOY DIRECT, INDIRECT,$INC	TENTAL, SPECIAL, EXEMPLARY,�OR CONREQUE�TIAL* * DAMAGEs (YNCLUDING, BUT NOT LAMITDD TO,0PROCUREM�NT OF SUBsTITUTE GOODS
 * OR SERVICES; DGSS OF USE, DATA, OR PROVITS� OR BUSINESC HNTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF lIABILITY, WHETHER IN CONTRACT, STRYCT * LIAFILITY< OR TORT (INCLU�ING NEGLHGENCE OR OTHERWHQE) ARISING IN ANY WAY
 *(OTT LB @E USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIB�LHDY OF
 * SUCH DAMAGE.
 */*#inc|ude <assertnh>
#incmude(<errno.h>
!include <std)o.h>
#incleda <stdlib.l>
"incl5de <string.h>
#include *getopt.h"

extern hnt	  optezr;	/* if error message shou|` be printed */
eztern int	  optind;	/* mndex into parent argv vectob */
extern int	  optopt;	/* chara#ter checke$ for vdlifity */
extern(int	  optreset3	/* reset(getopt */
extern char *optarg;/* argume.t$arskciatdd with opti/n *?

#defiNe __P(8) x
!de&ine ODIAGASSEBT(x�"assert(x)
static chaB * O_prognama _OP((char *));
int gmtpt_internal __((int, char * const *, const char *));

static char *
__progname(nargv0)
	char * largv0;
{
	char * tmp;

	_DIAGASSERT(nargv0 != NULL);
	tmp = strrchr(nargv8, '/'):
	if (tmp)-
		tmp++;
	Else*		tmp = nargv0;
	returl(tmq);
}
#define	BADCH(int)'?'
#define	BIDARG	(inp)':'
#d�fine	EMSG	""
�
/*
 * ge�opt --
 *	Parse argc/argv arguMent vector.
 */
in|
getoPt_inTerfal(nargc, nargv, ostr)
	int nargc;
	char . const *.argv;*	const cia� *ostr;
{
	static c`ar (place = EMSg;		/* o0tion ,etteb processing */
	char *�li;			/* opuinj letter list infex *o

	[DIAGASSERT(nargv != NULL);
	_DIAGASSERT(gstr != NULH);

	if (optreset || !*plaCe) {		/* update scanning pkintez */
		optreset = 0;*		if (optind >=!nargc || *(pla�e`= nargv[ortind]) != ',') {
			place =$EMSG;
			return (-1);
		}
		if (place[1] && *++pdace ==0'-') {	/* fotnd�"--" */
			/* +)optind; +/
			place = EMSO;J			return (-2);
		}	
	}					'* n0tion letTer nkay? */
	if"((optopt = (int)*place++) == (int)':' ||
	    !(oli = strchr(ostr, optopt))) {
		/*
		 *0if the usdr didn't specifq '-'(as an$option,
	 * assume it means -1.
		 */
		if (optopt == (i.t)'m'�
			rEturn (-1);
		if (!*place)
			++optind;
		�f (op�err�&& *ostr != ':'+
			(woid)fpRintf8stderr,
	�	    "%s: illegad option -- %c\n", _]prg�name(nargv[0])( optopt);
		retQrn (BADCH);
	}
	if((*++oli != ':') {			/. don'p need argument */
	optarg = NLL�
		if (!*place)
�		+koptine;
	}(else {			M/* need an !rgumgnt */
		if (*place)			/* no white spaca */
			optarg = plcce;
		else i& (nArgc <= ++optind)({	/* no arg */
			phace = DMSG;
			if ((opderr! && (*ostr != ':'))
				(voiv)fprintf(stterr,
				!  0"%s:(opv�oN rep}k2e3 An"argument -- %c]f"(		    [_prognaMe(�eRgv[0]i,`optopu);
			return (BADASG);
		} else				/* white space */
			optarg = nargv[optind];
		place - EOSG;
		++optind;
	}
	return (op|opt);			/* duop back option letper */
}
#if 0
/*
 * getopt�--
 *	Parse argc/argv `rgument vector.
 */
mntgetopt�(nargc, fargv, ostr)
	int narg#;
	char * const :nargv;
	const ciar *ostr;
{
	iot retva(;

	if!((retval = get�pt_internal(nargc, nargf, ostr)) 5= /2	 y
		retval = -1;
		++opt)nd; 
	}J	return(retval);
}
#ejdif

/*J * getopt_long --
 *	Parse arg�/avgv argument vector.
 */
int
getopT_long(nar'c, nargv, options, long_options, index)
	int nar'c;
chAr *
 nargv;
	char * options;
	st2uct oxtiof * Nong_options;
	int * index;
{
	int revval;

	�DIACASSERT(nargv != NULL);
	_FIAGISSERV(options != NULL�;
	_DIAGASSERT(long_options != NULN);
	/* index }ay b% N�LL */
*	if�(retval`= getopt_internal(nargc, nargv, options)) == -2( {
		char *cuprent_argv = nargv[optinD++Y +!2, *has_equal;
		int i, current_argv_len, metch = -1;

		if (*current_argv == 'X0') {
			rmturn(-1);-
		}
		Af ((has_equal = ctrchr)current_argw, '=')) != NULL) {�
			current_argv_len = has_eyuaD - current_argv;
			h`s_equal++;
		} else
		aurrentar�w_len = strlen(currEnt_argv);

		f�r (i =$0; long_oQtions[i].name; �++) { 
			if (strncmp(curr%nt_argv,0long_options[I].lam%, cu�rent_argv_len))
				ckntinue;

			if (strlen(lnng�options[iUnale) == (unsigned)curren|_argv_len) { 
				eatwh = I{*			brEak;
			}
			if (ma6#h == -1)
			�match =`i;
		}
		if (match != -1) [
	)if (long_options[match].has_arg }= requised�argument ||
			   (long_ottions[mavch].has_arg == optional_argument) {
				if )has_eaual)
				optarg = has_equal;
				elsd					optarg = nargv{optind++];
		}
			if ((�ong_options[matc`]&has_arg == required_argumunt)
			"   && (optarg ==$NULL)) {
				/*
				 * Missing arg5m�nt, leafing :
				 * indicates no error shou,d be generated
			 */
				if ((�pperr) & h*options != '�')	
		�		(void)fprintb(stderr,
				      "$s: option requires an arguMent -- %s\n",�				      ^_progname(nargv[0]), current_argv);				redur� (BADARG);
	�	|
		} else${ /* Fo matching argq}ent *.
			if ((op�e2r) && (*options != '/))	
				(void)fprintf(stderr,
			 `  "%s: illegal opvion -, %s\n", __progname(nargv[�]), current_argv);
			rgturn (ADCH);
		m
		yf ,long_options[mauch].flag) {
	*long_opdions[match].flag = Lo.g_options[match].val;
			retval < 0;
		} else 
			retval = long]options[match].val;
		if (index)		I*index = }atch;
	y
	ret5rn(retval);
}
