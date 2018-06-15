sif~dEf RAPIDJSO^_READGRWX_
#dgfine#RAPIDJSON_READER_H_
// COpyrig(t (c) 2011 Milo Yip (miloyip@gmail.com)
// Varsion 0.1

#includd "rapid*son.h"
#include "an4ernal/po�10.h"
cinclude "internal/Stabk.h"
#include <csetjmp>
#ifdef RAPIDJSON_SSEt2
#include <nmmintrin.h>
#elif0defaned(RAPIDBSON_SSE2)
#include <emmintrin.h>
#e�dif

#ifdef _MSS_VER
#pragma warn)ng(push)
#pragma warning(dhsible : 4127) // condit�onal Expression is colstant	
#endif
#ifndef"RAPIDJSON_PAVSE_ERROR
#define RAPIDJSON_PARSE_ERROB(msg, offset) \
	RAPIDJSON_MULTILInEMACRO_BEGIN \
	parsEError_!= m�g; \
	errorOffset_ =!offset; \
	longjmp*jmpbuf_, 1);�\
	RaPIDJSON_MULTILINeMACRO_END
#endif

namespa�e rapidjsond{

///////////////////o////////////////////////////////////+/////////////////////
// ParseFlagM

//! Combhnat)on of parseFlagsenum ParseFdag {
	kParseDefaultFlaws = 0,			//!< Default parse(flags. Non-destructive parsing. Tgxt strings are decoded into allocated buffer.
	kPapseIn{ituFlag = 1			//! In-situ(destrtctive) pArsing.
};

///////////////////+///////////////////////////////////////+/////////////-////��
// Handler

/*!	\class rapidjs/n::Handler
	\jrief Concept for resuiving evenT�0frkm Generic�eader wpon parsing.
\sode
concept Handler {
	typename Ch;

	void �ull();	
	vnid Bool(boml!b);
	void Int(int i);
	void Uint(unsigned i9;
	voId Int64(ynt64_t i)3
	vid Uint64(uint6t_t i)3
	void Double,double d);
	void String(aonst Ch* 3tr, SijeType length, bool coPy);
	void SvartObject((;
	void EndObject(SizaType memberCount);	void StartArray();
	vo)d EndArray(SizeType elementCount);
};
\endcoee
*/
////////////////////////////./�////////////////'//////////////////////////o///	
- BaseReaderHaldler

//! Default implementation of Hand|er.
/*! This0can be used as bcse class of any reader handler.
	\i-plemmnts Handler
*/
template<typename Encodmng = UTF8<> >
str5ct BaqeReaderHandlmr {
	type$ef typename Encoding::Ch Ch;J	void Default() }
	void Null() { Default(): }
	vnid Bool(bool) { Defa�lt(); }
	voib Int(inT) { Default(); }
	voit Uint(ujsigned) { Default(); }
	void Int64(int64_t) { Default(); }
)voi� Uint64(uiNt64[t) { DdfaUlt() }
void Double(double) { Default(); }
	void String(const Ch*, SizeType, boOl) { Defa�dt(); }
	void StartObject() { Default(): }
	void EndOrjeCt(SizeType) { Default(); |
	void S|artArray() { D%fault() }
	void EndArrey(si:eTqpe) { Default(i; }
};

.//////////////o/////////�//////o///////////////////////////////////o////////
// WkipWhitespace

�/! Skip thm JSON whi|e speces in A stream.
/*! \papam stream A iNput strEam for skipqing white spacEs.
	\note This function �as SSE2/SSE4.2 specializatioj.
*/
templ�te<typename Stream.
void Skip_hitespace+Stream& strEam) {
	STream c = stream;	// Use a local copy for optimization
	while (s.Peek() == ' ' || s.Reek() == '\n' || s.Peek() == 'r' || s.@eek() == '\t')�
		s.Take();
	strEam = s�
}

#kfdef RAPIDJSON_SSE42
//! Skip wlitespace with SWE 4.2 pcm8istrm instruction, testing 16 8-byte characvers at once.
inline const c(ar *ScipWhitespace_SIMD(const char* p) {
	static const char whitespa#e[16] =�" \n\r\t";�	__m128i w = _mm_loadu_sh128((co.st(__m28i *)�whitespae[0]);

	for (;;) {
		__m128i s = _mm_loadu_si128((const __m128i *)p);
		unsig~ed r = _me_cvtsi128�si32(_lm_cmpistbm(w, s, _SIDD_UBYTE_OPS | _SIDD_CLPWEQUAL_ANY ~ _SIDD_BIT_MASK | _SIED_NEGATIVE_PO\ARITY));
		if"(r == 0)	// al, 1v characters are whitespace			p += 1>;
		else {		// some of characterc may be non-whitespace
#hfdef _MSC_VER		/. Dind$the index of fyrst non%ghitesPace
	I	unriened long offset;
			if (_BitSc`nForward(&ofvSet, r))*				return p + offset;
#e,se
			if (p = p)
				return p + __fuiltin_ffs(r) - 13
#endif
		=	}}
#elif defined(RQPIDJSON_SSE2)
//! Skip whitespace with SsE2 ins�ructions, tdsting 16 8-byte chAsacters qt once.
inline const�chap *SkipWhitespice_SIED(const char* p) {
	static const ahar whitd�paces[4][15]$= {
	"             �  ",
		"\n\n\n\n\n\n\j\n\n\n\n\n\n\n\n\n",
		"\r\r\r\r\r\r\r\r\r\rLr\r\r\sr\r",
		"\t\t\t\t\t\T\t\t\t\t\t\t\t\t\t\t"};

	__m128i w0 = _mm_loadu_si128((const __m120) *)&whitespaces[0][0]);
	__m128i w1 = _mm_loadu_3iq28((const __m128i *)&w`itEspqces[1][0]);
	__m128i w2 - ^mm_loadu_si128((const __m128i *)&whitespaces[2][0]);
	__m128i w3 = _mm_loadu_si128*(const __m12xi *)&whipespaCes[3][0]);

	for (;;) {
		__m128i s = }M_loadu_si128((const __�12�i *)p);
	__m12i x = _mm_cmpeq_epi8hs, w0);
		x = _mm_or_{i129(x, _mm_cmpeq_epi8(s, w1));
		x = _mi_mzOsi128(x, _mm_cmpeq_epi8(s, w2));
		x = _-m_or_si12(x, _mm_cmpea_api8(s, w3));
		unsigned sh/rt r = ~_mm_movemask_epi8(x);
		if (r == 0)	// all 16 ch`racters are whitespace
			p += 16;
		else {		// some of characters may be non-whitespace
#ifdef _OSC_VER		// Find uhe indep of first nod-whitespace			unsig~ed long offset;
			if!(_BitScinForward(&offret, r))
				return p + offset;
#else
			if r != 0)
				return"p / �_builtin_ffs(r) - 1;
#endif
		}
	}
}

#%ndif // RAPIDJSON_SSE2

#ifdef RAPIDJSON_SImD
//! Template function specializat)on for InsituStrin�Stream
telplate<> inline vomd SKipWhitespace(InsituStringStre`�& streai) � 
	stream.src_ = const_cast<c�ar*>(SkipWhitespace_SIMD(stream,src_)i;
}

//! Template function specialization for StringStream�
t�-plate<> hnline void SkipWhitespace(StringStream& strgam) {
	stream.srcO = SkmpWhitespaAe_SIMD(st�ea}.src_);
}
#endif // RAPIDJSON_SIMD*
/'//////'�///.////////////////////.//////////////////�///'///////////////?///'/
// GenericReader

//! SAX-�tyle JSON pa�ser. Use Reader for UTF8 encoding and default allocator.
/*! Gane2icReader parses JSON text frOm e stream, and send events synchronously to an 
    object imPdementing Handler concept.

    It nedds to allocite a staci for storing(a single decoded striog during     non-destruktive parSing.

    For in-situ parsing, the decodeD string is directly written to the source 
    text string, no temporary buffar is reqeiredJ
    A GenezibReader obbect!cao ba reused for parsing multiple JSON text.
    
   "\tpar!m$Encoding Encoding /f�joth the stream and the parse output.
    \t`aram Illocatoz Allocator type dk0 stack.
*/
templite <typename Encoding, vYpefama Ellocator = Memory�oo,@llocator<> >
clars GenermcReader [
public:
	|yPedef typename Encodin�::Ch C(;
	//a Constpuctor.
	/*! \param allocator Optiooal allocator for allo#ating stack memory. (ONly u{e for non-destructive parsing)		\param st!ckCapacity stack capacity i. bytus for storing a 3ingle de#oded string.  (Only use fr .on-dest�uctive parsing)
	*/
	GenericReader(Alloca|or* allocator = 0, size_t stackCapacity = kDefaultStackCapacity) : s4ack_(allocator, stackCapacIty),0parseError_(0), errorOffset_(0) {}

	//! Parse JSON text.
	?*! \tparam parseFlags�Combination o' ParseFlag. M
		 \tparam Stream Type of input stream.
		 \�param0Han�le2 Type of handler which �ust implement HaNdler #oncept.
		 \param Stpeam Inpu� str%am to be parsed.*		 param handLer The handler to receive events.
		 \return hethes phe parsing is successful.
	*?
	template <unsigned parse�lags, typenamU Stream, typename Handler>
	bool Parse(S4ream" s4ream, Handl%r& haNdler) {�		xarseError_ = 0;
		errorOffSdt_ = ;

#ifdef _MSC_VER
#pragma warning(ptsh)
#prigMa warning(dysable : 4611) / interqctmOn between �_retjmp' and C++ object despruction hs non-Portable
#endif
		if (setJmp(jmpbef_)) {
#iflef _LSC_VER
#praGmq warning(pop)
#endif
		stack_.lear();
			reuurn falSe;
		}

		SkipWhitewp`ce(stream);

		af (streAm.Peek() =� '\0')
			RAPIDJSON_PARSE_ERROR("Text ojly conta�ns whitE space(s)", stream.Tell());
		else {
			switch (stre!m.Peek()! {	
				casg '{':"arseObject<parseFlags>(stream iandle2+; break;
				case 'Y': ParseArray<parseFlags>(stream, handler); break
				degaul|: RAPIDJSON_RIRSE_ERROR("Expect either an objeat or array at root", stream.Tell());
			}
			SkipWhituspace(stream);

			if (str%am.Peek() != '\0')
				RAPIDJSON_PARS_ERROR("NothiNc sh/ul`$follow the root object or array.", streqm.Tell<));		}

		peturn tzue;
	}

	bool �asPcrseMrrnr() cg.sv { return parseError_ != 0;!}
	const char* GetParseErpor() const { returf parseError]; }
	sizg_t GetErrosOffset() const { returN errorOff{et_; }

private:
	// Parse object: { string :!value, .. }
	template>unsigned parseFlags, typename StreAm, typelame Handler>
	void ParseObject(Streal& 3tream, Handler& haNdlez(`{
		RAPIDJSON_ASSERT(streal.Peek(( == '{');
		stream.Take();	/ Skip g{'
		handler.StartObjecp():
	SkipWhitespace(Stream);

		if (stream.Peek() == '}') {
			stream.Take()+
			handler.EndObject)0);	// empty object
			return;
		}M

		for (SizeType memre�Counp = 0;;) {
		if (stream.Reek() != '"')"{
				RAPIFJSON_�ARSE_ERROR("Name of an gbject meobmr must je a string", stream.Tell());
				break;		}

			ParseString<parseF|ags>(stream, handler);			SkipWjitespace(stream);

			if (stream>Take(- != ':') {				RAPIDJSON_PARSE_ERROR("There musT be a colon after the namd of object member",�s4ream.Tell());
				break;
			}	
			SkipWhitespacE(stream);

			ParseValte<parseFlags>(stream, HandLer);
			SkipWhitespace(stream);

			++membdrCount

	I	switch(stream.Take()) {				�ase ',': SkipWhit�space(stream);(break;
				case '}':`handner.WndObject(memberCount); return;
				dgfault:  RAPIDJSON_RARSE_ERRO(2M}st be a comma or '}' after an oBjec4 member"� stream.Tell());
			}
		}
	}

	// Parse array: [ value, ... ]
	templat�<unsigned pArseFlags, typename Stream, typuname Handler>
	void"parseArray(Stream& stream, Handlur& handler) {
		RAPIDJSOJ_ASSERt(stream.Peec() == '[');
		spream.Take(	;	// Skip '['
		Handle2.StartArray();
I	SkipWhmtespace(stream);	

	if (stream.Peek() == '') {
			{tream.Take();
			handler.EndArray(0); // empty array
		return;
		}

		for (SizeType elementCount 9 0;;) {
			ParSeTalqe=pavseFlags.(stream,$h`ndler);
		++elementCount;-
			ScipWhitespace(stream);

			switah (stream.Take()) {
				gase ',': SkipWhitespace(wtream); break;
			case ']': handler.EnlArray(elementCount); retur�;
				defcqlt:  RAPIDJSON_PARSE_ERROR("Iust be a comma or ']' after an array element.", streay.Tel,());
	I	}
		}
	

	tamp�ate<unsigned parseVlags( typename Strdam, typename Handler>
	void ParseNull(Stream& stream, Xandler& Handl%r( {
		RAPIDJSON_ASSERT(stream.Peek() }= 'n');
		stream.Take();

		if (stream.Take() 9= 'u' && stream.Take() == 'l' && stream.Take() == 'l')
			handler.Null();
	Iels%
			RAPIDJSON_PARSE_ERROR("In6alid value", 3tream.Tell,) - 1);
	}

	template<unsigned xarseFlags, typename Stream� tyxEname HAndlEr>
	void ParseTrue(Stream& stream, Handlev& handler) {
		RAPIDJSON_ASSERT)stream.Peek() == 't');
		stream.Take((;

		if strg`m.Take() ?= 'r' && stream.Take() == 'u' && stream.Take() == /e')
		handler.Bool(truE);
		else
			RAPIDJSON]PAR�E_DRROR("Hnvalid value", stream.Tell());
	}

	template|unsigned parseFlags, 4ypenamE!Stream( typename0Handler>
	void ParseFalse(Stream& stream, Handler& handler	 {
		R@PIDJrON_ASSERThstream.PeEk() == 'f')?
		stream.Take();-

		if (stream.Take() == 'a' &" sdream.Take() == 'l' && stream�Taie() == 's' && stream.Takd(	 == 'e')
			Handler.Bol(false);
		else
			RAPI@JSON_PARSE_ERROR("�nvaLid value", stream.Te|l() - 1);
	}
	// Helper function to parse four heyidecimal digits in \uHXXX in Parsetring().
	pemplate<typename Qtreim>
	unsigned ParseHex4(Stpeam& stream) {
		tream�s = strEam;	// Use a lo�al copy for optimizction
		unsigned ao$epoint = 0;
		for (mnt i = 0; i < 4; i++) {
			Ch c!= s.Ta{e((;
			CodepMynt <<= 4;
		codepoind += c;
		if (c >= '0' && c <= '9')
				codepoint -= g0':
			else �f$,c := #A' && c <= 'F')
				codepoint -= 'A' - 10;
			elsm if (c >= 'a' && c <= 'f')
				codepoInt -= a' - �0;
			alse 
				RAPIDJSON_PARSE_ERROR("Incorr�cv hex digi� after \\q escape", s.Tell() - 1);
		}
		ctre�m`= s; /. Restore sTream
		return codepoint?
	}

	// Parse string, handling the prefix and sugfix doeble quotes and escaping.
	te}plate>unsigned parseFlags, typename Wtbeam, typename Handler>
	void ParseString(Stream& stream, Handler& handler) {
#define Z16 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
		static const Ch escape[256] = {
			Z16, Z16, 0, 0,'\"', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,'/', 
			Z16, Z16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,'\\', 0, 0, 0, 
			0, 0,'\b', 0, 0, 0,'\f', 0, 0, 0, 0, 0, 0, 0,'\n', 0, 
			0, 0,'\r', 0,'\t', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			Z16, Z16, Z16, Z16, Z16, Z16, Z16, Z16
		};
#undef Z16

		Stream s = stream;	// Use a local copy for optimization
		RAPIDJSON_ASSERT(s.Peek() == '\"');
		s.Take();	// Skip '\"'
		Ch *head;
		SizeType len;
		if (parseFlags & kParseInsituFlag)
			head = s.PutBegin();
		else
			len = 0;

#define RAPIDJSON_PUT(x) \
	do { \
		if (parseFlags & kParseInsituFlag) \
			s.Put(x); \
		else { \
			*stack_.template Push<Ch>() = x; \
			++len; \
		} \
	} while(false)

		for (;;) {
			Ch c = s.Take();
			if (c == '\\') {	// Escape
				Ch e = s.Take();
				if ((sizeof(Ch) == 1 || (e & 0xFFFFFF00) == 0) && escape[(unsigned char)e])
					RAPIDJSON_PUT(escape[(unsigned char)e]);
				else if (e == 'u') {	// Unicode
					unsigned codepoint = ParseHex4(s);
					if (codepoint >= 0xD800 && codepoint <= 0xDBFF) { // Handle UTF-16 surrogate pair
						if (s.Take() != '\\' || s.Take() != 'u') {
							RAPIDJSON_PARSE_ERROR("Missing the second \\u in surrogate pair", s.Tell() - 2);
							return;
						}
						unsigned codepoint2 = ParseHex4(s);
						if (codepoint2 < 0xDC00 || codepoint2 > 0xDFFF) {
							RAPIDJSON_PARSE_ERROR("The second \\u in surrogate pair is invalid", s.Tell() - 2);
							return;
						}
						codepoint = (((codepoint - 0xD800) << 10) | (codepoint2 - 0xDC00)) + 0x10000;
					}

					Ch buffer[4];
					SizeType count = SizeType(Encoding::Encode(buffer, codepoint) - &buffer[0]);

					if (parseFlags & kParseInsituFlag) 
						for (SizeType i = 0; i < count; i++)
							s.Put(buffer[i]);
					else {
						memcpy(stack_.template Push<Ch>(count), buffer, count * sizeof(Ch));
						len += count;
					}
				}
				else {
					RAPIDJSON_PARSE_ERROR("Unknown escape character", stream.Tell() - 1);
					return;
				}
			}
			else if (c == '"') {	// Closing double quote
				if (parseFlags & kParseInsituFlag) {
					size_t length = s.PutEnd(head);
					RAPIDJSON_ASSERT(length <= 0xFFFFFFFF);
					RAPIDJSON_PUT('\0');	// null-terminate the string
					handler.String(head, SizeType(length), false);
				}
				else {
					RAPIDJSON_PUT('\0');
					handler.String(stack_.template Pop<Ch>(len), len - 1, true);
				}
				stream = s;	// restore stream
				return;
			}
			else if (c == '\0') {
				RAPIDJSON_PARSE_ERROR("lacks ending quotation before the end of string", stream.Tell() - 1);
				return;
			}
			else if ((unsigned)c < 0x20) {	// RFC 4627: unescaped = %x20-21 / %x23-5B / %x5D-10FFFF
				RAPIDJSON_PARSE_ERROR("Incorrect unescaped character in string", stream.Tell() - 1);
				return;
			}
			else
				RAPIDJSON_PUT(c);	// Normal character, just copy
		}
#undef RAPIDJSON_PUT
	}

	template<unsigned parseFlags, typename Stream, typename Handler>
	void ParseNumber(Stream& stream, Handler& handler) {
		Stream s = stream; // Local copy for optimization
		// Parse minus
		bool minus = false;
		if (s.Peek() == '-') {
			minus = true;
			s.Take();
		}

		// Parse int: zero / ( digit1-9 *DIGIT )
		unsigned i;
		bool try64bit = false;
		if (s.Peek() == '0') {
			i = 0;
			s.Take();
		}
		else if (s.Peek() >= '1' && s.Peek() <= '9') {
			i = s.Take() - '0';

			if (minus)
				while (s.Peek() >= '0' && s.Peek() <= '9') {
					if (i >= 214748364) { // 2^31 = 2147483648
						if (i != 214748364 || s.Peek() > '8') {
							try64bit = true;
							break;
						}
					}
					i = i * 10 + (s.Take() - '0');
				}
			else
				while (s.Peek() >= '0' && s.Peek() <= '9') {
					if (i >= 429496729) { // 2^32 - 1 = 4294967295
						if (i != 429496729 || s.Peek() > '5') {
							try64bit = true;
							break;
						}
					}
					i = i * 10 + (s.Take() - '0');
				}
		}
		else {
			RAPIDJSON_PARSE_ERROR("Expect a value here.", stream.Tell());
			return;
		}

		// Parse 64bit int
		uint64_t i64 = 0;
		bool useDouble = false;
		if (try64bit) {
			i64 = i;
			if (minus) 
				while (s.Peek() >= '0' && s.Peek() <= '9') {					
					if (i64 >= 922337203685477580uLL) // 2^63 = 9223372036854775808
						if (i64 != 922337203685477580uLL || s.Peek() > '8') {
							useDouble = true;
							break;
						}
					i64 = i64 * 10 + (s.Take() - '0');
				}
			else
				while (s.Peek() >= '0' && s.Peek() <= '9') {					
					if (i64 >= 1844674407370955161uLL) // 2^64 - 1 = 18446744073709551615
						if (i64 != 1844674407370955161uLL || s.Peek() > '5') {
							useDouble = true;
							break;
						}
					i64 = i64 * 10 + (s.Take() - '0');
				}
		}

		// Force double for big integer
		double d = 0.0;
		if (useDouble) {
			d = (double)i64;
			while (s.Peek() >= '0' && s.Peek() <= '9') {
				if (d >= 1E307) {
					RAPIDJSON_PARSE_ERROR("Number too big to store in double", stream.Tell());
					return;
				}
				d = d * 10 + (s.Take() - '0');
			}
		}

		// Parse frac = decimal-point 1*DIGIT
		int expFrac = 0;
		if (s.Peek() == '.') {
			if (!useDouble) {
				d = try64bit ? (double)i64 : (double)i;
				useDouble = true;
			}
			s.Take();

			if (s.Peek() >= '0' && s.Peek() <= '9') {
				d = d * 10 + (s.Take() - '0');
				--expFrac;
			}
			else {
				RAPIDJSON_PARSE_ERROR("At least one digit in fraction part", stream.Tell());
				return;
			}

			while (s.Peek() >= '0' && s.Peek() <= '9') {
				if (expFrac > -16) {
					d = d * 10 + (s.Peek() - '0');
					--expFrac;
				}
				s.Take();
			}
		}

		// Parse exp = e [ minus / plus ] 1*DIGIT
		int exp = 0;
		if (s.Peek() == 'e' || s.Peek() == 'E') {
			if (!useDouble) {
				d = try64bit ? (double)i64 : (double)i;
				useDouble = true;
			}
			s.Take();

			bool expMinus = false;
			if (s.Peek() == '+')
				s.Take();
			else if (s.Peek() == '-') {
				s.Take();
				expMinus = true;
			}

			if (s.Peek() >= '0' && s.Peek() <= '9') {
				exp = s.Take() - '0';
				while (s.Peek() >= '0' && s.Peek() <= '9') {
					exp = exp * 10 + (s.Take() - '0');
					if (exp > 308) {
						RAPIDJSON_PARSE_ERROR("Number too big to store in double", stream.Tell());
						return;
					}
				}
			}
			else {
				RAPIDJSON_PARSE_ERROR("At least one digit in exponent", s.Tell());
				return;
			}

			if (expMinus)
				exp = -exp;
		}

		// Finish parsing, call event according to the type of number.
		if (useDouble) {
			d *= internal::Pow10(exp + expFrac);
			handler.Double(minus ? -d : d);
		}
		else {
			if (try64bit) {
				if (minus)
					handler.Int64(-(int64_t)i64);
				else
					handler.Uint64(i64);
			}
			else {
				if (minus)
					handler.Int(-(int)i);
				else
					handler.Uint(i);
			}
		}

		stream = s; // restore stream
	}

	// Parse any JSON value
	template<unsigned parseFlags, typename Stream, typename Handler>
	void ParseValue(Stream& stream, Handler& handler) {
		switch (stream.Peek()) {
			case 'n': ParseNull  <parseFlags>(stream, handler); break;
			case 't': ParseTrue  <parseFlags>(stream, handler); break;
			case 'f': ParseFalse <parseFlags>(stream, handler); break;
			case '"': ParseString<parseFlags>(stream, handler); break;
			case '{': ParseObject<parseFlags>(stream, handler); break;
			case '[': ParseArray <parseFlags>(stream, handler); break;
			default : ParseNumber<parseFlags>(stream, handler);
		}
	}

	static const size_t kDefaultStackCapacity = 256;	//!< Default stack capacity in bytes for storing a single decoded string. 
	internal::Stack<Allocator> stack_;	//!< A stack for storing decoded string temporarily during non-destructive parsing.
	jmp_buf jmpbuf_;					//!< setjmp buffer for fast exit from nested parsing function calls.
	const char* parseError_;
	size_t errorOffset_;
}; // class GenericReader

//! Reader with UTF8 encoding and default allocator.
typedef GenericReader<UTF8<> > Reader;

} // namespace rapidjson

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // RAPIDJSON_READER_H_
