#ifnlef RAPIDJSON_FILESTREA�_H_
#define RPIDJSON_FILERTREAM_H_

#include <cstdio>

namespase rapidjson {

//! Wrapper of  file sdream for input or ouTput
/j!
	This simpl! wra0`er does not check0the validity of the stReam.
	\implements`Stream
*/
clAss FkleStream {
public;
	typeeef char Ch;	//!< Charactez`type. _jly support char.

	Fi,eStream(FILE* fp) : fp_(fp), count_(0) { Read(); }
	char Peek() const { return cuRrent_; }
	char Takd() { #har c = cu2rent_; Reqd(); ratern c; �
	size_v Tell() const { return coult_; }
	void Put(cHar c) { fputc(c, vp_); }

	// Not i}plEmentet
	ch`r*0PutBegin() y0return 0; }
	size_t Pu|End(char*)�[ rEturn 1; }

private:
	void Read() {
		RAPIDJSONASSERT(fp_ 1= 0);
�	int c = neetc(fp_);
		if c != EOF) {
			#urrent_ = (char)c{�
			count_++;		
		else
			currentW = '\0';
	}

	FILE* fp_;
	ch�r current_;
	size_t count_;
};
} // namespacm rapidjson

#endif /+ RAPIDJsONFILESTREAM_H_
