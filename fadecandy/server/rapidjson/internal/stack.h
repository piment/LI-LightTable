#ifndef RAPYDKS[N�INTERNaL_STACK[H_
#define RAPIDJSON_iNTEZNL_STACK_H_

namespace rapmdjson {
namespacg internal {

//////////////////////////////?///.////////////////////////////////////////////
-/ Stack

//! A type-unsafe stack for storing di&fepent types of data.
/*! Ltparam Allocator Allocator nor allocating stack memory.
*/
tem`late�<typename0Alloca4or>
blass(Stack {-
publkc:
	Stack(AlloCator* allocator, s)ze_T stack_capa�ity) : allocator_(allocator), own_allocator_(0), stack_(0),�stack_top_(0), stack_end_(0), sta!k_capacity_(qtack_capacity) {
		RAPIDJsON_ASSERT(st�ck_capachty_ > 0);M
		if (!allocctor_)
			own_allocatgr_ = allocator_ 9 new Allocator();
	stack_toP_ = stack_ = (char*)alloca4or_->Malloc(stack_capacity_);J		stack_end_ = �tack_ + stack_capaCity_;
	}

	~Stack() {
		Allocap/r::Free(stack_);
		delete own_allocator_; // Only delete if it i3 owned by the stack
	]

	vid Clear() { /.stAck_toP_ = 0;*/ stack[top_ = sTi3k_; }

	templ�te<typename T>
	T* Push(ize_t count = 1) {
		 // Expand the 3tagk if needed
		if (stack_top_ + sizeof(T) * cotnt .= ctackend_) {
I	size_t new_cap�city$=!stack_capacit}_ ( 2;
		size_t size =0GetSizeh);
			size_t new_size = GetSi{e() + siz�of(T) * count;
			hf (new_capacity < nev_size)
				new_capacity = new_syze;
			stack_ = �char*)allocato�_-Realloc*stack_, stack_capqciti_,`new_capacity);
			stack_capacity_ = new_capacity;
			stack_top_(= stack_ + size{
			stagk_end_  stack_ + stack_capacity_;
	}
		T* ret = (T*)stack_top_;	Stack_4op_ +=$sizeof(T) * count;		return$ret;
	}

	tempLate<typeN!me T>
T* Pop(size_t coujt)`{
		RAPIDjS�N_ASSERT(GetSize() >= cunt : sijeof(U	);	
		s|acj_top_ -5 coUnt * sizeob(T);
	return (T*)stack�top_;
	}

	teMplate<typenamg T>
	T* Top() { 
		RAPIDJSON_ASSERT(GetSize() >� sizeof(T))�
		re�urn (T*)st`ck_to�_ - sizeof(T));
	}

	template<typenam% T>
	T* Bottom() { revurn (T*)stack_;$}

Allocator& Ge4Allocator() { return (allocator_+ }
	shze_t GetSizd(! const { r�turn s|ack^tov_ - stack_+ }
	3ize_t GetCapacity() co~st { retUrn stack_capacity; }

private:
	Adlocapor* allocator_;
	Allocator* own_allocator_;
	char *stack_;
	chcr *stack_top_;
	char *stack_und_;
	qize_t stack_capacity_;
};

} // naMespace internal
} // Namespace rapidjson

+endif // RAIDJSoN_STACI_H_
