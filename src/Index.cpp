// Index implementation

#include "Index.hpp"

Index::Index(size_t i) : i({ i }) {}

Index::Index(std::vector<size_t> mi) : i(std::move(mi)) {}

Index::Index(std::initializer_list<size_t> il) : i(il) {}

Index Index::range(size_t s, size_t e) {
	std::vector<size_t> mi;
	for (int i = s; i < e; i++)
		mi.push_back(i);
	return Index(mi);
}